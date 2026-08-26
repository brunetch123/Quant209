"""Behavioural tests for the PySide6 app shell: stepper gating, navigation, theme,
state signals, worker threading, and the cold-start import discipline.

Why ``QT_QPA_PLATFORM=offscreen`` is set here, at module scope
--------------------------------------------------------------
It has to be in the environment before the *first* ``QApplication`` is constructed, not
before PySide6 is imported. pytest-qt builds that application lazily in its session-scoped
``qapp`` fixture, so setting it at import time of this module (collection, long before any
test body runs) is early enough and keeps the requirement next to the tests that need it.
This is the only Qt-touching module in the suite. ``setdefault`` so a developer running
with a real display (``QT_QPA_PLATFORM=windows``) to *look* at the window keeps their choice.

The QSettings discipline
-------------------------
:class:`~quant209.app.state.AppState` persists the msconvert path. Every test here injects a
throwaway INI-file ``QSettings`` under ``tmp_path`` so the suite never reads or writes the
developer's real registry hive. A test that silently changed the app's remembered msconvert
path would be difficult to debug.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import subprocess
import sys
import threading
from pathlib import Path

import pytest

from PySide6.QtCore import QSettings, Qt, QThreadPool

from quant209.app import theme
from quant209.app.main import MainWindow, configure_app
from quant209.app.state import (
    STEP_ASSIGN,
    STEP_CONVERT,
    STEP_EXPORT,
    STEP_LABELS,
    STEP_PROCESS,
    STEP_REVIEW,
    STEP_SELECT,
    AppState,
)
from quant209.app.workers import ConvertWorker, PipelineWorker
from quant209.core.sequence import FileEntry


# --------------------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------------------


@pytest.fixture
def settings(tmp_path: Path) -> QSettings:
    """A throwaway INI-backed QSettings, so no test touches the real registry."""
    return QSettings(str(tmp_path / "quant209.ini"), QSettings.Format.IniFormat)


@pytest.fixture
def state(settings: QSettings) -> AppState:
    return AppState(settings=settings)


@pytest.fixture
def window(qtbot, qapp, state: AppState) -> MainWindow:
    configure_app(qapp)
    win = MainWindow(state)
    qtbot.addWidget(win)
    return win


def _entry(tmp_path: Path, name: str, file_type: str = "sample", *, create: bool = True) -> FileEntry:
    path = tmp_path / name
    if create:
        if path.suffix.lower() == ".mzml":
            path.write_text("<mzML/>", encoding="utf-8")
        else:
            path.mkdir(exist_ok=True)
    return FileEntry(path=path, file_type=file_type)


def _mzml_run(tmp_path: Path) -> list[FileEntry]:
    """A run whose entries are already mzML: nothing left to convert."""
    return [
        _entry(tmp_path, "AN_Cal01_20231002.mzML", "cal"),
        _entry(tmp_path, "AN_SMP16_20231002.mzML", "sample"),
    ]


def _raw_run(tmp_path: Path) -> list[FileEntry]:
    """A run of raw Agilent ``.D`` acquisitions: conversion still owed."""
    return [
        _entry(tmp_path, "AN_Cal01_20231002.D", "cal"),
        _entry(tmp_path, "AN_SMP16_20231002.D", "sample"),
    ]


# --------------------------------------------------------------------------------------
# Window chrome
# --------------------------------------------------------------------------------------


def test_window_shows_the_six_workflow_steps_in_order(window: MainWindow):
    assert [row.step_label for row in window.stepper.steps] == list(STEP_LABELS)
    assert [row.step_label for row in window.stepper.steps] == [
        "Select", "Assign", "Convert", "Process", "Review", "Export"
    ]


def test_window_opens_at_the_specified_minimum_size(window: MainWindow):
    assert window.minimumWidth() >= 1200
    assert window.minimumHeight() >= 760
    assert window.windowTitle() == "Quant 209"


def test_a_fresh_run_unlocks_select_only(window: MainWindow):
    assert window.state.furthest_step() == STEP_SELECT
    assert window.current_step == STEP_SELECT
    assert window.stepper.state_of(STEP_SELECT) == "active"
    for step in (STEP_ASSIGN, STEP_CONVERT, STEP_PROCESS, STEP_REVIEW, STEP_EXPORT):
        assert window.stepper.state_of(step) == "locked", STEP_LABELS[step]
        assert not window.stepper.steps[step].isEnabled()


def test_entries_unlock_assign_and_convert(window: MainWindow, tmp_path: Path):
    window.state.set_entries(_raw_run(tmp_path))

    assert window.state.furthest_step() == STEP_CONVERT
    assert window.state.is_step_available(STEP_ASSIGN)
    assert window.state.is_step_available(STEP_CONVERT)
    # Assign is *complete* the moment files are selected. Roles are auto-detected from
    # the filenames, so its output already exists. Convert still owes work.
    assert window.stepper.state_of(STEP_ASSIGN) == "complete"
    assert window.stepper.state_of(STEP_CONVERT) == "available"
    assert window.stepper.state_of(STEP_PROCESS) == "locked"
    assert window.stepper.state_of(STEP_SELECT) == "active"  # active outranks complete


def test_entries_that_are_already_mzml_unlock_process(window: MainWindow, tmp_path: Path):
    window.state.set_entries(_mzml_run(tmp_path))

    assert window.state.mzml_ready()
    assert window.state.furthest_step() == STEP_PROCESS
    assert window.stepper.state_of(STEP_PROCESS) == "available"
    assert window.stepper.state_of(STEP_REVIEW) == "locked"


def test_converted_output_unlocks_process(window: MainWindow, tmp_path: Path):
    window.state.set_run_folder(tmp_path)
    window.state.set_entries(_raw_run(tmp_path))
    assert window.state.furthest_step() == STEP_CONVERT

    converted = window.state.converted_dir
    converted.mkdir(parents=True, exist_ok=True)
    for entry in window.state.entries:
        (converted / f"{entry.path.stem}.mzML").write_text("<mzML/>", encoding="utf-8")
    window.state.refresh()

    assert window.state.furthest_step() == STEP_PROCESS
    assert window.stepper.state_of(STEP_CONVERT) == "complete"


def test_a_result_unlocks_review_and_export(window: MainWindow, tmp_path: Path):
    window.state.set_entries(_mzml_run(tmp_path))
    window.state.set_result(object())

    assert window.state.furthest_step() == STEP_EXPORT
    assert window.stepper.state_of(STEP_REVIEW) == "available"
    assert window.stepper.state_of(STEP_EXPORT) == "available"
    assert window.stepper.state_of(STEP_PROCESS) == "complete"


def test_review_is_only_complete_once_it_has_been_opened(window: MainWindow, tmp_path: Path):
    """A check mark against a screen nobody looked at would be a claim the app cannot make."""
    window.state.set_entries(_mzml_run(tmp_path))
    window.state.set_result(object())
    assert STEP_REVIEW not in window.state.completed_steps()

    window.goto_step(STEP_REVIEW)
    window.goto_step(STEP_EXPORT)

    assert window.stepper.state_of(STEP_REVIEW) == "complete"


def test_a_new_result_un_reviews_the_run(window: MainWindow, tmp_path: Path):
    window.state.set_entries(_mzml_run(tmp_path))
    window.state.set_result(object())
    window.goto_step(STEP_REVIEW)
    assert STEP_REVIEW in window.state.completed_steps()

    window.state.set_result(object())

    assert STEP_REVIEW not in window.state.completed_steps()


# --------------------------------------------------------------------------------------
# Navigation
# --------------------------------------------------------------------------------------


def test_clicking_an_unlocked_step_navigates(window: MainWindow, qtbot, tmp_path: Path):
    window.state.set_entries(_raw_run(tmp_path))

    qtbot.mouseClick(window.stepper.steps[STEP_CONVERT], Qt.MouseButton.LeftButton)

    assert window.current_step == STEP_CONVERT
    assert window.stack.currentWidget() is window.page(STEP_CONVERT)
    assert window.stepper.state_of(STEP_CONVERT) == "active"


def test_clicking_a_locked_step_does_nothing(window: MainWindow, qtbot):
    qtbot.mouseClick(window.stepper.steps[STEP_REVIEW], Qt.MouseButton.LeftButton)

    assert window.current_step == STEP_SELECT
    assert window.stack.currentWidget() is window.page(STEP_SELECT)


def test_back_navigation_is_always_allowed(window: MainWindow, qtbot, tmp_path: Path):
    window.state.set_entries(_mzml_run(tmp_path))
    window.state.set_result(object())
    assert window.goto_step(STEP_EXPORT)

    qtbot.mouseClick(window.stepper.steps[STEP_SELECT], Qt.MouseButton.LeftButton)

    assert window.current_step == STEP_SELECT


def test_focus_follows_navigation_only_when_it_is_already_in_the_stepper(
    window: MainWindow, qtbot, tmp_path: Path
):
    """Two rows must never look highlighted at once. See ``Stepper.follow_focus``."""
    window.show()
    qtbot.waitExposed(window)
    window.state.set_entries(_mzml_run(tmp_path))

    window.stepper.steps[STEP_SELECT].setFocus()
    window.goto_step(STEP_CONVERT)
    assert window.focusWidget() is window.stepper.steps[STEP_CONVERT]

    # Focus parked on the page instead: navigation must leave it alone. It has to be a
    # widget that is actually on screen. ``setFocus`` on a hidden one is a no-op, and the
    # test would then be asserting about focus that never moved off the stepper.
    page_widget = window.page(STEP_CONVERT).continue_button
    assert page_widget.isVisible()
    page_widget.setFocus()
    window.goto_step(STEP_ASSIGN)
    assert window.focusWidget() is not window.stepper.steps[STEP_ASSIGN]


def test_goto_step_refuses_a_locked_step(window: MainWindow):
    assert window.goto_step(STEP_PROCESS) is False
    assert window.current_step == STEP_SELECT


def test_the_header_follows_the_active_page(window: MainWindow, tmp_path: Path):
    window.state.set_entries(_raw_run(tmp_path))
    select_title = window.header.title_label.text()
    window.goto_step(STEP_ASSIGN)

    assert window.header.title_label.text() == window.page(STEP_ASSIGN).title
    assert window.header.title_label.text() != select_title


def test_pages_are_built_lazily(window: MainWindow, tmp_path: Path):
    """Only the Select page exists before anything is navigated to. The whole
    cold-start budget rests on Review/Export (pyqtgraph, pandas) staying unimported."""
    assert window.built_steps() == [STEP_SELECT]

    window.state.set_entries(_mzml_run(tmp_path))
    window.goto_step(STEP_CONVERT)

    assert window.built_steps() == [STEP_SELECT, STEP_CONVERT]


def test_the_sidebar_names_the_run_folder_once_chosen(window: MainWindow, tmp_path: Path):
    run_folder = tmp_path / "SMP16-30_20231002"
    run_folder.mkdir()
    assert not window.sidebar.run_folder_label.isVisible() or not window.sidebar.run_folder_label.text()

    window.state.set_run_folder(run_folder)

    assert run_folder.name in window.sidebar.run_folder_label.text()
    assert str(run_folder) in window.sidebar.run_folder_label.toolTip()


# --------------------------------------------------------------------------------------
# Theme
# --------------------------------------------------------------------------------------


def test_theme_is_applied_to_the_application(qapp):
    configure_app(qapp)
    qss = qapp.styleSheet()

    assert qss.strip()
    assert theme.PALETTE["accent"] in qss
    assert theme.PALETTE["accent"] == "#0F6E6B"


def test_palette_exports_the_design_system_colors():
    for key in (
        "background", "surface", "border", "text", "text_secondary",
        "accent", "accent_hover", "accent_pressed", "accent_contrast",
        "success", "warning", "error",
    ):
        assert theme.PALETTE[key].startswith("#"), key
    assert theme.PALETTE["background"] == "#FAFBFC"
    assert theme.PALETTE["surface"] == "#FFFFFF"
    assert theme.PALETTE["border"] == "#E3E6E8"
    assert theme.PALETTE["text"] == "#1A1D1F"
    assert theme.PALETTE["text_secondary"] == "#6A7178"


def test_no_default_qt_gray_dialog_backgrounds(qapp):
    """A spot-check that the stock Qt gray never shows through: every surface the theme
    paints is named in the stylesheet."""
    configure_app(qapp)
    qss = qapp.styleSheet()
    for token in ("QScrollBar", "QToolTip", "QTableView", "QHeaderView", "QProgressBar", "QPushButton"):
        assert token in qss, token


# --------------------------------------------------------------------------------------
# Cold start: import discipline
# --------------------------------------------------------------------------------------


def test_importing_the_app_entry_point_pulls_no_scientific_stack():
    """``quant209.app.main`` must reach a painted window without numpy/pandas/pyqtgraph.

    Checked in a subprocess: this test session has already imported pandas via the core
    tests, so an in-process ``sys.modules`` assertion would be meaningless.
    """
    probe = (
        "import sys, quant209.app.main;"
        "banned=('numpy','pandas','pyqtgraph','matplotlib','openpyxl',"
        "'quant209.core.pipeline','quant209.report.excel');"
        "print(','.join(m for m in banned if m in sys.modules))"
    )
    proc = subprocess.run([sys.executable, "-c", probe], capture_output=True, text=True, timeout=120)

    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "", f"eagerly imported: {proc.stdout.strip()}"


# --------------------------------------------------------------------------------------
# State
# --------------------------------------------------------------------------------------


def test_state_emits_on_every_transition(state: AppState, qtbot, tmp_path: Path):
    with qtbot.waitSignal(state.entries_changed, timeout=1000):
        state.set_entries(_mzml_run(tmp_path))
    with qtbot.waitSignal(state.result_changed, timeout=1000):
        state.set_result(object())
    with qtbot.waitSignal(state.busy_changed, timeout=1000):
        state.set_busy(True)
    assert state.busy is True


def test_choosing_a_new_run_folder_clears_the_previous_run(state: AppState, tmp_path: Path):
    first, second = tmp_path / "run_a", tmp_path / "run_b"
    first.mkdir()
    second.mkdir()
    state.set_run_folder(first)
    state.set_entries(_mzml_run(tmp_path))
    state.set_result(object())

    state.set_run_folder(second)

    assert state.entries == []
    assert state.result is None
    assert state.converted_dir == second / "Converted"
    assert state.furthest_step() == STEP_SELECT


def test_msconvert_path_is_autodetected_from_the_run_folder(state: AppState, tmp_path: Path):
    exe = tmp_path / "MSconvert" / "pwiz-bin" / "msconvert.exe"
    exe.parent.mkdir(parents=True)
    exe.write_text("", encoding="utf-8")
    run_folder = tmp_path / "SMP16-30_20231002"
    run_folder.mkdir()

    state.set_run_folder(run_folder)

    assert state.msconvert_path == exe


def test_msconvert_path_persists_between_sessions(settings: QSettings, tmp_path: Path):
    exe = tmp_path / "msconvert.exe"
    exe.write_text("", encoding="utf-8")
    first = AppState(settings=settings)
    first.set_msconvert_path(exe)

    assert AppState(settings=settings).msconvert_path == exe


# --------------------------------------------------------------------------------------
# Workers
# --------------------------------------------------------------------------------------


@pytest.fixture
def pool() -> QThreadPool:
    return QThreadPool.globalInstance()


def test_pipeline_worker_reports_progress_then_the_result(monkeypatch, qtbot, pool, tmp_path: Path):
    sentinel = object()

    def stub(entries, converted_dir, params=None, progress=None, **kwargs):
        progress("calibrate", 0, 3, "read calibration injection")
        progress("integrate", 0, 2, "SMP16")
        return sentinel

    monkeypatch.setattr("quant209.core.pipeline.run_analysis", stub)
    worker = PipelineWorker(_mzml_run(tmp_path), tmp_path / "Converted")
    seen: list[tuple] = []
    threads: set[int] = set()
    worker.signals.progress.connect(
        lambda *a: (seen.append(a), threads.add(threading.get_ident()))
    )
    finished: list[object] = []
    worker.signals.finished.connect(finished.append)

    with qtbot.waitSignal(worker.signals.finished, timeout=10_000):
        pool.start(worker)

    assert [a[0] for a in seen] == ["calibrate", "integrate"]
    assert seen[1] == ("integrate", 0, 2, "SMP16")
    assert finished == [sentinel]
    assert threads == {threading.get_ident()}, "signals must be delivered on the UI thread"


def test_pipeline_worker_turns_an_exception_into_an_error_signal(monkeypatch, qtbot, pool, tmp_path: Path):
    def stub(*args, **kwargs):
        raise ValueError("this run has no calibration injection")

    monkeypatch.setattr("quant209.core.pipeline.run_analysis", stub)
    worker = PipelineWorker(_mzml_run(tmp_path), tmp_path / "Converted")
    finished: list[object] = []
    worker.signals.finished.connect(finished.append)

    with qtbot.waitSignal(worker.signals.error, timeout=10_000) as blocker:
        pool.start(worker)

    message, trace = blocker.args
    assert "no calibration injection" in message
    assert "ValueError" in trace
    assert finished == []


def test_pipeline_worker_cancellation_abandons_the_run(monkeypatch, qtbot, pool, tmp_path: Path):
    gate = threading.Event()
    reached: list[str] = []

    def stub(entries, converted_dir, params=None, progress=None, **kwargs):
        progress("integrate", 0, 2, "first")
        reached.append("first")
        gate.wait(10)
        progress("integrate", 1, 2, "second")  # must not return from here
        reached.append("second")
        return object()

    monkeypatch.setattr("quant209.core.pipeline.run_analysis", stub)
    worker = PipelineWorker(_mzml_run(tmp_path), tmp_path / "Converted")
    finished: list[object] = []
    errors: list[tuple] = []
    worker.signals.finished.connect(finished.append)
    worker.signals.error.connect(lambda *a: errors.append(a))

    def on_progress(*_args):
        worker.cancel()
        gate.set()

    worker.signals.progress.connect(on_progress)

    with qtbot.waitSignal(worker.signals.cancelled, timeout=10_000):
        pool.start(worker)

    assert reached == ["first"], "the pipeline must abort at the next checkpoint"
    assert finished == []
    assert errors == []


def test_convert_worker_wraps_convert_many(monkeypatch, qtbot, pool, tmp_path: Path):
    results = ["converted-a", "converted-b"]

    def stub(sources, msconvert_exe, out_dir, extra_args=(), progress=None):
        for i, src in enumerate(sources):
            progress(i, len(sources), Path(src).name)
        return results

    monkeypatch.setattr("quant209.convert.msconvert.convert_many", stub)
    sources = [e.path for e in _raw_run(tmp_path)]
    worker = ConvertWorker(sources, tmp_path / "msconvert.exe", tmp_path / "Converted")
    seen: list[tuple] = []
    worker.signals.progress.connect(lambda *a: seen.append(a))

    with qtbot.waitSignal(worker.signals.finished, timeout=10_000) as blocker:
        pool.start(worker)

    assert blocker.args[0] == results
    assert [a[0] for a in seen] == ["convert", "convert"]
    assert seen[0][3] == sources[0].name
    assert seen[0][2] == 2
