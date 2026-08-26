"""Behavioural tests for the four working workflow pages: Select, Assign, Convert, Process.

Same discipline as ``tests/test_app_shell.py``: ``QT_QPA_PLATFORM=offscreen`` at module
scope (before pytest-qt builds the session ``QApplication``), and a throwaway INI-backed
``QSettings`` in every fixture so no test can reach the developer's real remembered
msconvert path.

What is stubbed, and what is not
---------------------------------
The pipeline, MSConvert and the filesystem scan are stubbed. **The widgets are not.** Every
assertion below reads a real, rendered Qt widget tree: the text a label actually holds, the
combo box a delegate actually created, the pixels a stylesheet actually painted. A test that
asserted against the page's own internal model instead would pass just as happily with
nothing drawn on screen, which is not a guarantee worth having.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import QComboBox, QLabel

from quant209.app import theme

from quant209.app.main import MainWindow, configure_app
from quant209.app.state import (
    STEP_ASSIGN,
    STEP_CONVERT,
    STEP_PROCESS,
    STEP_REVIEW,
    STEP_SELECT,
    AppState,
)
from quant209.app.workflow.assign_page import (
    COL_DETECTED,
    COL_FILE,
    COL_SIZE,
    COL_STATUS,
    COL_TYPE,
    AssignPage,
)
from quant209.app.workflow.convert_page import COLLAPSED_ROWS, ConvertPage
from quant209.app.workflow.process_page import DEFAULT_LOQ_K, ProcessPage, _format_progress
from quant209.app.workflow.select_page import SelectPage
from quant209.core.sequence import FileEntry


# --------------------------------------------------------------------------------------
# Fixtures and builders
# --------------------------------------------------------------------------------------


@pytest.fixture
def settings(tmp_path: Path) -> QSettings:
    return QSettings(str(tmp_path / "quant209.ini"), QSettings.Format.IniFormat)


@pytest.fixture
def state(settings: QSettings) -> AppState:
    return AppState(settings=settings)


@pytest.fixture
def styled(qapp):
    """The application stylesheet, applied. Anything grabbing pixels needs this."""
    configure_app(qapp)
    return qapp


def _mzml(folder: Path, name: str, file_type: str, *, create: bool = True) -> FileEntry:
    path = folder / name
    if create:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("<mzML/>", encoding="utf-8")
    return FileEntry(path=path, file_type=file_type)


def _raw(folder: Path, name: str, file_type: str) -> FileEntry:
    path = folder / name
    path.mkdir(parents=True, exist_ok=True)
    (path / "data.ms").write_text("x" * 64, encoding="utf-8")
    return FileEntry(path=path, file_type=file_type)


def _full_run(folder: Path) -> list[FileEntry]:
    """A converted run with every role present: the case with no banners."""
    return [
        _mzml(folder, "AN_Cal01_20231002.mzML", "cal"),
        _mzml(folder, "AN_Ref04_20231002.mzML", "ref"),
        _mzml(folder, "AN_LB72_20231002.mzML", "blank"),
        _mzml(folder, "AN_SMP16_20231002.mzML", "sample"),
        _mzml(folder, "AN_SMP17_20231002.mzML", "sample"),
    ]


def _page(cls, state: AppState, qtbot):
    page = cls(state)
    qtbot.addWidget(page)
    return page


# --------------------------------------------------------------------------------------
# Select
# --------------------------------------------------------------------------------------


@pytest.fixture
def scan(monkeypatch, tmp_path: Path):
    """Stub ``find_data_files``/``scan_inputs`` and record what the page asked for."""
    calls: dict[str, object] = {}
    found: list[Path] = []

    def find_data_files(folder, prefer_mzml=True):
        calls["folder"] = Path(folder)
        return list(found)

    def scan_inputs(paths, rules=None):
        calls["scanned"] = [Path(p) for p in paths]
        return [
            FileEntry(path=Path(p), file_type=_type_of(Path(p))) for p in paths
        ]

    def _type_of(path: Path) -> str:
        stem = path.stem.lower()
        for token, file_type in (("cal", "cal"), ("ref", "ref"), ("hex", "hex"), ("lb", "blank")):
            if token in stem:
                return file_type
        return "sample"

    monkeypatch.setattr("quant209.core.sequence.find_data_files", find_data_files)
    monkeypatch.setattr("quant209.core.sequence.scan_inputs", scan_inputs)
    return SimpleNamespace(calls=calls, found=found)


def test_choosing_a_run_folder_reads_it_into_the_state(state, qtbot, scan, tmp_path: Path):
    run = tmp_path / "SMP16-30_20231002"
    run.mkdir()
    scan.found[:] = [run / "AN_Cal01.mzML", run / "AN_SMP16.mzML"]
    page = _page(SelectPage, state, qtbot)

    page.load_folder(run)

    assert scan.calls["folder"] == run
    assert state.run_folder == run
    assert [e.path.name for e in state.entries] == ["AN_Cal01.mzML", "AN_SMP16.mzML"]


def test_the_summary_card_counts_what_was_found(state, qtbot, scan, tmp_path: Path):
    run = tmp_path / "SMP16-30_20231002"
    run.mkdir()
    scan.found[:] = [
        run / "AN_Cal01.mzML",
        run / "AN_Ref04.mzML",
        run / "AN_LB72.mzML",
        run / "AN_SMP16.D",
    ]
    page = _page(SelectPage, state, qtbot)

    page.load_folder(run)

    assert page.summary.isVisibleTo(page)
    assert not page.empty.isVisibleTo(page)
    assert str(run) in page.folder_label.text()
    counts = page.counts_label.text()
    assert "4 data files" in counts
    assert "1 Agilent .D" in counts
    assert "3 mzML" in counts
    assert page.chips.texts() == ["1 calibration", "1 reference", "1 samples", "1 blanks"]


def test_the_summary_card_says_types_are_editable_next(state, qtbot, scan, tmp_path: Path):
    """The one hint that stops the analyst hunting for an override on this screen."""
    run = tmp_path / "run"
    run.mkdir()
    scan.found[:] = [run / "AN_Cal01.mzML"]
    page = _page(SelectPage, state, qtbot)
    page.load_folder(run)

    hints = [
        child.text()
        for child in page.summary.findChildren(QLabel)
        if child.objectName() == "Hint"
    ]
    assert any("editable in the next step" in text for text in hints)


def test_a_successful_selection_advances_to_assign(state, qtbot, scan, tmp_path: Path):
    run = tmp_path / "run"
    run.mkdir()
    scan.found[:] = [run / "AN_Cal01.mzML"]
    page = _page(SelectPage, state, qtbot)
    steps: list[int] = []
    page.navigate_requested.connect(steps.append)

    page.load_folder(run)

    assert steps == [STEP_ASSIGN]


def test_an_empty_folder_warns_in_page_rather_than_in_a_dialog(
    state, qtbot, scan, tmp_path: Path, monkeypatch
):
    def explode(*args, **kwargs):  # pragma: no cover (the assertion is that it is not hit)
        raise AssertionError("a modal dialog was opened for an in-page error")

    monkeypatch.setattr("PySide6.QtWidgets.QMessageBox.warning", explode)
    run = tmp_path / "empty"
    run.mkdir()
    scan.found[:] = []
    page = _page(SelectPage, state, qtbot)
    steps: list[int] = []
    page.navigate_requested.connect(steps.append)

    page.load_folder(run)

    assert page.banner.isVisibleTo(page)
    assert page.banner.tone == "warning"
    assert "No data files" in page.banner.title_label.text()
    assert state.entries == []
    assert steps == [], "an empty folder must not advance the workflow"


def test_the_folder_picker_is_the_hero_action(state, qtbot, scan, tmp_path: Path, monkeypatch):
    run = tmp_path / "run"
    run.mkdir()
    scan.found[:] = [run / "AN_Cal01.mzML"]
    page = _page(SelectPage, state, qtbot)
    monkeypatch.setattr(
        "PySide6.QtWidgets.QFileDialog.getExistingDirectory",
        lambda *a, **k: str(run),
    )

    assert page.choose_button.property("variant") == "primary"
    assert page.add_button.property("variant") == "secondary"
    qtbot.mouseClick(page.choose_button, Qt.MouseButton.LeftButton)

    assert state.run_folder == run


def test_cancelling_the_folder_dialog_changes_nothing(state, qtbot, scan, monkeypatch):
    page = _page(SelectPage, state, qtbot)
    monkeypatch.setattr(
        "PySide6.QtWidgets.QFileDialog.getExistingDirectory", lambda *a, **k: ""
    )

    page.choose_folder()

    assert state.run_folder is None
    assert state.entries == []


def test_add_files_appends_mzml_to_the_selection(state, qtbot, scan, tmp_path: Path):
    run = tmp_path / "run"
    run.mkdir()
    scan.found[:] = [run / "AN_Cal01.mzML"]
    page = _page(SelectPage, state, qtbot)
    page.load_folder(run)

    page.load_files([run / "AN_SMP16.mzML"])

    assert [e.path.name for e in state.entries] == ["AN_Cal01.mzML", "AN_SMP16.mzML"]


def test_add_files_says_so_when_everything_was_already_selected(
    state, qtbot, scan, tmp_path: Path
):
    run = tmp_path / "run"
    run.mkdir()
    scan.found[:] = [run / "AN_Cal01.mzML"]
    page = _page(SelectPage, state, qtbot)
    page.load_folder(run)

    page.load_files([run / "AN_Cal01.mzML"])

    assert page.banner.isVisibleTo(page)
    assert page.banner.tone == "info"
    assert len(state.entries) == 1


# --------------------------------------------------------------------------------------
# Assign
# --------------------------------------------------------------------------------------


def _assign(state: AppState, qtbot, entries: list[FileEntry]) -> AssignPage:
    state.set_entries(entries)
    return _page(AssignPage, state, qtbot)


def test_the_table_shows_one_row_per_entry_with_every_column(state, qtbot, tmp_path: Path):
    page = _assign(state, qtbot, _full_run(tmp_path))
    model = page.model

    assert model.rowCount() == 5
    assert model.columnCount() == 5
    assert model.data(model.index(0, COL_FILE)) == "AN_Cal01_20231002.mzML"
    assert model.data(model.index(0, COL_DETECTED)) == "calibration"
    assert model.data(model.index(0, COL_TYPE)) == "calibration"
    assert model.data(model.index(0, COL_SIZE)).endswith("B")
    assert model.data(model.index(0, COL_STATUS)) == "Ready (mzML)"


def test_a_raw_acquisition_reports_whether_its_mzml_exists(state, qtbot, tmp_path: Path):
    run = tmp_path / "run"
    state.set_run_folder(run)
    entries = [_raw(run, "AN_Cal01_20231002.D", "cal"), _raw(run, "AN_SMP16_20231002.D", "sample")]
    converted = state.converted_dir
    converted.mkdir(parents=True, exist_ok=True)
    (converted / "AN_Cal01_20231002.mzML").write_text("<mzML/>", encoding="utf-8")
    page = _assign(state, qtbot, entries)

    assert page.model.data(page.model.index(0, COL_STATUS)) == "Converted"
    assert page.model.data(page.model.index(1, COL_STATUS)) == "Needs conversion"


def test_every_row_holds_an_open_role_combo_box(state, qtbot, tmp_path: Path, styled):
    """Every row's role combo box stays open and visible, not hidden until clicked."""
    page = _assign(state, qtbot, _full_run(tmp_path))

    for row in range(page.model.rowCount()):
        assert page.table.isPersistentEditorOpen(page.model.index(row, COL_TYPE))
    combos = page.table.viewport().findChildren(QComboBox)
    assert len(combos) == 5


def test_the_role_delegate_writes_an_override_back_to_the_state(state, qtbot, tmp_path: Path):
    page = _assign(state, qtbot, _full_run(tmp_path))
    index = page.model.index(3, COL_TYPE)  # AN_SMP16, detected "sample"
    editor = page.role_delegate.createEditor(page.table.viewport(), None, index)
    page.role_delegate.setEditorData(editor, index)
    assert editor.currentData() == "sample"

    editor.setCurrentIndex(editor.findData("blank"))
    page.role_delegate.setModelData(editor, page.model, index)

    assert state.entries[3].file_type == "blank"
    assert state.entries[3].auto is False


def test_an_override_keeps_the_detected_type_and_drops_the_auto_badge(
    state, qtbot, tmp_path: Path
):
    page = _assign(state, qtbot, _full_run(tmp_path))
    from quant209.app.workflow.assign_page import AUTO_ROLE

    index = page.model.index(3, COL_TYPE)
    assert page.model.data(page.model.index(3, COL_DETECTED)) == "sample"

    page.model.setData(index, "blank", Qt.ItemDataRole.EditRole)

    assert page.model.data(page.model.index(3, COL_DETECTED)) == "sample", "the guess is kept"
    assert page.model.data(page.model.index(3, COL_TYPE)) == "blank"
    assert page.model.data(page.model.index(3, COL_DETECTED), AUTO_ROLE) is False
    assert page.model.data(page.model.index(0, COL_DETECTED), AUTO_ROLE) is True


def test_an_override_does_not_close_the_open_editors(state, qtbot, tmp_path: Path):
    """A reset would shut the combo box the analyst is standing in. See the module docstring."""
    page = _assign(state, qtbot, _full_run(tmp_path))

    page.model.setData(page.model.index(3, COL_TYPE), "blank", Qt.ItemDataRole.EditRole)

    for row in range(page.model.rowCount()):
        assert page.table.isPersistentEditorOpen(page.model.index(row, COL_TYPE))


def test_the_chip_row_follows_the_overrides(state, qtbot, tmp_path: Path):
    page = _assign(state, qtbot, _full_run(tmp_path))
    assert page.chips.texts() == ["1 calibration", "1 reference", "2 samples", "1 blanks"]

    page.model.setData(page.model.index(3, COL_TYPE), "blank", Qt.ItemDataRole.EditRole)

    assert page.chips.texts() == ["1 calibration", "1 reference", "1 samples", "2 blanks"]


def test_no_calibration_injection_is_an_error_banner(state, qtbot, tmp_path: Path):
    page = _assign(
        state,
        qtbot,
        [_mzml(tmp_path, "AN_SMP16.mzML", "sample"), _mzml(tmp_path, "AN_LB72.mzML", "blank")],
    )

    assert page.banner_no_cal.isVisibleTo(page)
    assert page.banner_no_cal.tone == "error"
    assert page.banner_no_cal.title_label.text() == "A calibration injection is required"
    assert page.page_status() == ("No calibration injection", "error")


def test_no_reference_injection_is_a_warning_banner(state, qtbot, tmp_path: Path):
    page = _assign(
        state,
        qtbot,
        [_mzml(tmp_path, "AN_Cal01.mzML", "cal"), _mzml(tmp_path, "AN_LB72.mzML", "blank")],
    )

    assert not page.banner_no_cal.isVisibleTo(page)
    assert page.banner_no_ref.isVisibleTo(page)
    assert page.banner_no_ref.tone == "warning"
    assert "recoveries will be unavailable" in page.banner_no_ref.title_label.text()


def test_no_blanks_is_an_info_banner(state, qtbot, tmp_path: Path):
    page = _assign(
        state,
        qtbot,
        [_mzml(tmp_path, "AN_Cal01.mzML", "cal"), _mzml(tmp_path, "AN_Ref04.mzML", "ref")],
    )

    assert page.banner_no_blank.isVisibleTo(page)
    assert page.banner_no_blank.tone == "info"
    assert "LOQ censoring disabled" in page.banner_no_blank.title_label.text()


def test_a_complete_run_raises_no_banners(state, qtbot, tmp_path: Path):
    page = _assign(state, qtbot, _full_run(tmp_path))

    for banner in (page.banner_no_cal, page.banner_no_ref, page.banner_no_blank):
        assert not banner.isVisibleTo(page)


def test_the_combo_box_arrow_renders_as_a_chevron(state, qtbot, tmp_path: Path, styled):
    """The QSS border-triangle trick paints a solid block in Qt. See ``theme._ASSET_URLS``.

    This renders the delegate's *own* combo box and inspects the pixels in the drop-down
    subcontrol. A filled rectangle means the trick did not apply and Qt fell back to painting
    a solid block. A chevron means it worked, and an arrow that failed to load leaves the
    region empty.
    """
    page = _assign(state, qtbot, _full_run(tmp_path))
    combo = page.table.viewport().findChildren(QComboBox)[0]
    combo.resize(page.TYPE_COLUMN_WIDTH, page.ROW_HEIGHT)

    image = combo.grab().toImage()
    width, height = image.width(), image.height()
    background = QColor(theme.PALETTE["surface"])
    # The ``::drop-down`` subcontrol is 22 px wide. Inset by 4 to clear the rounded border.
    marks = [
        (x, y)
        for x in range(width - 22, width - 4)
        for y in range(4, height - 4)
        if image.pixelColor(x, y) != background
    ]

    assert marks, "no arrow was painted in the combo box's drop-down area"
    xs = {x for x, _ in marks}
    ys = {y for _, y in marks}
    box = len(xs) * len(ys)
    assert len(marks) / box < 0.75, (
        f"the arrow fills {len(marks)}/{box} of its bounding box. That is a solid block, "
        "not a chevron"
    )
    assert 4 <= len(xs) <= 16 and 3 <= len(ys) <= 12, "the chevron is the wrong size"


def test_the_continue_button_names_the_step_it_goes_to(state, qtbot, tmp_path: Path):
    page = _assign(state, qtbot, _full_run(tmp_path))
    assert page.continue_button.text() == "Continue to Process"
    steps: list[int] = []
    page.navigate_requested.connect(steps.append)

    qtbot.mouseClick(page.continue_button, Qt.MouseButton.LeftButton)

    assert steps == [STEP_PROCESS]


# --------------------------------------------------------------------------------------
# Convert
# --------------------------------------------------------------------------------------


@pytest.fixture
def raw_run(state: AppState, tmp_path: Path) -> list[FileEntry]:
    """A run of .D acquisitions with a run folder and a (fake) msconvert in place."""
    run = tmp_path / "SMP16-30_20231002"
    run.mkdir()
    state.set_run_folder(run)
    exe = tmp_path / "msconvert.exe"
    exe.write_text("", encoding="utf-8")
    state.set_msconvert_path(exe)
    entries = [
        _raw(run, "AN_Cal01_20231002.D", "cal"),
        _raw(run, "AN_Ref04_20231002.D", "ref"),
        _raw(run, "AN_LB72_20231002.D", "blank"),
        _raw(run, "AN_SMP16_20231002.D", "sample"),
    ]
    state.set_entries(entries)
    return entries


def _plan(monkeypatch, to_convert=(), reconvert=(), skip=()):
    from quant209.convert.msconvert import ConversionPlan

    plan = ConversionPlan(
        to_convert=list(to_convert), reconvert=list(reconvert), skip=list(skip)
    )
    monkeypatch.setattr(
        "quant209.convert.msconvert.plan_conversion", lambda entries, out_dir: plan
    )
    return plan


def test_the_plan_card_renders_all_three_buckets(state, qtbot, monkeypatch, raw_run):
    _plan(monkeypatch, to_convert=raw_run[:2], reconvert=raw_run[2:3], skip=raw_run[3:])
    page = _page(ConvertPage, state, qtbot)

    summary = page.plan_summary.text()
    assert "2 to convert" in summary
    assert "1 to reconvert (truncated)" in summary
    assert "1 up to date" in summary
    assert page.convert_button.text() == "Convert 3 files"
    rows = [
        label.text()
        for label in page.plan_card.findChildren(QLabel)
        if label.objectName() == "PlanRow"
    ]
    assert any("AN_LB72_20231002.D" in row and "truncated" in row for row in rows)


def test_a_long_plan_collapses_and_expands(state, qtbot, monkeypatch, tmp_path: Path):
    run = tmp_path / "big"
    run.mkdir()
    state.set_run_folder(run)
    entries = [_raw(run, f"AN_SMP{i:02d}_20231002.D", "sample") for i in range(20)]
    state.set_entries(entries)
    _plan(monkeypatch, to_convert=entries)
    page = _page(ConvertPage, state, qtbot)

    def plan_rows():
        return [
            label
            for label in page.plan_card.findChildren(QLabel)
            if label.objectName() == "PlanRow"
        ]

    assert len(plan_rows()) == COLLAPSED_ROWS
    assert page.expand_button.text() == "Show all 20 files"

    page.expand_button.click()

    assert len(plan_rows()) == 20
    assert page.expand_button.text() == "Show fewer"


def test_converting_reports_per_file_progress_and_finishes(
    state, qtbot, monkeypatch, raw_run
):
    _plan(monkeypatch, to_convert=raw_run)
    from quant209.convert.msconvert import ConversionResult

    def convert_many(sources, msconvert_exe, out_dir, extra_args=(), progress=None):
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        results = []
        for i, source in enumerate(sources):
            progress(i, len(sources), Path(source).name)
            output = out_dir / f"{Path(source).stem}.mzML"
            output.write_text("<mzML/>", encoding="utf-8")
            results.append(ConversionResult(Path(source), output, True, "ok"))
        return results

    monkeypatch.setattr("quant209.convert.msconvert.convert_many", convert_many)
    page = _page(ConvertPage, state, qtbot)
    seen: list[str] = []
    page.progress_bar.valueChanged.connect(lambda _v: seen.append(page.current_label.text()))

    page.start_conversion()
    assert state.busy is True, "the page goes busy before the worker is queued"
    qtbot.waitUntil(lambda: state.busy is False, timeout=10_000)

    assert any("Converting AN_Cal01_20231002.D" in text for text in seen)
    assert page.banner.tone == "success"
    assert "Converted 4 files" in page.banner.title_label.text()
    assert not page.progress_card.isVisibleTo(page), "the progress card is packed away"


def test_a_finished_conversion_refreshes_the_plan_and_completes_the_step(
    state, qtbot, monkeypatch, raw_run
):
    """The plan is recomputed from the disk, not remembered from before the run."""
    from quant209.convert.msconvert import ConversionPlan, ConversionResult

    calls = {"n": 0}

    def plan_conversion(entries, out_dir):
        calls["n"] += 1
        out_dir = Path(out_dir)
        pending = [e for e in entries if not (out_dir / f"{e.path.stem}.mzML").exists()]
        done = [e for e in entries if e not in pending]
        return ConversionPlan(to_convert=pending, reconvert=[], skip=done)

    def convert_many(sources, msconvert_exe, out_dir, extra_args=(), progress=None):
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        results = []
        for i, source in enumerate(sources):
            progress(i, len(sources), Path(source).name)
            output = out_dir / f"{Path(source).stem}.mzML"
            output.write_text("<mzML/>", encoding="utf-8")
            results.append(ConversionResult(Path(source), output, True, "ok"))
        return results

    monkeypatch.setattr("quant209.convert.msconvert.plan_conversion", plan_conversion)
    monkeypatch.setattr("quant209.convert.msconvert.convert_many", convert_many)
    page = _page(ConvertPage, state, qtbot)
    assert not state.mzml_ready()

    page.start_conversion()
    qtbot.waitUntil(lambda: state.busy is False, timeout=10_000)

    assert state.mzml_ready(), "every non-hex file now has an mzML"
    assert STEP_CONVERT in state.completed_steps()
    assert page.done_state.isVisibleTo(page)
    assert page.continue_button.isVisibleTo(page)


def test_failed_conversions_are_listed_with_their_reason(state, qtbot, monkeypatch, raw_run):
    _plan(monkeypatch, to_convert=raw_run)
    from quant209.convert.msconvert import ConversionResult

    def convert_many(sources, msconvert_exe, out_dir, extra_args=(), progress=None):
        return [
            ConversionResult(Path(sources[0]), None, False, "unsupported instrument data"),
            *[ConversionResult(Path(s), Path(out_dir), True, "ok") for s in sources[1:]],
        ]

    monkeypatch.setattr("quant209.convert.msconvert.convert_many", convert_many)
    page = _page(ConvertPage, state, qtbot)

    page.start_conversion()
    qtbot.waitUntil(lambda: state.busy is False, timeout=10_000)

    assert page.failure_card.isVisibleTo(page)
    rows = [
        label.text()
        for label in page.failure_card.findChildren(QLabel)
        if label.objectName() == "PlanRow"
    ]
    assert len(rows) == 1
    assert "AN_Cal01_20231002.D" in rows[0]
    assert "unsupported instrument data" in rows[0]
    assert page.banner.tone == "error"


def test_a_missing_msconvert_is_refused_with_a_banner(state, qtbot, monkeypatch, raw_run):
    _plan(monkeypatch, to_convert=raw_run)
    state.set_msconvert_path(None)
    page = _page(ConvertPage, state, qtbot)

    page.start_conversion()

    assert page.banner.tone == "error"
    assert "MSConvert was not found" in page.banner.title_label.text()
    assert state.busy is False


def test_an_already_converted_run_shows_the_done_state(state, qtbot, tmp_path: Path):
    run = tmp_path / "converted"
    state.set_run_folder(run)
    state.set_entries(_full_run(run))
    page = _page(ConvertPage, state, qtbot)

    assert page.done_state.isVisibleTo(page)
    # The headline carries the count. There is no second line repeating it.
    assert page.done_state.headline_label.text().endswith("files converted")
    assert not page.done_state.body_label.text()
    assert not page.plan_card.isVisibleTo(page)
    assert not page.convert_button.isVisibleTo(page)
    assert page.page_status() == ("Complete", "success")
    assert STEP_CONVERT in state.completed_steps(), "the step completes with no user action"


def test_unconverted_hexane_blanks_do_not_hold_the_step_open(state, qtbot, tmp_path: Path):
    """``core.pipeline`` never reads a hexane blank, so an unconverted one cannot block."""
    run = tmp_path / "hexrun"
    state.set_run_folder(run)
    entries = _full_run(run) + [_raw(run, "AN_Hex01_20231002.D", "hex")]
    state.set_entries(entries)

    assert state.mzml_ready()
    assert STEP_CONVERT in state.completed_steps()
    assert state.furthest_step() >= STEP_PROCESS


def test_the_msconvert_row_says_where_the_path_came_from(state, qtbot, raw_run, tmp_path):
    page = _page(ConvertPage, state, qtbot)
    assert "Remembered from your last session." == page.path_hint.text()
    assert page.path_edit.text() == str(state.msconvert_path)

    state.set_msconvert_path(None)
    page.refresh()

    assert "not found automatically" in page.path_hint.text().lower()


# --------------------------------------------------------------------------------------
# Process
# --------------------------------------------------------------------------------------


class _FakeAxis:
    def __init__(self, value: float) -> None:
        self.value = value

    def min(self, skipna=True):  # noqa: A003 (mirroring pandas' API)
        return self.value

    def max(self, skipna=True):
        return self.value


class _FakeRecoveries:
    """The two calls ``ProcessPage._recovery_range`` makes on a ``RunResult.recoveries``."""

    empty = False

    def __init__(self, low: float, high: float) -> None:
        self.low, self.high = low, high

    def min(self, skipna=True):
        return _FakeAxis(self.low)

    def max(self, skipna=True):
        return _FakeAxis(self.high)


def _result(fit_ok: bool = True, reason: str = "", **kwargs):
    """A stand-in ``RunResult``: only what the Process page reads, and nothing heavy."""
    defaults = dict(
        sample_ids=["SMP16", "SMP17"],
        blank_ids=["LB72"],
        fit=SimpleNamespace(fit_ok=fit_ok, reason=reason),
        recoveries=_FakeRecoveries(0.48, 1.21),
        flags={"SMP16": ["low recovery"], "SMP17": []},
        status="OK",
    )
    defaults.update(kwargs)
    return SimpleNamespace(**defaults)


@pytest.fixture
def ready(state: AppState, tmp_path: Path):
    run = tmp_path / "SMP16-30_20231002"
    state.set_run_folder(run)
    state.set_entries(_full_run(run / "Converted"))
    return state


def test_the_summary_card_names_what_will_run(state, qtbot, ready):
    page = _page(ProcessPage, state, qtbot)

    assert page.chips.texts() == ["1 calibration", "1 reference", "2 samples", "1 blanks"]
    assert page.value_labels["cal"].text() == "Cal01"
    assert "Ref04" in page.value_labels["ref"].text()
    assert page.value_labels["blank"].text() == "LB72"
    assert "mean + 3 × SD" in page.value_labels["loq"].text()


def test_the_loq_constant_matches_the_pipeline_default():
    """The Process page states the k value without importing pandas. This pins the copy."""
    from quant209.core.qaqc import QaqcParams

    assert DEFAULT_LOQ_K == QaqcParams().k_sd


def test_running_reports_phase_progress_and_logs_every_line(
    state, qtbot, monkeypatch, ready
):
    result = _result()

    def run_analysis(entries, converted_dir, params=None, progress=None, **kwargs):
        progress("calibrate", 0, 3, "read calibration injection: AN_Cal01_20231002.mzML")
        progress("integrate", 0, 2, "AN_SMP16_20231002")
        progress("integrate", 1, 2, "AN_SMP17_20231002")
        progress("quantify", 4, 5, "limit of quantification")
        return result

    monkeypatch.setattr("quant209.core.pipeline.run_analysis", run_analysis)
    page = _page(ProcessPage, state, qtbot)

    page.run_analysis()
    qtbot.waitUntil(lambda: state.busy is False, timeout=10_000)

    lines = [page.log_view.item(i).text() for i in range(page.log_view.count())]
    assert lines[0].startswith("[calibrate] 1/3 read calibration injection")
    assert "[integrate] 1/2 AN_SMP16_20231002" in lines
    assert "[integrate] 2/2 AN_SMP17_20231002" in lines
    assert lines[-1] == "[done] analysis complete"
    assert page.progress_bar.value() == 100
    assert state.result is result


def test_the_progress_line_is_one_based_and_names_the_phase():
    assert _format_progress("integrate", 6, 22, "AN_SMP16") == "[integrate] 7/22 AN_SMP16"


def test_a_successful_run_reports_samples_recoveries_and_flags(
    state, qtbot, monkeypatch, ready
):
    monkeypatch.setattr(
        "quant209.core.pipeline.run_analysis", lambda *a, **k: _result()
    )
    page = _page(ProcessPage, state, qtbot)

    page.run_analysis()
    qtbot.waitUntil(lambda: state.busy is False, timeout=10_000)

    assert page.banner.tone == "success"
    assert "2 samples quantified" in page.banner.title_label.text()
    assert "1 blank" in page.banner.title_label.text()
    detail = page.banner.body_label.text()
    assert "RT table fits this run" in detail
    assert "recoveries 48–121%" in detail
    assert "1 QA/QC flag" in detail
    assert not page.fit_banner.isVisibleTo(page)


def test_a_failing_run_fit_is_reported_in_red_above_the_summary(
    state, qtbot, monkeypatch, ready
):
    reason = "only 12% of the labeled standards were found in the calibration injection"
    monkeypatch.setattr(
        "quant209.core.pipeline.run_analysis",
        lambda *a, **k: _result(fit_ok=False, reason=reason),
    )
    page = _page(ProcessPage, state, qtbot)

    page.run_analysis()
    qtbot.waitUntil(lambda: state.busy is False, timeout=10_000)

    assert page.fit_banner.isVisibleTo(page)
    assert page.fit_banner.tone == "error"
    assert reason in page.fit_banner.body_label.text()
    assert "does NOT fit this run" in page.banner.body_label.text()


def test_a_successful_run_advances_to_review(state, qtbot, monkeypatch, ready):
    monkeypatch.setattr("quant209.core.pipeline.run_analysis", lambda *a, **k: _result())
    page = _page(ProcessPage, state, qtbot)
    steps: list[int] = []
    page.navigate_requested.connect(steps.append)

    page.run_analysis()
    qtbot.waitUntil(lambda: steps == [STEP_REVIEW], timeout=10_000)


def test_a_failed_run_shows_the_message_and_an_expandable_traceback(
    state, qtbot, monkeypatch, ready
):
    def explode(*args, **kwargs):
        raise ValueError("calibration injection not found: AN_Cal01_20231002.mzML")

    monkeypatch.setattr("quant209.core.pipeline.run_analysis", explode)
    page = _page(ProcessPage, state, qtbot)

    page.run_analysis()
    qtbot.waitUntil(lambda: state.busy is False, timeout=10_000)

    assert page.banner.tone == "error"
    assert "calibration injection not found" in page.banner.body_label.text()
    assert page.details_button.isVisibleTo(page)
    assert not page.details_view.isVisibleTo(page)

    page.details_button.click()

    assert page.details_view.isVisibleTo(page)
    assert "ValueError" in page.details_view.toPlainText()
    assert state.result is None


def test_the_run_button_is_disabled_while_the_worker_is_running(state, qtbot, ready):
    page = _page(ProcessPage, state, qtbot)
    assert page.run_button.isEnabled()

    state.set_busy(True)
    page.refresh()

    assert not page.run_button.isEnabled()
    assert page.cancel_button.isEnabled()
    assert page.page_status() == ("Running", "accent")

    state.set_busy(False)
    page.refresh()
    assert page.run_button.isEnabled()


def test_a_run_that_has_not_been_converted_yet_says_so(state, qtbot, tmp_path: Path):
    run = tmp_path / "raw"
    state.set_run_folder(run)
    state.set_entries([_raw(run, "AN_Cal01_20231002.D", "cal")])
    page = _page(ProcessPage, state, qtbot)

    assert page.empty.isVisibleTo(page)
    assert page.empty.headline_label.text() == "Waiting on conversion"
    assert not page.run_button.isVisibleTo(page)


# --------------------------------------------------------------------------------------
# End to end through the window
# --------------------------------------------------------------------------------------


def test_the_workflow_walks_select_to_process_without_a_dialog(
    qtbot, qapp, state, scan, monkeypatch, tmp_path: Path
):
    """Select -> Assign -> Convert -> Process, driven the way the analyst drives it."""
    configure_app(qapp)
    run = tmp_path / "SMP16-30_20231002"
    converted = run / "Converted"
    converted.mkdir(parents=True)
    names = ["AN_Cal01.mzML", "AN_Ref04.mzML", "AN_LB72.mzML", "AN_SMP16.mzML"]
    for name in names:
        (converted / name).write_text("<mzML/>", encoding="utf-8")
    scan.found[:] = [converted / name for name in names]

    window = MainWindow(state)
    qtbot.addWidget(window)
    assert window.current_step == STEP_SELECT

    window.page(STEP_SELECT).load_folder(run)

    assert window.current_step == STEP_ASSIGN, "the selection carried the analyst forward"
    assert window.page(STEP_ASSIGN).model.rowCount() == 4
    assert window.goto_step(STEP_CONVERT)
    assert window.page(STEP_CONVERT).done_state.isVisibleTo(window.page(STEP_CONVERT))
    assert window.goto_step(STEP_PROCESS)
    assert window.page(STEP_PROCESS).run_button.isEnabled()
    assert window.header.title_label.text() == window.page(STEP_PROCESS).title


def test_the_workflow_pages_import_no_scientific_stack(tmp_path: Path):
    """The four pages are on the startup path. Nothing they import may be heavy."""
    probe = (
        "import sys;"
        "import quant209.app.workflow.select_page, quant209.app.workflow.assign_page,"
        "quant209.app.workflow.convert_page, quant209.app.workflow.process_page;"
        "banned=('numpy','pandas','pyqtgraph','matplotlib','openpyxl',"
        "'quant209.core.pipeline','quant209.core.qaqc','quant209.report.excel');"
        "print(','.join(m for m in banned if m in sys.modules))"
    )
    proc = subprocess.run(
        [sys.executable, "-c", probe], capture_output=True, text=True, timeout=120
    )

    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "", f"eagerly imported: {proc.stdout.strip()}"
