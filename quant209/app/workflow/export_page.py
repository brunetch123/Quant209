"""Step 6, Export: write the analyst-format workbook, and the validation report beside it.

Two outputs, two cards, one folder
-----------------------------------
The Excel workbook is the product: one sheet per reported file, laid out the way an
analyst's own quantification workbooks are, so it can be opened next to last month's and
read without a legend. The validation report is the *evidence*: the automated numbers joined
against the analyst's own workbooks, with the agreement statistics and the verdict. It can
only be written when those workbooks are in the run folder.

So the second card is disabled most of the time, and **a disabled control that does not say
why is a dead end**: it names the folder it looked in and what it was looking for. That is
the same rule the empty states on every other page follow.

Why the folder is remembered
------------------------------
The output folder is the one thing on this screen an analyst sets once and wants forever.
It is persisted through :attr:`~quant209.app.state.AppState.settings` (injected, so the test
suite never writes to a real registry) and defaults to ``<run folder>/output``, beside the
data rather than in the user's documents, because a report is part of the run.

The import discipline
----------------------
``report.excel`` pulls openpyxl and ``report.validation`` pulls matplotlib. Both belong in a
worker's :meth:`~quant209.app.workers._BaseWorker._work`, which is where
:class:`~quant209.app.workers.ReportWorker` and
:class:`~quant209.app.workers.ValidationWorker` put them. This module imports neither, and
:data:`_SKIP_WORKBOOK_TOKEN` is the one fact it duplicates from ``report.truth`` rather than
importing it, pinned by ``tests/test_app_review.py`` so the two cannot drift.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from PySide6.QtCore import QThreadPool, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import QFileDialog, QHBoxLayout, QLabel, QVBoxLayout

from quant209.app import theme
from quant209.app.state import STEP_EXPORT
from quant209.app.workers import ReportWorker, ValidationWorker
from quant209.app.workflow.base import (
    Banner,
    Card,
    EmptyState,
    WorkflowPage,
    ghost_button,
    hint,
    meta_label,
    path_label,
    primary_button,
    secondary_button,
)

#: Where reports go inside a run folder, matching ``quant209.cli.DEFAULT_OUT_DIR`` so a run
#: exported from the GUI and one exported from the CLI land in the same place.
DEFAULT_OUT_DIRNAME = "output"

#: The workbook name pattern, matching ``quant209.cli.cmd_run``: ``<run folder>_report.xlsx``.
REPORT_SUFFIX = "_report.xlsx"

#: Where the validation report's eight files go, under the output folder. Mirrors the CLI's
#: ``<out>/validation/<run name>/``.
VALIDATION_DIRNAME = "validation"

#: ``report.truth._find_quant_workbooks``'s two exclusions, duplicated rather than imported
#: (importing it pulls openpyxl and pandas onto a page that only wants to know whether a
#: button should be enabled). ``tests/test_app_review.py`` pins this scan to that function.
_LOCKFILE_PREFIX = "~$"
_SKIP_WORKBOOK_TOKEN = "sequence"

_SETTINGS_EXPORT_DIR = "paths/export_dir"


def truth_workbooks(run_folder: Path | None) -> list[Path]:
    """The analyst quantification workbooks in ``run_folder``, which is what enables validation.

    A filename scan, not a parse: this runs on every refresh, and opening twenty workbooks
    to decide whether a button is enabled would make the page feel broken. The worker does
    the real read and reports honestly if the scan was optimistic.
    """
    if run_folder is None:
        return []
    try:
        candidates = sorted(run_folder.glob("*.xlsx"))
    except OSError:  # pragma: no cover (an unreachable network path)
        return []
    return [
        path
        for path in candidates
        if not path.name.startswith(_LOCKFILE_PREFIX)
        and _SKIP_WORKBOOK_TOKEN not in path.name.lower()
    ]


class ExportCard(Card):
    """One output: what it is, what it will be called, one button, and its outcome.

    A card rather than a row of buttons because each output has a *state* (ready, running,
    written, failed), and a button that has to carry all four in its label is the control an
    analyst presses twice.
    """

    def __init__(self, title: str, body: str, action_text: str, primary: bool = True) -> None:
        super().__init__(title)
        self.description = hint(body)
        self.body.addWidget(self.description)

        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(theme.SPACING)
        row.addWidget(meta_label("WRITES"))
        self.target_label = QLabel("")
        self.target_label.setObjectName("MetaValue")
        self.target_label.setWordWrap(True)
        row.addWidget(self.target_label, 1)
        self.body.addLayout(row)

        self.banner = Banner()
        self.body.addWidget(self.banner)

        actions = QHBoxLayout()
        actions.setContentsMargins(0, theme.SPACING // 2, 0, 0)
        actions.setSpacing(theme.SPACING)
        self.button = (primary_button if primary else secondary_button)(action_text)
        actions.addWidget(self.button)
        self.open_button = ghost_button("Open folder")
        self.open_button.setVisible(False)
        actions.addWidget(self.open_button)
        actions.addStretch(1)
        self.body.addLayout(actions)

        self.written: Path | None = None

    def set_ready(self, target: str, enabled: bool, reason: str = "") -> None:
        self.target_label.setText(target)
        self.button.setEnabled(enabled)
        if reason:
            self.description.setText(reason)

    def set_running(self, message: str) -> None:
        self.button.setEnabled(False)
        self.open_button.setVisible(False)
        self.banner.show_message(message, "", tone="info")

    def set_written(self, headline: str, path: Path) -> None:
        self.written = path
        self.banner.show_message(headline, str(path), tone="success")
        self.open_button.setVisible(True)

    def set_failed(self, message: str) -> None:
        self.written = None
        self.open_button.setVisible(False)
        self.banner.show_message("That did not work", message, tone="error")


class ExportPage(WorkflowPage):
    step = STEP_EXPORT
    #: The two cards and their success banners grow past the window on a laptop screen once
    #: both have run. This is the one page whose content is genuinely taller than its frame.
    scrollable = True

    def build(self) -> None:
        self._report_worker: ReportWorker | None = None
        self._validation_worker: ValidationWorker | None = None

        self.folder_card = Card("Output folder")
        self.folder_label = path_label("")
        self.folder_card.body.addWidget(self.folder_label)
        self.folder_card.body.addWidget(
            hint(
                "Remembered between runs. Both outputs are written here. The validation "
                f"report goes in a {VALIDATION_DIRNAME}/ subfolder beside the workbook."
            )
        )
        folder_row = QHBoxLayout()
        folder_row.setContentsMargins(0, 0, 0, 0)
        folder_row.setSpacing(theme.SPACING)
        self.choose_button = secondary_button("Change folder…")
        self.choose_button.clicked.connect(self.choose_folder)
        folder_row.addWidget(self.choose_button)
        folder_row.addStretch(1)
        self.folder_card.body.addLayout(folder_row)
        self.content.addWidget(self.folder_card)

        self.report_card = ExportCard(
            "Excel report",
            "One sheet per sample and blank, plus cal, ref, Summary and "
            "SamplesList_Standards sheets. A run_info sheet records how the numbers "
            "were produced.",
            "Write Excel report",
        )
        self.report_card.button.clicked.connect(self.write_report)
        self.report_card.open_button.clicked.connect(
            lambda: self.open_folder(self.report_card.written)
        )
        self.content.addWidget(self.report_card)

        self.validation_card = ExportCard(
            "Validation report",
            "",
            "Write validation report",
            primary=False,
        )
        self.validation_card.button.clicked.connect(self.write_validation)
        self.validation_card.open_button.clicked.connect(
            lambda: self.open_folder(self.validation_card.written)
        )
        self.content.addWidget(self.validation_card)

        self.empty = EmptyState(
            "Nothing to export yet",
            "Process the run first.",
        )
        self.content.addWidget(self.empty)
        self.content.addStretch(1)

        self.state.result_changed.connect(lambda _r: self.refresh())
        self.refresh()

    # ---- the output folder ------------------------------------------------------------

    def output_dir(self) -> Path:
        """Where reports go: the remembered folder, else ``<run folder>/output``."""
        stored = self.state.settings.value(_SETTINGS_EXPORT_DIR, "")
        if stored:
            return Path(str(stored))
        run_folder = self.state.run_folder
        base = run_folder if run_folder is not None else Path.cwd()
        return base / DEFAULT_OUT_DIRNAME

    def set_output_dir(self, folder: Path | str) -> None:
        self.state.settings.setValue(_SETTINGS_EXPORT_DIR, str(folder))
        self.state.settings.sync()
        self.refresh()

    def choose_folder(self) -> None:
        chosen = QFileDialog.getExistingDirectory(
            self, "Choose the output folder", str(self.output_dir())
        )
        if chosen:
            self.set_output_dir(chosen)

    def report_path(self) -> Path:
        run_folder = self.state.run_folder
        name = run_folder.name if run_folder is not None else "run"
        return self.output_dir() / f"{name}{REPORT_SUFFIX}"

    def validation_dir(self) -> Path:
        run_folder = self.state.run_folder
        name = run_folder.name if run_folder is not None else "run"
        return self.output_dir() / VALIDATION_DIRNAME / name

    @staticmethod
    def open_folder(path: Path | None) -> None:
        """Reveal a written output in the platform file manager.

        The *folder*, never the file: opening the workbook itself would hand a half-read
        Excel process a file the analyst may be about to overwrite, and the reason to press
        this is almost always "where did that go".
        """
        if path is None:
            return
        folder = path if path.is_dir() else path.parent
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(folder)))

    # ---- writing --------------------------------------------------------------------------

    def write_report(self) -> None:
        result = self.state.result
        if result is None or self.state.busy:
            return
        out_path = self.report_path()
        run_folder = self.state.run_folder
        worker = ReportWorker(result, out_path, run_folder.name if run_folder else None)
        worker.signals.finished.connect(self._on_report_written)
        worker.signals.error.connect(lambda message, _tb: self.report_card.set_failed(message))
        worker.signals.done.connect(self._on_done)
        self._report_worker = worker
        self.report_card.set_running(f"Writing {out_path.name}…")
        self.state.set_busy(True)
        self.refresh()
        QThreadPool.globalInstance().start(worker)

    def _on_report_written(self, path: Any) -> None:
        written = Path(path) if path is not None else self.report_path()
        self.report_card.set_written("Report written", written)
        # Recorded on the state, which is what marks the Export step complete in the
        # stepper, the one place in the application that says the run is finished.
        self.state.set_export_path(written)

    def write_validation(self) -> None:
        result = self.state.result
        run_folder = self.state.run_folder
        if result is None or run_folder is None or self.state.busy:
            return
        out_dir = self.validation_dir()
        worker = ValidationWorker(result, run_folder, out_dir, run_folder.name)
        worker.signals.finished.connect(self._on_validation_written)
        worker.signals.error.connect(
            lambda message, _tb: self.validation_card.set_failed(message)
        )
        worker.signals.done.connect(self._on_done)
        self._validation_worker = worker
        self.validation_card.set_running("Reading the analyst's workbooks…")
        self.state.set_busy(True)
        self.refresh()
        QThreadPool.globalInstance().start(worker)

    def _on_validation_written(self, payload: Any) -> None:
        out_dir, report = payload if isinstance(payload, tuple) else (payload, None)
        passed = bool(getattr(report, "passed", False))
        self.validation_card.set_written(
            f"Validation report written: {'PASS' if passed else 'FAIL'}"
            if report is not None
            else "Validation report written",
            Path(out_dir),
        )

    def _on_done(self) -> None:
        self._report_worker = None
        self._validation_worker = None
        self.state.set_busy(False)
        self.refresh()

    # ---- rendering --------------------------------------------------------------------------

    def refresh(self) -> None:
        result = self.state.result
        has_result = result is not None
        busy = self.state.busy

        self.folder_card.setVisible(has_result)
        self.report_card.setVisible(has_result)
        self.validation_card.setVisible(has_result)
        self.empty.setVisible(not has_result)
        if not has_result:
            self.status_changed.emit()
            return

        self.folder_label.setText(str(self.output_dir()))
        self.choose_button.setEnabled(not busy)
        self.report_card.set_ready(self.report_path().name, enabled=not busy)

        workbooks = truth_workbooks(self.state.run_folder)
        run_name = self.state.run_folder.name if self.state.run_folder else "this run"
        if workbooks:
            reason = (
                f"Joins this run's numbers against the analyst's own quantification "
                f"{'workbook' if len(workbooks) == 1 else 'workbooks'} "
                f"({', '.join(p.name for p in workbooks[:3])}"
                f"{', …' if len(workbooks) > 3 else ''}) and reports the agreement. This covers "
                "median |RPD| on areas, masses, recoveries and response factors, plus "
                "the verdict."
            )
        else:
            reason = (
                f"Unavailable: no analyst quantification workbook was found in {run_name}. "
                "The validation report compares this run against the analyst's own numbers, "
                "so it needs their .xlsx workbooks in the run folder."
            )
        self.validation_card.set_ready(
            f"{VALIDATION_DIRNAME}/{run_name}/summary.md + 5 CSV/PDF files",
            enabled=bool(workbooks) and not busy,
            reason=reason,
        )
        self.status_changed.emit()

    def page_status(self) -> tuple[str, str] | None:
        if self.state.busy:
            return ("Writing", "accent")
        if self.state.export_path is not None:
            return ("Written", "success")
        if self.state.result is not None:
            return ("Ready", "accent")
        return None
