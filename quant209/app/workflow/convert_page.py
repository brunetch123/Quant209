"""Step 3, Convert: run MSConvert over the raw Agilent acquisitions.

The step that is usually already done
--------------------------------------
Most run folders reaching this application have a ``Converted/`` directory in them already,
because conversion often happens as files are acquired, and ``core.sequence.find_data_files``
prefers mzML when it finds any. So the common case for this screen is *nothing to do*, and it
says so in one line and offers the way forward, rather than presenting an empty progress bar
and a disabled button and leaving the analyst to work out that they are finished. The
uncommon case, a folder of raw ``.D`` files or a previous conversion that was interrupted, is
where the plan card earns its place.

Why the plan is shown before anything runs
-------------------------------------------
``convert.msconvert.plan_conversion`` already answers "which of these actually need
converting", in three buckets, and the interesting one is ``reconvert``: an mzML that exists
but is truncated, from a conversion that was killed mid-write. The failure mode it causes (a
file that parses, decodes to nothing, and quietly produces no peaks) is exactly the kind of
thing an analyst should be told about *before* a twenty-minute run rather than after. So the
plan is rendered as a list with a reason against each file, and the primary button says how
much work it is about to do.

What "complete" means here
---------------------------
Every selected file that the pipeline will actually read must have a usable mzML, which
excludes the hexane carryover blanks, since ``core.pipeline.QUANTIFIED_TYPES`` never opens
one. That rule lives in :meth:`~quant209.app.state.AppState.mzml_ready`, not here, so the
stepper's check mark and this page's banner cannot disagree.

MSConvert is a subprocess and can fail per file
------------------------------------------------
``convert_many`` returns a :class:`~quant209.convert.msconvert.ConversionResult` per source,
not an exception, and a run of thirty can come back with two failures and twenty-eight
successes. Those two are listed with the message msconvert gave, because "conversion failed"
without the reason sends the analyst to a log file this application could simply have read.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from PySide6.QtCore import Qt, QThreadPool
from PySide6.QtWidgets import (
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QProgressBar,
    QVBoxLayout,
)

from quant209.app import theme
from quant209.app.state import HEX_TYPE, STEP_CONVERT, STEP_PROCESS
from quant209.app.workers import ConvertWorker
from quant209.app.workflow.base import (
    Banner,
    Card,
    EmptyState,
    WorkflowPage,
    clear_layout,
    ghost_button,
    hint,
    primary_button,
    secondary_button,
)

#: How many plan rows are shown before the list collapses behind "Show all N files".
#: Eight is about a screenful at this row height, enough to see the shape of the plan
#: without the button below it disappearing off the bottom of the page.
COLLAPSED_ROWS = 8

#: Per-bucket wording and tone for a plan row. The reason travels with the filename because
#: "reconvert" without "truncated" would look like the application being indecisive.
PLAN_REASONS: dict[str, tuple[str, str]] = {
    "to_convert": ("to convert", "neutral"),
    "reconvert": ("truncated, reconvert", "warning"),
    "skip": ("up to date", "success"),
}


class ConvertPage(WorkflowPage):
    step = STEP_CONVERT

    def build(self) -> None:
        self._plan: Any = None
        self._rows: list[tuple[str, str, str]] = []  # (name, reason, tone)
        self._expanded = False
        self._worker: ConvertWorker | None = None
        self._failures: list[Any] = []

        self.banner = Banner()
        self.content.addWidget(self.banner)

        self.content.addWidget(self._build_msconvert_card())
        self.content.addWidget(self._build_plan_card())
        self.content.addWidget(self._build_progress_card())
        self.content.addWidget(self._build_failure_card())

        # The fast path this run folder actually takes: nothing to do, said in one line, with
        # the one action that is still meaningful (re-convert on purpose) attached to it.
        self.reconvert_button = secondary_button("Re-convert all")
        self.reconvert_button.clicked.connect(self.reconvert_all)
        self.done_state = EmptyState(
            "All files converted",
            "",
            action=self.reconvert_button,
        )
        self.content.addWidget(self.done_state)

        self.empty = EmptyState(
            "No files to convert",
            "Select a run folder first.",
        )
        self.content.addWidget(self.empty)

        self.content.addStretch(1)
        self.content.addLayout(self._build_actions())

        self.state.msconvert_path_changed.connect(self._on_msconvert_path_changed)
        self.state.entries_changed.connect(lambda _e: self.refresh())
        self.refresh()

    # ---- construction -------------------------------------------------------------------

    def _build_msconvert_card(self) -> Card:
        card = Card("MSConvert")
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(theme.SPACING)
        self.path_edit = QLineEdit()
        self.path_edit.setPlaceholderText("Path to msconvert.exe")
        self.path_edit.setClearButtonEnabled(False)
        # editingFinished, not textChanged: writing a half-typed path to QSettings on every
        # keystroke would persist "C:\Prog" as the remembered executable.
        self.path_edit.editingFinished.connect(self._commit_path)
        row.addWidget(self.path_edit, 1)
        self.browse_button = secondary_button("Browse…")
        self.browse_button.clicked.connect(self.browse_msconvert)
        row.addWidget(self.browse_button)
        card.body.addLayout(row)
        self.path_hint = hint("")
        card.body.addWidget(self.path_hint)
        self.msconvert_card = card
        return card

    def _build_plan_card(self) -> Card:
        card = Card("Conversion plan")
        self.plan_summary = QLabel("")
        self.plan_summary.setObjectName("MetaValue")
        self.plan_summary.setWordWrap(True)
        card.body.addWidget(self.plan_summary)

        self.plan_rows_layout = QVBoxLayout()
        self.plan_rows_layout.setContentsMargins(0, 0, 0, 0)
        self.plan_rows_layout.setSpacing(0)
        card.body.addLayout(self.plan_rows_layout)

        self.expand_button = ghost_button("")
        self.expand_button.clicked.connect(self._toggle_expanded)
        expand_row = QHBoxLayout()
        expand_row.setContentsMargins(0, 0, 0, 0)
        expand_row.addWidget(self.expand_button)
        expand_row.addStretch(1)
        card.body.addLayout(expand_row)

        self.plan_card = card
        return card

    def _build_progress_card(self) -> Card:
        card = Card("Converting")
        self.progress_bar = QProgressBar()
        self.progress_bar.setTextVisible(False)  # the bar is 8 px tall, so the counts are text
        self.progress_bar.setRange(0, 1)
        card.body.addWidget(self.progress_bar)

        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(theme.SPACING)
        self.current_label = QLabel("")
        self.current_label.setObjectName("RunningLabel")
        row.addWidget(self.current_label, 1)
        self.counts_label = QLabel("")
        self.counts_label.setObjectName("MetaLabel")
        row.addWidget(self.counts_label, 0, Qt.AlignmentFlag.AlignRight)
        card.body.addLayout(row)

        cancel_row = QHBoxLayout()
        cancel_row.setContentsMargins(0, 0, 0, 0)
        cancel_row.addStretch(1)
        self.cancel_button = secondary_button("Cancel")
        self.cancel_button.clicked.connect(self.cancel)
        cancel_row.addWidget(self.cancel_button)
        card.body.addLayout(cancel_row)

        card.setVisible(False)
        self.progress_card = card
        return card

    def _build_failure_card(self) -> Card:
        card = Card("Failed")
        self.failure_layout = QVBoxLayout()
        self.failure_layout.setContentsMargins(0, 0, 0, 0)
        self.failure_layout.setSpacing(0)
        card.body.addLayout(self.failure_layout)
        card.setVisible(False)
        self.failure_card = card
        return card

    def _build_actions(self) -> QHBoxLayout:
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(theme.SPACING)
        self.convert_button = primary_button("Convert")
        self.convert_button.clicked.connect(self.start_conversion)
        row.addWidget(self.convert_button)
        row.addStretch(1)
        self.continue_button = primary_button("Continue to Process")
        self.continue_button.clicked.connect(lambda: self.advance_to(STEP_PROCESS))
        row.addWidget(self.continue_button)
        return row

    # ---- msconvert path -----------------------------------------------------------------

    def _commit_path(self) -> None:
        text = self.path_edit.text().strip()
        self.state.set_msconvert_path(Path(text) if text else None)
        self.refresh()

    def browse_msconvert(self) -> None:
        current = self.state.msconvert_path
        start = str(current.parent) if current is not None else str(Path.cwd())
        chosen, _filter = QFileDialog.getOpenFileName(
            self, "Locate msconvert.exe", start, "msconvert (msconvert.exe);;All files (*)"
        )
        if not chosen:
            return
        self.state.set_msconvert_path(Path(chosen))
        self.refresh()

    def _on_msconvert_path_changed(self, _path) -> None:
        self.refresh()

    # ---- the plan -----------------------------------------------------------------------

    def compute_plan(self) -> Any:
        """Ask ``convert.msconvert`` what this run still owes. Lazy import, per the budget."""
        from quant209.convert.msconvert import plan_conversion

        entries = self.state.entries
        out_dir = self.state.converted_dir
        if not entries or out_dir is None:
            return None
        return plan_conversion(entries, out_dir)

    def _plan_row_data(self, plan) -> list[tuple[str, str, str]]:
        """Pending files first, each with the reason it is in the plan."""
        rows: list[tuple[str, str, str]] = []
        for bucket in ("to_convert", "reconvert", "skip"):
            reason, tone = PLAN_REASONS[bucket]
            for entry in getattr(plan, bucket):
                rows.append((entry.path.name, reason, tone))
        return rows

    def _toggle_expanded(self) -> None:
        self._expanded = not self._expanded
        self._render_plan_rows()

    def _render_plan_rows(self) -> None:
        clear_layout(self.plan_rows_layout)
        shown = self._rows if self._expanded else self._rows[:COLLAPSED_ROWS]
        for name, reason, tone in shown:
            label = QLabel(f"{name}   ·   {reason}")
            label.setObjectName("PlanRow")
            label.setProperty("tone", tone)
            self.plan_rows_layout.addWidget(label)
        hidden = len(self._rows) - len(shown)
        if hidden > 0:
            self.expand_button.setText(f"Show all {len(self._rows)} files")
            self.expand_button.setVisible(True)
        elif self._expanded and len(self._rows) > COLLAPSED_ROWS:
            self.expand_button.setText("Show fewer")
            self.expand_button.setVisible(True)
        else:
            self.expand_button.setVisible(False)

    # ---- running ------------------------------------------------------------------------

    def _sources(self, everything: bool = False) -> list[Path]:
        """The ``.D`` paths to hand msconvert: the plan's pending list, or all of them."""
        entries = self.state.entries
        if everything:
            return [e.path for e in entries if e.path.suffix.lower() != ".mzml"]
        if self._plan is None:
            return []
        return [e.path for e in self._plan.pending]

    def start_conversion(self) -> None:
        self._run(self._sources())

    def reconvert_all(self) -> None:
        self._run(self._sources(everything=True))

    def _run(self, sources: list[Path]) -> None:
        if self.state.busy or not sources:
            return
        exe = self.state.msconvert_path
        out_dir = self.state.converted_dir
        if exe is None or not Path(exe).is_file():
            self.banner.show_message(
                "MSConvert was not found",
                "Point Quant 209 at msconvert.exe above. The path is remembered "
                "between sessions.",
                tone="error",
            )
            return
        if out_dir is None:
            self.banner.show_message(
                "There is nowhere to write the converted files",
                "Choose a run folder on the Select step first.",
                tone="error",
            )
            return

        self.banner.clear()
        self._failures = []
        self.failure_card.setVisible(False)
        self.progress_bar.setRange(0, len(sources))
        self.progress_bar.setValue(0)
        self.current_label.setText("Starting MSConvert…")
        self.counts_label.setText(f"0 of {len(sources)}")
        self.progress_card.setVisible(True)

        worker = ConvertWorker(sources, exe, out_dir)
        worker.signals.progress.connect(self._on_progress)
        worker.signals.finished.connect(self._on_finished)
        worker.signals.error.connect(self._on_error)
        worker.signals.cancelled.connect(self._on_cancelled)
        worker.signals.done.connect(self._on_done)
        self._worker = worker
        self.state.set_busy(True)
        self.refresh()
        QThreadPool.globalInstance().start(worker)

    def cancel(self) -> None:
        if self._worker is not None:
            self._worker.cancel()
            self.current_label.setText("Finishing the current file, then stopping…")

    def _on_progress(self, _phase: str, index: int, total: int, name: str) -> None:
        self.progress_bar.setRange(0, total)
        self.progress_bar.setValue(index)
        self.current_label.setText(f"Converting {name}")
        self.counts_label.setText(f"{index + 1} of {total}")

    def _on_finished(self, results) -> None:
        results = list(results or [])
        self._failures = [r for r in results if not getattr(r, "ok", False)]
        converted = len(results) - len(self._failures)
        self.progress_bar.setValue(self.progress_bar.maximum())
        if self._failures:
            self.banner.show_message(
                f"{len(self._failures)} of {len(results)} files could not be converted",
                f"{converted} converted successfully. The failures are listed below with the "
                "message MSConvert gave.",
                tone="error",
            )
        else:
            self.banner.show_message(
                f"Converted {converted} file{'' if converted == 1 else 's'}",
                "Every selected acquisition now has an mzML.",
                tone="success",
            )

    def _on_error(self, message: str, _traceback: str) -> None:
        self.banner.show_message("Conversion failed", message, tone="error")

    def _on_cancelled(self) -> None:
        self.banner.show_message(
            "Conversion cancelled",
            "Files converted before you stopped are finished and were kept.",
            tone="info",
        )

    def _on_done(self) -> None:
        self._worker = None
        self.state.set_busy(False)
        self.progress_card.setVisible(False)
        # The conversion wrote files. ``mzml_ready`` reads the disk, so the workflow gate has
        # to be asked again rather than assumed.
        self.state.refresh()
        self.refresh()

    # ---- rendering ----------------------------------------------------------------------

    def refresh(self) -> None:
        entries = self.state.entries
        self._render_msconvert_row()

        if not entries:
            self.empty.setVisible(True)
            for widget in (
                self.plan_card, self.failure_card, self.msconvert_card, self.done_state
            ):
                widget.setVisible(False)
            for button in (self.convert_button, self.continue_button):
                button.setVisible(False)
            self.status_changed.emit()
            return

        self.empty.setVisible(False)
        self._plan = self.compute_plan()
        pending = list(self._plan.pending) if self._plan is not None else []
        self._rows = self._plan_row_data(self._plan) if self._plan is not None else []

        raw = [e for e in entries if e.path.suffix.lower() != ".mzml"]
        self.msconvert_card.setVisible(bool(raw))
        self.plan_card.setVisible(bool(pending))
        self._render_plan_summary()
        self._render_plan_rows()
        self._render_failures()

        self.done_state.setVisible(not pending)
        self.done_state.headline_label.setText(f"All {len(entries)} files converted")
        self.reconvert_button.setVisible(bool(raw))
        self.reconvert_button.setEnabled(bool(raw) and not self.state.busy)

        busy = self.state.busy
        self.convert_button.setVisible(bool(pending))
        self.convert_button.setEnabled(bool(pending) and not busy)
        self.convert_button.setText(
            f"Convert {len(pending)} file{'' if len(pending) == 1 else 's'}"
        )
        self.continue_button.setVisible(self.state.mzml_ready())
        self.continue_button.setEnabled(not busy)
        self.status_changed.emit()

    def _render_msconvert_row(self) -> None:
        path = self.state.msconvert_path
        text = str(path) if path is not None else ""
        if self.path_edit.text() != text and not self.path_edit.hasFocus():
            self.path_edit.setText(text)
        if path is None:
            self.path_hint.setText(
                "Not found automatically. Browse to msconvert.exe."
            )
        elif self.state.msconvert_auto:
            self.path_hint.setText(f"Found automatically at {path}.")
        else:
            self.path_hint.setText("Remembered from your last session.")

    def _render_plan_summary(self) -> None:
        if self._plan is None:
            self.plan_summary.setText("")
            return
        parts = [
            f"{len(self._plan.to_convert)} to convert",
            f"{len(self._plan.reconvert)} to reconvert (truncated)",
            f"{len(self._plan.skip)} up to date",
        ]
        self.plan_summary.setText("   ·   ".join(parts))

    def _render_failures(self) -> None:
        clear_layout(self.failure_layout)
        for result in self._failures:
            source = getattr(result, "source", None)
            name = Path(source).name if source is not None else "unknown file"
            label = QLabel(f"{name}   ·   {getattr(result, 'message', '')}")
            label.setObjectName("PlanRow")
            label.setProperty("tone", "error")
            label.setWordWrap(True)
            self.failure_layout.addWidget(label)
        self.failure_card.setVisible(bool(self._failures))

    def page_status(self) -> tuple[str, str] | None:
        if not self.state.entries:
            return None
        if self.state.busy:
            return ("Converting", "accent")
        if self.state.mzml_ready():
            return ("Complete", "success")
        pending = [
            e
            for e in self.state.entries
            if e.file_type != HEX_TYPE and not self._entry_ready(e)
        ]
        if self.state.msconvert_path is None:
            return ("MSConvert not found", "warning")
        return (f"{len(pending)} pending", "accent")

    def _entry_ready(self, entry) -> bool:
        mzml = self.state.mzml_for(entry)
        return mzml is not None and mzml.is_file()
