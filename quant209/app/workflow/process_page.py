"""Step 4, Process: integrate and quantify the run.

The one screen that narrates while it runs
--------------------------------------------
A run is twenty-odd mzML files, each read, integrated at this run's anchors and released
again, on a host with ~2 GB of RAM. It takes minutes, and a program that appears to have
stopped is slower than the same program narrating itself, whatever the clock says. So this
page reports three things at once, and none of them is a spinner:

* **An overall bar**, weighted across the pipeline's three phases rather than restarting at
  zero for each one. A bar that reaches 100% twice is a bar nobody trusts again.
* **A phase line** in the pipeline's own vocabulary (``[integrate] 7/22 LM_BRT16``), so the
  number moving is a file the analyst can name.
* **A log** that keeps every line, because "which file was it on when it slowed down" is a
  question asked afterwards, and a label that only ever shows the current value cannot
  answer it.

The run itself is on a :class:`~quant209.app.workers.PipelineWorker`, so the window stays
responsive throughout and Cancel is a button rather than the task manager.

Reading the result without importing pandas
--------------------------------------------
This module must not import ``core.pipeline`` at module scope. That is the cold-start rule
the whole ``workflow`` package is built around. The ``RunResult`` that comes back is
therefore read through ``getattr`` with fallbacks: by the time it exists, pandas is loaded
(the worker did it), but this file still must not be the reason it loads. The same discipline
is why :data:`DEFAULT_LOQ_K` is a literal here instead of ``QaqcParams().k_sd``. A test pins
the two together so they cannot drift.
"""

from __future__ import annotations

from typing import Any

from PySide6.QtCore import Qt, QThreadPool
from PySide6.QtWidgets import (
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPlainTextEdit,
    QProgressBar,
    QSizePolicy,
)

from quant209.app import theme
from quant209.app.state import STEP_PROCESS, STEP_REVIEW
from quant209.app.workers import PipelineWorker
from quant209.app.workflow.base import (
    Banner,
    Card,
    ChipRow,
    EmptyState,
    WorkflowPage,
    ghost_button,
    meta_label,
    meta_value,
    primary_button,
    role_chips,
    secondary_button,
    set_variant,
)

#: ``core.qaqc.QaqcParams.k_sd``, repeated rather than imported (importing ``core.qaqc``
#: pulls pandas into the startup path). ``tests/test_app_pages.py`` asserts they agree.
DEFAULT_LOQ_K = 3.0

#: Where each pipeline phase sits on the overall bar, as ``(start, end)`` percentages.
#: Integration is four fifths of the bar because it is very nearly all of the wall clock.
#: One mzML read and integrated per file, against a fixed handful of steps either side.
PHASE_SPANS: dict[str, tuple[float, float]] = {
    "calibrate": (0.0, 8.0),
    "integrate": (8.0, 88.0),
    "quantify": (88.0, 100.0),
}

#: How long the success banner is left on screen before the application moves to Review.
#: Long enough to register that the run finished and read the verdict's colour. Short enough
#: that it is not a pause the analyst has to sit through.
ADVANCE_DELAY_MS = 600

#: The primary button's two labels. Named because the button reserves the width of the
#: longer one at build time. See ``ProcessPage._build_actions``.
RUN_TEXT = "Run analysis"
RUN_AGAIN_TEXT = "Run analysis again"

#: The run log keeps every line. This is only a ceiling against a pathological progress
#: callback, not a working limit (a 30-file run produces about forty lines).
MAX_LOG_LINES = 2000


def _format_progress(phase: str, index: int, total: int, label: str) -> str:
    """The pipeline's own callback (``[integrate] 7/22 LM_BRT16``), rendered for a human.

    ``index`` arrives 0-based (``core.pipeline``'s contract) and is shown 1-based, because
    "0/22" next to a file that is currently being read is simply wrong to a reader.
    """
    return f"[{phase}] {min(index + 1, total)}/{total} {label}".rstrip()


class ProcessPage(WorkflowPage):
    step = STEP_PROCESS

    def build(self) -> None:
        self._worker: PipelineWorker | None = None
        self._shown_result: Any = None
        self._traceback = ""

        self.banner = Banner()
        self.fit_banner = Banner()
        self.content.addWidget(self.fit_banner)
        self.content.addWidget(self.banner)

        self.details_view = QPlainTextEdit()
        self.details_view.setObjectName("Traceback")
        self.details_view.setReadOnly(True)
        self.details_view.setFont(theme.mono_font())
        self.details_view.setMaximumHeight(200)
        self.details_view.setVisible(False)
        self.content.addWidget(self.details_view)

        self.content.addWidget(self._build_summary_card())
        self.content.addWidget(self._build_run_card(), 1)

        self.empty = EmptyState("Nothing to process", "")
        self.content.addWidget(self.empty)
        self.content.addStretch(1)
        self.content.addLayout(self._build_actions())

        self.state.entries_changed.connect(lambda _e: self.refresh())
        self.state.result_changed.connect(lambda _r: self.refresh())
        self.refresh()

    # ---- construction -------------------------------------------------------------------

    def _build_summary_card(self) -> Card:
        card = Card("Run summary")
        self.chips = ChipRow()
        card.body.addWidget(self.chips)

        grid = QGridLayout()
        grid.setContentsMargins(0, theme.SPACING // 2, 0, 0)
        grid.setHorizontalSpacing(theme.SPACING * 2)
        grid.setVerticalSpacing(theme.SPACING // 2 + 2)
        grid.setColumnStretch(1, 1)
        self.value_labels: dict[str, QLabel] = {}
        for row, (key, caption) in enumerate(
            (
                ("cal", "CALIBRATION"),
                ("ref", "REFERENCE"),
                ("blank", "BLANKS"),
                ("loq", "LIMIT OF QUANTIFICATION"),
            )
        ):
            grid.addWidget(meta_label(caption), row, 0, Qt.AlignmentFlag.AlignTop)
            value = meta_value("")
            self.value_labels[key] = value
            grid.addWidget(value, row, 1)
        card.body.addLayout(grid)
        self.summary_card = card
        return card

    def _build_run_card(self) -> Card:
        card = Card("Running")
        self.progress_bar = QProgressBar()
        self.progress_bar.setTextVisible(False)
        self.progress_bar.setRange(0, 100)
        card.body.addWidget(self.progress_bar)

        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(theme.SPACING)
        self.phase_label = QLabel("")
        self.phase_label.setObjectName("RunningLabel")
        self.phase_label.setFont(theme.mono_font(theme.FONT_SIZE_PT))
        row.addWidget(self.phase_label, 1)
        self.cancel_button = secondary_button("Cancel")
        self.cancel_button.clicked.connect(self.cancel)
        row.addWidget(self.cancel_button)
        card.body.addLayout(row)

        self.log_view = QListWidget()
        self.log_view.setObjectName("RunLog")
        self.log_view.setFont(theme.mono_font())
        self.log_view.setSelectionMode(QListWidget.SelectionMode.NoSelection)
        self.log_view.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        # 140, not 180: at the minimum window size the banner, the summary card and this one
        # together came to a few pixels more than the content area, and the overflow was paid
        # by the bottom of this card, clipping the newest log line, the one being read.
        # The list is the elastic part of the page, so it also takes any spare height.
        self.log_view.setMinimumHeight(140)
        self.log_view.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        card.body.addWidget(self.log_view, 1)

        card.setVisible(False)
        self.run_card = card
        return card

    def _build_actions(self) -> QHBoxLayout:
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(theme.SPACING)
        self.run_button = primary_button(RUN_AGAIN_TEXT)
        # Reserve the width of the longer of the two labels before showing the shorter one.
        # Without it the button visibly re-flows (and, for the frame between the text
        # changing and the layout catching up, clips its own label) at the exact moment the
        # run finishes, which is the moment the analyst is looking straight at it.
        self.run_button.setMinimumWidth(self.run_button.sizeHint().width())
        self.run_button.setText(RUN_TEXT)
        self.run_button.clicked.connect(self.run_analysis)
        row.addWidget(self.run_button)
        self.details_button = ghost_button("Show details")
        self.details_button.clicked.connect(self._toggle_details)
        self.details_button.setVisible(False)
        row.addWidget(self.details_button)
        row.addStretch(1)
        self.review_button = primary_button("Go to Review")
        self.review_button.clicked.connect(lambda: self.advance_to(STEP_REVIEW))
        self.review_button.setVisible(False)
        row.addWidget(self.review_button)
        return row

    # ---- running ------------------------------------------------------------------------

    def run_analysis(self) -> None:
        """Start the pipeline on a worker thread. The window stays live throughout."""
        if self.state.busy:
            return
        entries = self.state.entries
        converted_dir = self.state.converted_dir or (
            self.state.run_folder if self.state.run_folder is not None else "."
        )
        self.banner.clear()
        self.fit_banner.clear()
        self.details_view.setVisible(False)
        self.details_button.setVisible(False)
        self.log_view.clear()
        self.progress_bar.setValue(0)
        self.phase_label.setText("[calibrate] starting")
        self.run_card.setVisible(True)

        worker = PipelineWorker(entries, converted_dir)
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
            self.phase_label.setText("Stopping at the next file…")

    def _on_progress(self, phase: str, index: int, total: int, label: str) -> None:
        line = _format_progress(phase, index, total, label)
        self.phase_label.setText(line)
        self.progress_bar.setValue(int(self._overall(phase, index, total)))
        self._append_log(line)

    @staticmethod
    def _overall(phase: str, index: int, total: int) -> float:
        start, end = PHASE_SPANS.get(phase, (0.0, 100.0))
        fraction = (index / total) if total else 0.0
        return start + fraction * (end - start)

    def _append_log(self, line: str) -> None:
        item = QListWidgetItem(line)
        self.log_view.addItem(item)
        while self.log_view.count() > MAX_LOG_LINES:
            self.log_view.takeItem(0)
        # ``scrollToItem``, not ``scrollToBottom``: the latter scrolls to the range the view
        # knew about before this item was laid out, which leaves the newest line (the only
        # one anybody is looking at) half below the fold.
        self.log_view.scrollToItem(item, QListWidget.ScrollHint.PositionAtBottom)

    def _on_finished(self, result) -> None:
        self.progress_bar.setValue(100)
        self._append_log("[done] analysis complete")
        self.state.set_result(result)
        self._render_result_banner(result)
        self._shown_result = result
        # A beat on the verdict, then Review. See ``WorkflowPage.advance_to``.
        self.advance_to(STEP_REVIEW, delay_ms=ADVANCE_DELAY_MS)

    def _on_error(self, message: str, traceback_text: str) -> None:
        self._shown_result = None
        self._traceback = traceback_text
        self.fit_banner.clear()
        self.banner.show_message("The run could not be completed", message, tone="error")
        self.details_view.setPlainText(traceback_text)
        self.details_button.setText("Show details")
        self.details_button.setVisible(True)
        self._append_log(f"[error] {message}")

    def _on_cancelled(self) -> None:
        self.banner.show_message(
            "Run cancelled",
            "Nothing was written and no result was kept. A partial run does not produce "
            "a partial result.",
            tone="info",
        )
        self._append_log("[cancelled] stopped at the analyst's request")

    def _on_done(self) -> None:
        self._worker = None
        self.state.set_busy(False)
        self.refresh()

    def _toggle_details(self) -> None:
        showing = not self.details_view.isVisible()
        self.details_view.setVisible(showing)
        self.details_button.setText("Hide details" if showing else "Show details")

    # ---- the verdict --------------------------------------------------------------------

    def _render_result_banner(self, result) -> None:
        """The green summary, and, if there was a fit failure, a red banner above it.

        The fit report is the one thing on this screen that can invalidate every number
        behind it, so it gets its own banner rather than a clause in a sentence. A run whose
        retention-time table does not describe its own calibration injection still produced
        numbers, and they are still shown. What must not happen is the analyst reading them
        without knowing.
        """
        fit = getattr(result, "fit", None)
        fit_ok = bool(getattr(fit, "fit_ok", True))
        if fit_ok:
            self.fit_banner.clear()
        else:
            reason = getattr(fit, "reason", "") or "the retention-time table does not fit this run"
            self.fit_banner.show_message(
                "The retention-time table does not fit this run", reason, tone="error"
            )
        self.banner.show_message(
            self._result_headline(result), self._result_detail(result, fit_ok), tone="success"
        )

    @staticmethod
    def _result_headline(result) -> str:
        samples = len(getattr(result, "sample_ids", []) or [])
        blanks = len(getattr(result, "blank_ids", []) or [])
        parts = [f"{samples} sample{'' if samples == 1 else 's'} quantified"]
        if blanks:
            parts.append(f"{blanks} blank{'' if blanks == 1 else 's'}")
        return "   ·   ".join(parts)

    def _result_detail(self, result, fit_ok: bool) -> str:
        parts = ["RT table fits this run" if fit_ok else "RT table does NOT fit this run"]
        recovery = self._recovery_range(result)
        if recovery:
            parts.append(f"recoveries {recovery}")
        flags = getattr(result, "flags", None) or {}
        try:
            count = sum(len(v) for v in flags.values())
        except (AttributeError, TypeError):  # a stub result, or flags of an unexpected shape
            count = 0
        parts.append(f"{count} QA/QC flag{'' if count == 1 else 's'}")
        return "   ·   ".join(parts)

    @staticmethod
    def _recovery_range(result) -> str:
        """The recovery range formatted as text, for example "48 to 121%", or "" if there
        are none.

        Wrapped in a broad except because this is a banner subtitle, and a display nicety
        must never be the thing that turns a finished run into an error dialog.
        """
        frame = getattr(result, "recoveries", None)
        if frame is None:
            return ""
        try:
            if getattr(frame, "empty", True):
                return ""
            low = float(frame.min(skipna=True).min(skipna=True))
            high = float(frame.max(skipna=True).max(skipna=True))
        except Exception:  # noqa: BLE001 (see the docstring)
            return ""
        if low != low or high != high:  # NaN
            return ""
        return f"{low * 100:.0f}–{high * 100:.0f}%"

    # ---- rendering ----------------------------------------------------------------------

    def refresh(self) -> None:
        entries = self.state.entries
        ready = bool(entries) and self.state.mzml_ready()
        busy = self.state.busy
        result = self.state.result

        self.summary_card.setVisible(ready)
        self.empty.setVisible(not ready)
        if not ready:
            self._render_empty(entries)
        else:
            self._render_summary(entries)

        self.run_button.setVisible(ready)
        self.run_button.setEnabled(ready and not busy)
        self.run_button.setText(RUN_AGAIN_TEXT if result is not None else RUN_TEXT)
        # Once there is a result the hero action is Review, not another run.
        set_variant(self.run_button, "secondary" if result is not None else "primary")
        self.cancel_button.setEnabled(busy)
        self.review_button.setVisible(result is not None and not busy)
        # The card is a live progress readout during a run and a record of one afterwards.
        self.run_card.title_label.setText("Running" if busy else "Run log")

        if result is not None and result is not self._shown_result and not busy:
            # A result that arrived while this page was not on screen (or a page rebuilt
            # around an existing one): show the verdict rather than a blank screen.
            self._render_result_banner(result)
            self._shown_result = result
        self.status_changed.emit()

    def _render_empty(self, entries) -> None:
        if not entries:
            self.empty.headline_label.setText("Nothing to process")
            self.empty.body_label.setText("Select a run folder first.")
        else:
            pending = [e for e in self.state.entries if not self._entry_ready(e)]
            self.empty.headline_label.setText("Waiting on conversion")
            self.empty.body_label.setText(
                f"{len(pending)} of {len(entries)} files not converted."
            )
        self.empty.body_label.setVisible(True)

    def _entry_ready(self, entry) -> bool:
        mzml = self.state.mzml_for(entry)
        return mzml is not None and mzml.is_file()

    def _render_summary(self, entries) -> None:
        self.chips.set_chips(role_chips(entries))
        by_role: dict[str, list[str]] = {}
        for entry in entries:
            by_role.setdefault(entry.file_type, []).append(entry.sample_id)

        cal = by_role.get("cal", [])
        self.value_labels["cal"].setText(
            ", ".join(cal) if cal else "none, so this run cannot be quantified"
        )
        ref = by_role.get("ref", [])
        self.value_labels["ref"].setText(
            ", ".join(ref)
            if ref
            else "none"
        )
        blanks = by_role.get("blank", [])
        self.value_labels["blank"].setText(
            ", ".join(blanks) if blanks else "none"
        )
        self.value_labels["loq"].setText(
            f"mean + {DEFAULT_LOQ_K:g} × SD of blanks"
        )

    def page_status(self) -> tuple[str, str] | None:
        if self.state.busy:
            return ("Running", "accent")
        if self.state.result is not None:
            return ("Complete", "success")
        if self.state.entries and self.state.mzml_ready():
            return ("Ready", "accent")
        return None
