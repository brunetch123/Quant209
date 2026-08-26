"""Step 1, Select: choose the analysis-run folder.

The one decision on this screen
--------------------------------
An analysis run is a *folder*. Everything the application needs, including the calibration
injections, the reference injections, the blanks, the samples, and (usually) the analyst's
own workbooks, is inside one, named for the batch and the date. So the hero action is a
folder picker and nothing else competes with it: one large primary button, one sentence
saying what a run folder is, and a quiet secondary path for the case the folder convention
does not cover.

Why "Add files..." is secondary and only takes .mzML
-----------------------------------------------------
A native open-file dialog cannot select an Agilent ``.D`` *directory*. It is a folder, and
a file dialog in file mode will not return one. Offering it as an equal choice would
therefore promise something it cannot do for the format most runs actually arrive in. It
stays for the genuine case it does serve: mzML files gathered from more than one place, or a
single file being re-checked.

Errors are banners, not dialogs
--------------------------------
"This folder has no data files in it" is information about the folder the user is looking
at, and it belongs on the screen next to the folder's name, where it can be read twice. A
modal would say it once, demand an acknowledgement, and leave the screen looking exactly as
it did before.

The heavy imports
------------------
``core.sequence`` is cheap (``re`` and dataclasses), but it is still imported inside the
handler rather than at module scope: this is the one page a launch always builds, and the
rule that keeps the cold-start budget intact is easier to keep when it has no exceptions.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QFileDialog, QHBoxLayout, QVBoxLayout, QWidget

from quant209.app import theme
from quant209.app.state import STEP_ASSIGN, STEP_SELECT
from quant209.app.workflow.base import (
    Banner,
    Card,
    ChipRow,
    EmptyState,
    WorkflowPage,
    hint,
    meta_label,
    meta_value,
    path_label,
    primary_button,
    role_chips,
    secondary_button,
)

#: What the analyst is told a run folder is, on the empty screen. One paragraph, because an
#: empty state should explain what to do next rather than shrugging.
EMPTY_BODY = ""

#: The mzML-only filter on the secondary picker. See the module docstring.
MZML_FILTER = "mzML files (*.mzML);;All files (*)"


class SelectPage(WorkflowPage):
    step = STEP_SELECT

    # ---- construction -------------------------------------------------------------------

    def build(self) -> None:
        self.banner = Banner()
        self.content.addWidget(self.banner)

        self.choose_button = primary_button("Choose run folder…")
        self.choose_button.clicked.connect(self.choose_folder)
        self.add_button = secondary_button("Add files…")
        self.add_button.clicked.connect(self.add_files)

        actions = QWidget()
        actions_row = QHBoxLayout(actions)
        actions_row.setContentsMargins(0, 0, 0, 0)
        actions_row.setSpacing(theme.SPACING + 4)
        actions_row.addWidget(self.choose_button)
        actions_row.addWidget(self.add_button)

        self.empty = EmptyState("No run folder selected", EMPTY_BODY, action=actions)
        self.content.addWidget(self.empty)

        self.content.addWidget(self._build_summary())
        self.content.addStretch(1)
        self.refresh()

    def _build_summary(self) -> Card:
        """The card that replaces the empty state once a folder has been read."""
        self.summary = Card("Run folder")
        self.summary.setVisible(False)

        self.folder_label = path_label("")
        self.summary.body.addWidget(self.folder_label)

        counts = QHBoxLayout()
        counts.setContentsMargins(0, 0, 0, 0)
        counts.setSpacing(theme.SPACING)
        counts.addWidget(meta_label("DATA FILES"))
        self.counts_label = meta_value("")
        counts.addWidget(self.counts_label, 1)
        self.summary.body.addLayout(counts)

        types = QVBoxLayout()
        types.setContentsMargins(0, 0, 0, 0)
        types.setSpacing(theme.SPACING // 2)
        types.addWidget(meta_label("DETECTED TYPES"))
        self.chips = ChipRow()
        types.addWidget(self.chips)
        self.summary.body.addLayout(types)

        self.summary.body.addWidget(hint("Types are editable in the next step."))

        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(theme.SPACING)
        self.change_button = secondary_button("Choose a different folder…")
        self.change_button.clicked.connect(self.choose_folder)
        row.addWidget(self.change_button)
        self.add_more_button = secondary_button("Add files…")
        self.add_more_button.clicked.connect(self.add_files)
        row.addWidget(self.add_more_button)
        row.addStretch(1)
        self.summary.body.addLayout(row)
        return self.summary

    # ---- actions ------------------------------------------------------------------------

    def _start_dir(self) -> str:
        folder = self.state.run_folder
        if folder is not None:
            return str(folder.parent if folder.parent.exists() else folder)
        return str(Path.cwd())

    def choose_folder(self) -> None:
        """Hero action: pick a run folder, read it, and move on if there was anything in it."""
        chosen = QFileDialog.getExistingDirectory(
            self, "Choose run folder", self._start_dir(), QFileDialog.Option.ShowDirsOnly
        )
        if not chosen:  # cancelled: leave the page exactly as it was
            return
        self.load_folder(Path(chosen))

    def load_folder(self, folder: Path) -> None:
        """Read ``folder`` into application state. The folder picker's whole body of work.

        Split out from :meth:`choose_folder` so the scan is reachable without a modal. The
        end-to-end smoke run drives it, and so does every test of this page.
        """
        from quant209.core.sequence import find_data_files, scan_inputs

        try:
            paths = find_data_files(folder)
        except OSError as exc:  # an unreadable share, or a folder deleted during the scan
            self.state.set_run_folder(folder)
            self.banner.show_message(
                "That folder could not be read", str(exc), tone="error"
            )
            self.refresh()
            return

        # Set the folder first: it clears any previous run (see AppState.set_run_folder), and
        # doing it after the entries would throw the entries away again.
        self.state.set_run_folder(folder)
        if not paths:
            self.banner.show_message(
                "No data files in this folder",
                f"{folder} holds no .D acquisitions and no .mzML files. Choose the folder "
                "one level up, or the Converted subfolder if this batch is already "
                "converted.",
                tone="warning",
            )
            self.refresh()
            return

        self.banner.clear()
        self.state.set_entries(scan_inputs(list(paths)))
        self.refresh()
        # The types are already assigned. The analyst's next job is to check them, so the
        # application takes them there rather than making them find the step themselves.
        self.advance_to(STEP_ASSIGN)

    def add_files(self) -> None:
        """Secondary action: add individual .mzML files to whatever is already selected."""
        chosen, _filter = QFileDialog.getOpenFileNames(
            self, "Add .mzML files", self._start_dir(), MZML_FILTER
        )
        if not chosen:
            return
        self.load_files([Path(p) for p in chosen])

    def load_files(self, paths: list[Path]) -> None:
        """Add ``paths`` to the selection, keeping the existing entries and their overrides."""
        from quant209.core.sequence import scan_inputs

        existing = self.state.entries
        if not existing and self.state.run_folder is None:
            # Give the run an identity even when it was assembled file by file: the sidebar,
            # the Convert step's output directory and the export's default name all read it.
            self.state.set_run_folder(paths[0].parent)
            existing = self.state.entries  # set_run_folder clears, so re-read rather than assume

        known = {entry.path for entry in existing}
        added = scan_inputs([p for p in paths if p not in known])
        if not added:
            self.banner.show_message(
                "Those files are already selected",
                "Nothing was added. Every file chosen is already in this run.",
                tone="info",
            )
            self.refresh()
            return

        self.banner.clear()
        self.state.set_entries(existing + added)
        self.refresh()
        self.advance_to(STEP_ASSIGN)

    # ---- rendering ----------------------------------------------------------------------

    def refresh(self) -> None:
        entries = self.state.entries
        folder = self.state.run_folder
        self.empty.setVisible(not entries)
        self.summary.setVisible(bool(entries))
        if entries:
            self.folder_label.setText(str(folder) if folder else "Files selected individually")
            self.counts_label.setText(self._counts_text(entries))
            self.chips.set_chips(role_chips(entries))
        self.status_changed.emit()

    @staticmethod
    def _counts_text(entries) -> str:
        """The two numbers Convert depends on, shown as "29 data files  ·  0 Agilent .D  ·  29
        mzML"."""
        mzml = sum(1 for e in entries if e.path.suffix.lower() == ".mzml")
        raw = len(entries) - mzml
        parts = [f"{len(entries)} data file{'' if len(entries) == 1 else 's'}"]
        if raw:
            parts.append(f"{raw} Agilent .D")
        if mzml:
            parts.append(f"{mzml} mzML")
        return "  ·  ".join(parts)

    def page_status(self) -> tuple[str, str] | None:
        entries = self.state.entries
        if not entries:
            return None
        return (f"{len(entries)} files", "accent")
