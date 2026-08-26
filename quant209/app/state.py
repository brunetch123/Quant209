"""The application's single mutable object: what run is open, how far it has got.

:class:`AppState` is the only place the GUI keeps anything, and it holds *no science*. It
knows which files were selected and what types they were given, where their converted mzML
should be, whether a :class:`~quant209.core.pipeline.RunResult` has come back, and whether
a worker is running. Every number in it came from somewhere else.

Why the workflow gate lives here and not in the sidebar
--------------------------------------------------------
The sidebar stepper and the pages both need the same answer to "how far can the user go?",
and if either computed it, the two would eventually disagree. The stepper would offer a
step the page then refused to run. :meth:`AppState.furthest_step` is therefore the single
authority, the stepper renders it, and :meth:`MainWindow.goto_step` enforces it.

The gate is a chain, evaluated from the far end back:

============  ==========================================================================
Select        always, since the user must be able to start over
Assign        once files have been selected
Convert       once files have been selected, and it may have nothing to do
Process       once every entry resolves to an mzML that exists on disk
Review        once a result has come back
Export        once a result has come back
============  ==========================================================================

Convert unlocks together with Assign rather than after it. Type assignment is auto-detected
from filenames and is *already done* when the files land, so making the user visit Assign to
acknowledge it would be a click that means nothing. A run folder that already holds converted
mzML skips Convert entirely. Process then unlocks the moment the files are selected, which is
the fast path from Select straight to Process.

Signals rather than polling
----------------------------
Each mutator emits a specific signal (``entries_changed``, ``result_changed``, ...) *and*
:attr:`availability_changed`. Pages listen to the specific one they care about. The stepper
listens only to ``availability_changed`` and re-renders. Nothing polls, and nothing has to
know which other widget also needed to hear about a change.

``msconvert_path`` persistence
-------------------------------
It is the one setting a user should never have to enter twice, so it is written to
``QSettings``. It is also auto-detected from ``MSconvert/pwiz-bin/msconvert.exe`` beside the
run folder or in the repository. A stored path that no longer exists, after a ProteoWizard
install has moved, is treated as absent so auto-detection can replace it. The alternative is
failing at conversion time with a stale path.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from PySide6.QtCore import QObject, QSettings, Signal

from quant209.core.sequence import FileEntry

if TYPE_CHECKING:  # pragma: no cover. Typing only. Importing this would pull pandas in
    from quant209.core.pipeline import RunResult

# --------------------------------------------------------------------------------------
# The six workflow steps
# --------------------------------------------------------------------------------------

STEP_SELECT = 0
STEP_ASSIGN = 1
STEP_CONVERT = 2
STEP_PROCESS = 3
STEP_REVIEW = 4
STEP_EXPORT = 5
STEP_COUNT = 6

#: Sidebar labels, one word each, so the stepper reads as a path rather than a paragraph.
STEP_LABELS: tuple[str, ...] = ("Select", "Assign", "Convert", "Process", "Review", "Export")

#: Page titles. One noun phrase each: the screen says what it is, and the README says
#: what to do with it. Explanatory copy on a working screen is read once and then
#: occupies space for the rest of the tool's life.
STEP_TITLES: tuple[str, ...] = (
    "Select a run",
    "Assign file types",
    "Convert to mzML",
    "Process the run",
    "Review the results",
    "Export the report",
)

#: Where this lab keeps ProteoWizard, relative to a run folder or the repository root.
MSCONVERT_RELATIVE = Path("MSconvert") / "pwiz-bin" / "msconvert.exe"

#: ``core.sequence``'s file type for a hexane carryover blank. Named here rather than
#: imported from ``core.pipeline`` (which would pull pandas into the startup path) and
#: pinned to it by ``tests/test_app_pages.py``.
HEX_TYPE = "hex"

#: msconvert's output directory name inside a run folder. Matches ``core.sequence.
#: find_data_files`` and the CLI's default, so a folder converted by the CLI is recognised
#: as converted by the GUI and vice versa.
CONVERTED_DIRNAME = "Converted"

_SETTINGS_MSCONVERT = "paths/msconvert"


def detect_msconvert(run_folder: Path | None = None) -> Path | None:
    """Find ``msconvert.exe`` in the usual places, nearest-first.

    Looked for beside the run folder before the repository so a lab that ships a
    ProteoWizard copy inside each run's parent directory wins over a developer checkout's.
    """
    candidates: list[Path] = []
    if run_folder is not None:
        candidates += [run_folder / MSCONVERT_RELATIVE, run_folder.parent / MSCONVERT_RELATIVE]
    # …/quant209/app/state.py -> repository root
    candidates.append(Path(__file__).resolve().parents[2] / MSCONVERT_RELATIVE)
    for candidate in candidates:
        try:
            if candidate.is_file():
                return candidate
        except OSError:  # pragma: no cover. An unreachable network path, for instance.
            continue
    return None


class AppState(QObject):
    """The open run, and how far through the workflow it has got."""

    #: ``Path | None``
    run_folder_changed = Signal(object)
    #: ``list[FileEntry]``
    entries_changed = Signal(list)
    #: ``Path | None``
    converted_dir_changed = Signal(object)
    #: ``Path | None``
    msconvert_path_changed = Signal(object)
    #: ``RunResult | None``, typed as ``object`` so this module never
    #: has to import ``core.pipeline`` (and therefore pandas) to declare its own signals.
    result_changed = Signal(object)
    #: ``Path | None``, the Excel workbook the Export step wrote.
    export_path_changed = Signal(object)
    busy_changed = Signal(bool)
    #: Emitted after *any* change that could move :meth:`furthest_step` or
    #: :meth:`completed_steps`. The stepper listens to this one signal and nothing else.
    availability_changed = Signal()

    def __init__(self, settings: QSettings | None = None, parent: QObject | None = None) -> None:
        super().__init__(parent)
        # Injectable so tests (and a future "reset all settings") can point persistence at
        # a throwaway file instead of the user's registry hive.
        self._settings = settings if settings is not None else QSettings()
        self._run_folder: Path | None = None
        self._entries: list[FileEntry] = []
        self._converted_dir: Path | None = None
        self._result: Any = None
        self._export_path: Path | None = None
        #: Steps the user has actually opened. Only Review needs it, as
        #: :meth:`completed_steps` shows, but it is recorded for every step so the rule stays
        #: one rule.
        self._visited: set[int] = set()
        self._busy = False
        #: True while the current msconvert path is a guess this module made rather than
        #: one the user or a previous session chose. See :meth:`set_msconvert_path`.
        self._msconvert_auto = True
        self._msconvert_path = self._load_msconvert_path()

    @property
    def settings(self) -> QSettings:
        """The store this state persists to. It is injected, so a test never touches the
        user's own settings.

        Exposed because more than this class needs it: the Export page remembers the output
        folder the analyst chose, and a page that reached for a bare ``QSettings()`` instead
        would write into the developer's real registry hive from inside the test suite.
        """
        return self._settings

    # ---- run folder -------------------------------------------------------------------

    @property
    def run_folder(self) -> Path | None:
        return self._run_folder

    def set_run_folder(self, folder: Path | str | None) -> None:
        """Open a run folder. A *different* folder discards the previous run entirely.

        Keeping the old entries or result around after the user points at another folder
        would let a Review page show numbers from a run whose name is no longer on screen,
        which is the worst kind of wrong: plausible.
        """
        folder = Path(folder) if folder is not None else None
        if folder == self._run_folder:
            return
        self._run_folder = folder
        self._entries = []
        self._result = None
        self._export_path = None
        self._visited.clear()
        self._converted_dir = folder / CONVERTED_DIRNAME if folder is not None else None
        # Re-detect for the new folder unless the user pinned a path themselves: a run
        # folder that ships its own ProteoWizard is the nearest and most likely right one,
        # but a deliberate choice outranks any guess.
        if self._msconvert_auto or self._msconvert_path is None or not self._msconvert_path.is_file():
            detected = detect_msconvert(folder)
            if detected is not None:
                self._apply_msconvert_path(detected, auto=True)
        self.run_folder_changed.emit(self._run_folder)
        self.entries_changed.emit(self._entries)
        self.result_changed.emit(None)
        self.export_path_changed.emit(None)
        self.converted_dir_changed.emit(self._converted_dir)
        self.availability_changed.emit()

    # ---- entries ----------------------------------------------------------------------

    @property
    def entries(self) -> list[FileEntry]:
        """The run's files. A copy: callers mutate their own list, not the state's."""
        return list(self._entries)

    def set_entries(self, entries: list[FileEntry]) -> None:
        self._entries = list(entries)
        self.entries_changed.emit(self.entries)
        self.availability_changed.emit()

    def set_entry_type(self, index: int, file_type: str) -> None:
        """Override one file's detected type (the Assign step's whole job).

        Sets ``auto=False`` on the entry so that downstream code, and the Assign table's own
        "auto" badge, can tell a human decision from a filename guess.
        """
        entry = self._entries[index]
        if entry.file_type == file_type:
            return
        entry.file_type = file_type
        entry.auto = False
        self.entries_changed.emit(self.entries)
        self.availability_changed.emit()

    # ---- converted mzML ---------------------------------------------------------------

    @property
    def converted_dir(self) -> Path | None:
        """Where msconvert output lives: ``<run folder>/Converted`` unless overridden."""
        return self._converted_dir

    def set_converted_dir(self, folder: Path | str | None) -> None:
        folder = Path(folder) if folder is not None else None
        if folder == self._converted_dir:
            return
        self._converted_dir = folder
        self.converted_dir_changed.emit(folder)
        self.availability_changed.emit()

    def mzml_for(self, entry: FileEntry) -> Path | None:
        """The mzML this entry would be read from, or ``None`` if there is nowhere to look.

        Mirrors ``core.pipeline._mzml_path``. An entry that is already an mzML is itself, and
        anything else is ``<converted_dir>/<stem>.mzML``. Existence is *not* checked here,
        :meth:`mzml_ready` is where that question is asked.
        """
        if entry.path.suffix.lower() == ".mzml":
            return entry.path
        if self._converted_dir is None:
            return None
        return self._converted_dir / f"{entry.path.stem}.mzML"

    def mzml_ready(self) -> bool:
        """Whether every file the pipeline will *read* has an mzML on disk.

        Hexane carryover blanks are excluded, because ``core.pipeline.QUANTIFIED_TYPES`` does
        not include them: the pipeline never opens a ``hex`` file, so an unconverted one
        cannot stop a run, and blocking Process on nineteen conversions whose output nothing
        will read would be the application inventing work. They are still *offered* for
        conversion, as :meth:`pending_conversions` shows, because the analyst selected them.
        """
        required = [e for e in self._entries if e.file_type != HEX_TYPE]
        if not required:
            return False
        for entry in required:
            mzml = self.mzml_for(entry)
            if mzml is None or not mzml.is_file():
                return False
        return True

    def pending_conversions(self) -> list[FileEntry]:
        """The entries with no mzML yet, which is what the Convert step still owes."""
        pending = []
        for entry in self._entries:
            mzml = self.mzml_for(entry)
            if mzml is None or not mzml.is_file():
                pending.append(entry)
        return pending

    # ---- msconvert --------------------------------------------------------------------

    @property
    def msconvert_path(self) -> Path | None:
        return self._msconvert_path

    @property
    def msconvert_auto(self) -> bool:
        """Whether the current path is a guess this module made rather than a user's choice.

        The Convert page says so on screen ("found automatically at ..."), because a path the
        analyst never typed is one they have no reason to trust until they are told where it
        came from.
        """
        return self._msconvert_auto

    def set_msconvert_path(self, path: Path | str | None) -> None:
        """Pin the msconvert executable. Persisted, and never overridden by detection."""
        self._apply_msconvert_path(Path(path) if path is not None else None, auto=False)

    def _apply_msconvert_path(self, path: Path | None, *, auto: bool) -> None:
        self._msconvert_auto = auto
        if path == self._msconvert_path:
            return
        self._msconvert_path = path
        # Only a deliberate choice is written to settings. Persisting a guess would turn it
        # into a pinned path on the next launch, and detection would stop happening at all.
        if not auto:
            if path is None:
                self._settings.remove(_SETTINGS_MSCONVERT)
            else:
                self._settings.setValue(_SETTINGS_MSCONVERT, str(path))
            self._settings.sync()
        self.msconvert_path_changed.emit(path)
        self.availability_changed.emit()

    def _load_msconvert_path(self) -> Path | None:
        stored = self._settings.value(_SETTINGS_MSCONVERT, None)
        if stored:
            path = Path(str(stored))
            if path.is_file():
                self._msconvert_auto = False
                return path
        self._msconvert_auto = True
        return detect_msconvert(None)

    # ---- result -----------------------------------------------------------------------

    @property
    def result(self) -> RunResult | None:
        return self._result

    def set_result(self, result: RunResult | None) -> None:
        self._result = result
        # A fresh result un-reviews the run: whatever the analyst read last time was about
        # different numbers.
        self._visited.discard(STEP_REVIEW)
        if result is None:
            self._export_path = None
            self.export_path_changed.emit(None)
        self.result_changed.emit(result)
        self.availability_changed.emit()

    @property
    def export_path(self) -> Path | None:
        return self._export_path

    def set_export_path(self, path: Path | str | None) -> None:
        path = Path(path) if path is not None else None
        self._export_path = path
        self.export_path_changed.emit(path)
        self.availability_changed.emit()

    # ---- busy -------------------------------------------------------------------------

    @property
    def busy(self) -> bool:
        """True while a worker is running. Pages disable their action buttons on this."""
        return self._busy

    def set_busy(self, busy: bool) -> None:
        busy = bool(busy)
        if busy == self._busy:
            return
        self._busy = busy
        self.busy_changed.emit(busy)

    # ---- workflow gate ----------------------------------------------------------------

    def furthest_step(self) -> int:
        """The highest step the user is allowed to open. See the module docstring."""
        if self._result is not None:
            return STEP_EXPORT
        if self.mzml_ready():
            return STEP_PROCESS
        if self._entries:
            return STEP_CONVERT
        return STEP_SELECT

    def is_step_available(self, step: int) -> bool:
        return 0 <= step <= self.furthest_step()

    def completed_steps(self) -> set[int]:
        """The steps whose output now exists.

        Distinct from availability: Convert is *available* as soon as files are selected,
        and *complete* only once every one of them has an mzML. A step is marked complete
        when the thing it produces is on disk or in memory, never merely because the user
        walked past it.
        """
        done: set[int] = set()
        if self._entries:
            done.add(STEP_SELECT)
            done.add(STEP_ASSIGN)  # types are auto-assigned the moment files are selected
        if self.mzml_ready():
            done.add(STEP_CONVERT)
        if self._result is not None:
            done.add(STEP_PROCESS)
            # Review is the one step whose output is a human judgement rather than data, so
            # it completes on having been *opened*, not on a result existing. A check mark
            # against a screen nobody has looked at is a lie, and on a scientific product
            # "reviewed" is precisely the claim that must not be made on the analyst's
            # behalf.
            if STEP_REVIEW in self._visited:
                done.add(STEP_REVIEW)
        if self._export_path is not None:
            done.add(STEP_EXPORT)
        return done

    def mark_visited(self, step: int) -> None:
        """Record that the user opened ``step``. Called by the window on navigation."""
        if step in self._visited:
            return
        self._visited.add(step)
        self.availability_changed.emit()

    def refresh(self) -> None:
        """Re-evaluate everything that depends on the filesystem and notify listeners.

        :meth:`mzml_ready` reads the disk, so a conversion finishing (or a folder being
        populated behind the app's back) changes the answer without any setter being
        called. The Convert worker calls this when it finishes, and so do the tests.
        """
        self.availability_changed.emit()

    def reset(self) -> None:
        """Back to a blank application, keeping only the persisted msconvert path."""
        self._run_folder = None
        self._entries = []
        self._converted_dir = None
        self._result = None
        self._export_path = None
        self._visited.clear()
        self.set_busy(False)
        self.run_folder_changed.emit(None)
        self.entries_changed.emit([])
        self.converted_dir_changed.emit(None)
        self.result_changed.emit(None)
        self.export_path_changed.emit(None)
        self.availability_changed.emit()
