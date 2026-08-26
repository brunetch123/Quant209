"""Quant 209 desktop entry point: the window, and the cold-start discipline that shapes it.

The three-second rule
----------------------
The budget is **under three seconds from launch to a usable window**, and it is spent
almost entirely on Qt itself. Nothing else may be spent before the window paints.
Concretely, and enforced by ``tests/test_app_shell.py``:

* This module imports **no numpy, no pandas, no pyqtgraph, no openpyxl and no
  ``core.pipeline``** at module scope, directly or transitively.
* ``quant209.app.state`` and ``quant209.app.sidebar`` are cheap by construction
  (``core.sequence`` is dataclasses and ``re``).
* The workflow pages are **built lazily**, one per navigation, so a launch pays for the
  Select page and nothing else. See :mod:`quant209.app.workflow`.
* The heavy imports happen inside a worker's ``run`` (:mod:`quant209.app.workers`) or a
  page's ``on_activated``, both of which are minutes after the window is on screen.

``--measure-startup`` makes the budget checkable. The process quits as soon as the first
frame has been painted and prints the elapsed milliseconds. ``tools/measure_startup.py``
drives it and reports the wall-clock median, which includes the interpreter start that this
module cannot see.

Navigation
-----------
The window is the only thing that navigates, and it asks
:meth:`quant209.app.state.AppState.is_step_available` first, every time. The stepper's
disabled rows are the *visible* half of that gate, not the gate itself. Back-navigation is
always allowed because a completed step is always re-openable: an analyst who wants to check
what they assigned in step 2 should not have to start over.
"""

from __future__ import annotations

import sys
import time

_IMPORT_T0 = time.perf_counter()  # before Qt: the earliest point this module can observe

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QIcon
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import (
    QApplication,
    QFrame,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

import quant209
from quant209.app import theme
from quant209.app.sidebar import Sidebar
from quant209.app.state import STEP_SELECT, AppState
from quant209.app.workflow import load_page_class
from quant209.app.workflow.base import StatusChip, WorkflowPage

APP_NAME = "Quant 209"
ORG_NAME = "Quant209"

USAGE = """\
Quant 209. PCB congener quantification from Agilent GC-MS/MS runs.

usage: quant209-gui [--measure-startup] [--version] [--help]

  --measure-startup  quit as soon as the window has painted and print the elapsed
                     milliseconds (the cold-start budget check; see tools/measure_startup.py)
  --version          print the version and exit
  --help             print this message and exit

For headless/scripted runs use the `quant209` command instead.\
"""

#: Watchdog for ``--measure-startup``: if no paint event ever arrives (a platform plugin
#: with no compositor, say), quit anyway rather than leaving a process behind for the
#: measurement harness to time out on.
_MEASURE_WATCHDOG_MS = 20_000


# --------------------------------------------------------------------------------------
# Header strip
# --------------------------------------------------------------------------------------


class PageHeader(QWidget):
    """Page title, and an optional right-aligned status chip.

    Owned by the window rather than by each page: one header, rendered once, means every
    screen's title sits on the same baseline at the same size. The pages supply the words.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("PageHeader")
        layout = QHBoxLayout(self)
        layout.setContentsMargins(
            theme.CONTENT_MARGIN, theme.CONTENT_MARGIN, theme.CONTENT_MARGIN, theme.SPACING * 2
        )
        layout.setSpacing(theme.SPACING * 2)

        text = QVBoxLayout()
        text.setContentsMargins(0, 0, 0, 0)
        text.setSpacing(theme.SPACING // 2)
        self.title_label = QLabel("")
        self.title_label.setObjectName("PageTitle")
        text.addWidget(self.title_label)
        layout.addLayout(text, 1)

        self.chip = StatusChip("")
        layout.addWidget(self.chip, 0, Qt.AlignmentFlag.AlignTop)

    def show_page(self, title: str, status: tuple[str, str] | None) -> None:
        self.title_label.setText(title)
        self.chip.set_status(*(status or ("", "neutral")))


# --------------------------------------------------------------------------------------
# The window
# --------------------------------------------------------------------------------------


class MainWindow(QMainWindow):
    """Sidebar stepper on the left, header + page stack on the right."""

    #: Emitted once, just after the first frame reaches the screen. The startup measurement
    #: hangs off this. Nothing else should.
    first_painted = Signal()

    def __init__(self, state: AppState | None = None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.state = state if state is not None else AppState(parent=self)
        self._pages: dict[int, WorkflowPage] = {}
        self._current_step = STEP_SELECT
        self._painted = False

        self.setWindowTitle(APP_NAME)
        self.setMinimumSize(theme.MIN_WINDOW_WIDTH, theme.MIN_WINDOW_HEIGHT)
        self.resize(theme.MIN_WINDOW_WIDTH, theme.MIN_WINDOW_HEIGHT)

        central = QWidget()
        central.setObjectName("Central")
        root = QHBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        self.sidebar = Sidebar()
        self.sidebar.step_activated.connect(self.goto_step)
        root.addWidget(self.sidebar)

        content = QWidget()
        content.setObjectName("ContentArea")
        content_layout = QVBoxLayout(content)
        content_layout.setContentsMargins(0, 0, 0, 0)
        content_layout.setSpacing(0)

        self.header = PageHeader()
        content_layout.addWidget(self.header)

        rule = QFrame()
        rule.setObjectName("HeaderRule")
        rule.setFrameShape(QFrame.Shape.HLine)
        content_layout.addWidget(rule)

        self.stack = QStackedWidget()
        content_layout.addWidget(self.stack, 1)
        root.addWidget(content, 1)

        self.setCentralWidget(central)

        # Anything that could move the workflow gate re-renders the stepper and the
        # visible page. The specific signals drive the parts that only care about one fact.
        self.state.availability_changed.connect(self._sync)
        self.state.run_folder_changed.connect(self.sidebar.set_run_folder)
        self.state.busy_changed.connect(lambda _busy: self._sync())

        # The Select page is the one page a launch pays for.
        self.goto_step(STEP_SELECT, force=True)
        self._center_on_screen()

    # ---- pages ------------------------------------------------------------------------

    @property
    def stepper(self):
        return self.sidebar.stepper

    @property
    def current_step(self) -> int:
        return self._current_step

    def page(self, step: int) -> WorkflowPage:
        """The page for ``step``, importing and constructing it on first use."""
        page = self._pages.get(step)
        if page is None:
            page = load_page_class(step)(self.state)
            page.status_changed.connect(self._refresh_header)
            # A page asks to navigate. The window decides. ``goto_step`` re-checks the gate,
            # so an auto-advance that arrives a moment after the state it was based on
            # changed (a worker finishing into a run the user has already replaced) is
            # refused rather than showing a page built for data that is gone.
            page.navigate_requested.connect(self.goto_step)
            self._pages[step] = page
            self.stack.addWidget(page)
        return page

    def built_steps(self) -> list[int]:
        """Which pages have actually been constructed. The lazy-import test reads this."""
        return sorted(self._pages)

    # ---- navigation -------------------------------------------------------------------

    def goto_step(self, step: int, force: bool = False) -> bool:
        """Show ``step``. Returns False, changing nothing, if it is still locked.

        ``force`` is for the initial navigation only, before any state exists.
        """
        if not force and not self.state.is_step_available(step):
            return False
        page = self.page(step)
        self._current_step = step
        self.state.mark_visited(step)
        self.stack.setCurrentWidget(page)
        page.on_activated()
        page.refresh()
        self._sync()
        return True

    # ---- rendering --------------------------------------------------------------------

    def _sync(self) -> None:
        """Re-render the chrome from application state."""
        furthest = self.state.furthest_step()
        if self._current_step > furthest:
            # The gate can close while this page is on screen, if the user opens a different
            # run folder while standing on Review. Fall back rather than leave a page showing
            # a run that no longer exists.
            self.goto_step(furthest)
            return
        self.stepper.update_states(self._current_step, furthest, self.state.completed_steps())
        self.stepper.follow_focus(self._current_step)
        self.sidebar.set_run_folder(self.state.run_folder)
        page = self._pages.get(self._current_step)
        if page is not None:
            page.refresh()
        self._refresh_header()

    def _refresh_header(self) -> None:
        page = self._pages.get(self._current_step)
        if page is None:  # pragma: no cover (the current page is always built)
            return
        self.header.show_page(page.title, page.page_status())

    def _center_on_screen(self) -> None:
        screen = QGuiApplication.primaryScreen()
        if screen is None:  # pragma: no cover (headless without even a virtual screen)
            return
        available = screen.availableGeometry()
        frame = self.frameGeometry()
        frame.moveCenter(available.center())
        self.move(frame.topLeft())

    # ---- first paint ------------------------------------------------------------------

    def paintEvent(self, event) -> None:  # noqa: N802 (Qt's spelling)
        super().paintEvent(event)
        if not self._painted:
            self._painted = True
            # Deferred by one event-loop turn so the children have painted too: the signal
            # means "there is a finished frame on screen", not "the window began painting".
            QTimer.singleShot(0, self.first_painted.emit)


# --------------------------------------------------------------------------------------
# Application entry point
# --------------------------------------------------------------------------------------


def configure_app(app: QApplication) -> None:
    """Name the application (QSettings reads this) and apply the design system."""
    app.setApplicationName(APP_NAME)
    app.setOrganizationName(ORG_NAME)
    app.setApplicationDisplayName(APP_NAME)
    app.setApplicationVersion(quant209.__version__)
    icon = theme.ASSETS / "quant209.ico"
    if icon.is_file():
        app.setWindowIcon(QIcon(str(icon)))
    theme.apply_theme(app)


def _write(line: str) -> None:
    """Print, tolerating the GUI-subsystem launcher's absent stdout.

    ``quant209-gui`` is a ``gui_scripts`` entry point, so it runs under ``pythonw.exe``:
    with no console and no redirection, ``sys.stdout`` is ``None`` and a bare ``print``
    raises. Under ``tools/measure_startup.py``'s pipe it is a real stream, which is the
    case that matters.
    """
    stream = sys.stdout
    if stream is None:  # pragma: no cover (only reachable under pythonw with no pipe)
        return
    stream.write(line + "\n")
    stream.flush()


def _install_startup_measurement(app: QApplication, window: MainWindow) -> None:
    def report() -> None:
        elapsed_ms = (time.perf_counter() - _IMPORT_T0) * 1000.0
        _write(f"startup_ms={elapsed_ms:.1f}")
        app.quit()

    window.first_painted.connect(report)
    QTimer.singleShot(_MEASURE_WATCHDOG_MS, app.quit)


def main(argv: list[str] | None = None) -> int:
    """Launch the application. Returns the process exit code."""
    args = list(sys.argv[1:] if argv is None else argv)
    if "--help" in args or "-h" in args:
        _write(USAGE)
        return 0
    if "--version" in args:
        _write(f"{APP_NAME} {quant209.__version__}")
        return 0
    measure = "--measure-startup" in args

    # Fractional display scaling (125% is a common Windows setting) rounds to 100% or 200%
    # by default, which is why so many Qt apps look either cramped or enormous there. Must
    # be set before the QApplication exists.
    QApplication.setHighDpiScaleFactorRoundingPolicy(
        Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
    )
    app = QApplication.instance() or QApplication(sys.argv[:1])
    configure_app(app)

    window = MainWindow()
    if measure:
        _install_startup_measurement(app, window)
    window.show()
    return app.exec()


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
