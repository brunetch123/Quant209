"""The left rail: the mark and wordmark, the six-step stepper, and the run indicator.

The stepper is the only navigation in the application. It is always on screen, so at any
moment the analyst can see the whole job, which part of it they are in, and what is left.
Four states, each visually distinct at a glance:

=========  ====================================================================
locked     not reachable yet. The row is disabled, so a click does nothing at
           the widget level rather than being swallowed by a handler
available  reachable. Hover and focus respond
active     the page currently on screen: accent fill, accent bold label
complete   its output exists. The number circle becomes a green check
=========  ====================================================================

Why the rows are ``QPushButton``\\ s
-------------------------------------
A ``QWidget`` with a ``mousePressEvent`` would need a hand-written ``paintEvent`` before a
stylesheet would paint its background at all, plus hand-written hover tracking, plus its own
disabled look, and it would still not be reachable by keyboard. ``QPushButton`` brings all
four for free: ``setEnabled(False)`` *is* the locked state (Qt refuses the click and the
``:disabled`` rules restyle it), ``:hover``/``:pressed``/``:focus`` are stylesheet states,
and Tab-then-Space navigates the workflow without a mouse. The number circle and the label
are child ``QLabel``\\ s inside the button's own layout, which is what allows the circle to
be a circle and the label to change weight independently.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSize, Qt, Signal
from PySide6.QtGui import QIcon, QPixmap
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

import quant209
from quant209.app import theme
from quant209.app.state import STEP_COUNT, STEP_LABELS

#: The four row states, in increasing precedence when several apply at once.
STATE_LOCKED = "locked"
STATE_AVAILABLE = "available"
STATE_COMPLETE = "complete"
STATE_ACTIVE = "active"

#: U+2713 CHECK MARK: present in Segoe UI, Segoe UI Variable and the DejaVu fallback the
#: offscreen platform uses, so the complete state renders everywhere the app runs.
CHECK_GLYPH = "✓"


class StepRow(QPushButton):
    """One step in the stepper: number circle, label, and a state."""

    def __init__(self, index: int, label: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("StepRow")
        self.index = index
        self.step_label = label
        self.setFlat(True)
        # The row's height is set by ``min-height`` in the stylesheet, not here. See the
        # note above ``QPushButton#StepRow`` in theme.py.
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

        row = QHBoxLayout(self)
        row.setContentsMargins(theme.SPACING, theme.SPACING, theme.SPACING, theme.SPACING)
        row.setSpacing(theme.SPACING + 4)

        self.number_label = QLabel(str(index + 1))
        self.number_label.setObjectName("StepNumber")
        self.number_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        row.addWidget(self.number_label)

        self.text_label = QLabel(label)
        self.text_label.setObjectName("StepLabel")
        row.addWidget(self.text_label)
        row.addStretch(1)

        self._state = ""
        # repolish=False: the widget has not been polished yet, so Qt will read these
        # properties when it first applies the stylesheet. Re-polishing six rows (three
        # widgets each) against a full stylesheet before any of them has been shown is pure
        # startup cost, and the cold-start budget is not big enough to waste it.
        self.set_state(STATE_LOCKED, repolish=False)

    @property
    def state(self) -> str:
        """One of ``locked`` / ``available`` / ``active`` / ``complete``."""
        return self._state

    def set_state(self, state: str, repolish: bool = True) -> None:
        if state == self._state:
            return
        self._state = state
        # The check replaces the number only once the step's output exists, so a glance
        # down the rail counts what is done without reading a word.
        self.number_label.setText(CHECK_GLYPH if state == STATE_COMPLETE else str(self.index + 1))
        self.setEnabled(state != STATE_LOCKED)
        self.setCursor(
            Qt.CursorShape.ArrowCursor if state == STATE_LOCKED
            else Qt.CursorShape.PointingHandCursor
        )
        for widget in (self, self.number_label, self.text_label):
            widget.setProperty("state", state)
        if repolish:
            theme.repolish(self, self.number_label, self.text_label)


class Stepper(QWidget):
    """The six :class:`StepRow`\\ s, and the mapping from application state to their states."""

    #: Emitted with the step index when an unlocked row is clicked.
    step_activated = Signal(int)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("Stepper")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)

        self.steps: list[StepRow] = []
        for index in range(STEP_COUNT):
            row = StepRow(index, STEP_LABELS[index])
            row.clicked.connect(lambda _checked=False, i=index: self.step_activated.emit(i))
            layout.addWidget(row)
            self.steps.append(row)

    def state_of(self, step: int) -> str:
        return self.steps[step].state

    def follow_focus(self, step: int) -> None:
        """Move the keyboard focus to ``step``, but only if it is already in the stepper.

        Two rows must never look highlighted at once. The focus ring and the active fill are
        different signals ("Space activates this" vs "you are here") and they normally sit on
        the same row, because navigating means clicking it. They come apart when the app
        navigates *for* the user (Process auto-advancing to Review on completion), which
        leaves the ring behind on the row they last touched.

        Conditional, because the alternative is worse: unconditionally focusing the new step
        would yank focus out of whatever the user was typing in on the page every time state
        advanced.
        """
        focused = self.window().focusWidget() if self.window() else None
        if focused in self.steps and focused is not self.steps[step]:
            self.steps[step].setFocus(Qt.FocusReason.OtherFocusReason)

    def update_states(self, current: int, furthest: int, completed: set[int]) -> None:
        """Re-render from the three facts that decide every row.

        Precedence is active > complete > available > locked: the page you are looking at
        reads as *here* even when its output already exists, because "where am I" is the
        question the stepper is on screen to answer.
        """
        for index, row in enumerate(self.steps):
            if index == current:
                state = STATE_ACTIVE
            elif index > furthest:
                state = STATE_LOCKED
            elif index in completed:
                state = STATE_COMPLETE
            else:
                state = STATE_AVAILABLE
            row.set_state(state)


def mark_pixmap(width: int) -> QPixmap:
    """The PCB 77 mark, rasterised to ``width`` pixels.

    ``QIcon`` renders an SVG at whatever size is asked for, so the mark stays sharp on a
    high-density display instead of being scaled up from a fixed bitmap.
    """
    source = QIcon(str(theme.ASSETS / "mark.svg"))
    # Asking for a square fits the drawing inside it at its own aspect ratio, so the height
    # comes from the artwork rather than from a number that could drift away from it.
    return source.pixmap(QSize(width, width))


class Sidebar(QWidget):
    """Mark, wordmark, stepper, and once a run is open, which folder is being worked on."""

    step_activated = Signal(int)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("Sidebar")
        self.setFixedWidth(theme.SIDEBAR_WIDTH)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(
            theme.SPACING * 2, theme.SPACING * 3, theme.SPACING * 2, theme.SPACING * 2
        )
        layout.setSpacing(theme.SPACING)

        self.mark = QLabel()
        self.mark.setObjectName("Mark")
        self.mark.setPixmap(mark_pixmap(theme.MARK_WIDTH))
        self.mark.setAccessibleName("Quant 209")
        layout.addWidget(self.mark)

        self.wordmark = QLabel("Quant 209")
        self.wordmark.setObjectName("Wordmark")
        layout.addWidget(self.wordmark)

        self.wordmark_tag = QLabel("PCB congener quantification")
        self.wordmark_tag.setObjectName("WordmarkTag")
        self.wordmark_tag.setWordWrap(True)
        layout.addWidget(self.wordmark_tag)

        layout.addSpacing(theme.SPACING * 3)

        self.stepper = Stepper()
        self.stepper.step_activated.connect(self.step_activated)
        layout.addWidget(self.stepper)

        layout.addStretch(1)

        self.run_folder_label = QLabel("")
        self.run_folder_label.setObjectName("RunFolderLabel")
        self.run_folder_label.setWordWrap(True)
        self.run_folder_label.setVisible(False)
        layout.addWidget(self.run_folder_label)

        divider = QFrame()
        divider.setObjectName("SidebarDivider")
        divider.setFrameShape(QFrame.Shape.HLine)
        layout.addWidget(divider)

        self.version_label = QLabel(f"Version {quant209.__version__}")
        self.version_label.setObjectName("SidebarFooter")
        layout.addWidget(self.version_label)

    def set_run_folder(self, folder: Path | None) -> None:
        """Name the open run at the foot of the rail. The full path is the tooltip.

        The folder name is the run's identity in the lab (a batch, a site and a date), so it
        is what belongs on screen. The full path is what the analyst needs when two runs share
        a name, and a tooltip is the right weight for it.
        """
        if folder is None:
            self.run_folder_label.setText("")
            self.run_folder_label.setToolTip("")
            self.run_folder_label.setVisible(False)
            return
        self.run_folder_label.setText(f"Run  ·  {folder.name}")
        self.run_folder_label.setToolTip(str(folder))
        self.run_folder_label.setVisible(True)
