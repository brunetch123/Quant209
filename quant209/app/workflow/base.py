"""The scaffold every workflow page is built on, plus the shared content widgets.

One page, one question
-----------------------
A :class:`WorkflowPage` owns a *step* of the six-step workflow and answers exactly one
question for the analyst. It declares three things about itself, :attr:`~WorkflowPage.step`
and :attr:`~WorkflowPage.title`, which the window renders in
the header strip above the page rather than the page drawing its own. That is what keeps
every screen's title in the same place, at the same size, with the same baseline. A header
each page drew for itself would drift within two pages.

The layout contract
--------------------
Subclasses never set their own margins. They add widgets to :attr:`~WorkflowPage.content`,
a vertical layout already inset by ``theme.CONTENT_MARGIN`` and already capped at
``theme.CONTENT_MAX_WIDTH`` and left-aligned, so a page is legible on a 3840-px monitor
instead of stretching one table across a metre of glass. A page that genuinely needs the
full width (the chromatogram browser) sets :attr:`max_content_width` to ``0``, and one whose
content grows past the window height sets :attr:`~WorkflowPage.scrollable`. Both are opt-in,
so a page that says nothing gets the capped, non-scrolling column every screen shares.

Lifecycle
----------
* :meth:`on_activated`: the page just became visible. **This is where a page does its
  expensive imports**: the cold-start budget is spent before the window paints, so pyqtgraph
  and pandas may only be reached for once the user has navigated somewhere that needs them.
* :meth:`refresh`: application state changed while this page was visible.

Both are no-ops by default, so a page implements only the one it needs.

The shared vocabulary
----------------------
Everything below the scaffold is here because more than one page needs it and two pages
inventing it separately would drift: the five file roles and their human names
(:data:`ROLE_LABELS` and friends), the :class:`StatusChip` / :class:`ChipRow` counters, the
in-page :class:`Banner`, the :class:`Card` and :class:`EmptyState` containers, and the three
button ranks. Adding a widget here is the right move exactly when a second page wants it,
not before.

Navigation is a *request*
--------------------------
A page never shows another page. It emits :attr:`~WorkflowPage.navigate_requested` (via
:meth:`~WorkflowPage.advance_to`) and the window decides, re-checking the workflow gate on
the way. That is what lets Select hand off to Assign, and Process to Review, without either
page having to know whether the step it is pointing at is even unlocked.
"""

from __future__ import annotations

from PySide6.QtCore import QSize, Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QLayout,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from quant209.app import theme
from quant209.app.state import STEP_TITLES, AppState

# --------------------------------------------------------------------------------------
# The file-role vocabulary, shared by every page that names one
# --------------------------------------------------------------------------------------
# ``core.sequence`` types are three-to-six-letter tokens meant for code (``cal``, ``ref``,
# ``hex``). Four pages show them to a human, and four independent translations of the same
# five words is precisely how an application ends up calling the same thing a "hexane blank"
# on one screen and a "hex" on the next.

#: The five file types, in the order an analyst thinks about them: what calibrates the run,
#: what measures recovery, what is being measured, and the two kinds of blank.
ROLE_ORDER: tuple[str, ...] = ("cal", "ref", "sample", "blank", "hex")

#: Singular, for a table cell or a dropdown item.
ROLE_LABELS: dict[str, str] = {
    "cal": "calibration",
    "ref": "reference",
    "sample": "sample",
    "blank": "blank",
    "hex": "hexane blank",
}

#: Counted, for a chip or a sentence: ``ROLE_COUNTED["sample"] % 15`` -> "15 samples".
ROLE_COUNTED: dict[str, str] = {
    "cal": "%d calibration",
    "ref": "%d reference",
    "sample": "%d samples",
    "blank": "%d blanks",
    "hex": "%d hexane",
}

#: Chip tone per role. The two that decide whether a run can be quantified at all carry the
#: accent. The rest stay quiet, so a glance at the chip row finds the calibration first.
ROLE_TONES: dict[str, str] = {
    "cal": "accent",
    "ref": "accent",
    "sample": "neutral",
    "blank": "info",
    "hex": "neutral",
}


def role_chips(entries) -> list[tuple[str, str]]:
    """``(text, tone)`` chips counting ``entries`` by role, in :data:`ROLE_ORDER`.

    Roles with no files are left out rather than shown as zero: "0 blanks" is a fact about
    nothing, and the pages that care say so in a banner instead.
    """
    counts: dict[str, int] = {}
    for entry in entries:
        counts[entry.file_type] = counts.get(entry.file_type, 0) + 1
    chips = [
        (ROLE_COUNTED[role] % counts[role], ROLE_TONES[role])
        for role in ROLE_ORDER
        if counts.get(role)
    ]
    # Anything core.sequence starts detecting that this vocabulary has not caught up with
    # still appears, rather than silently vanishing from the count.
    chips += [
        (f"{count} {role}", "neutral")
        for role, count in sorted(counts.items())
        if role not in ROLE_LABELS
    ]
    return chips


# --------------------------------------------------------------------------------------
# Small shared widgets
# --------------------------------------------------------------------------------------


class StatusChip(QLabel):
    """A pill carrying one short verdict: "12 files", "Ready to process", "3 flags".

    Hidden when its text is empty, so a page can leave the chip permanently in the header
    and simply not have anything to say yet. An empty pill outline is visual noise.
    """

    def __init__(self, text: str = "", tone: str = "neutral", parent: QWidget | None = None) -> None:
        super().__init__(text, parent)
        self.setObjectName("StatusChip")
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed)
        self._tone = ""
        # repolish=False: this widget has never been polished, so Qt will read the ``tone``
        # property when it first applies the stylesheet. Re-polishing against a
        # stylesheet this size costs ~65 ms per widget on the target host, measured. A chip
        # row is rebuilt on every refresh, so paying it here turned a type override into a
        # third of a second of frozen window. Same defect (and same fix) as ``StepRow``.
        self.set_status(text, tone, repolish=False)

    def set_status(self, text: str, tone: str = "neutral", repolish: bool = True) -> None:
        changed = text != self.text()
        self.setText(text)
        self.setVisible(bool(text))
        if changed and self.parent() is not None:
            # The header chip's text changes on every state transition (Ready -> Running ->
            # Complete). ``setText`` only *requests* a re-layout, so for one frame the new,
            # longer word is painted inside the old, shorter pill and reads as clipped.
            self.adjustSize()
        if tone != self._tone:
            self._tone = tone
            self.setProperty("tone", tone)
            if repolish:
                theme.repolish(self)

    @property
    def tone(self) -> str:
        return self._tone


class ChipRow(QWidget):
    """A left-aligned row of :class:`StatusChip`\\ s, rebuilt from a list of ``(text, tone)``.

    Rebuilt rather than mutated because the *number* of chips changes with the data (a run
    with no blanks has no blank chip), and a row of empty pills held in reserve is exactly the
    kind of skeleton that reads as unfinished software.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("ChipRow")
        self._layout = QHBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)
        self._layout.setSpacing(theme.SPACING)
        self._layout.addStretch(1)
        self.chips: list[StatusChip] = []

    def set_chips(self, chips: list[tuple[str, str]]) -> None:
        # Pages refresh on every state change, and most of those changes leave the counts
        # exactly as they were. Rebuilding identical widgets would churn the layout (and the
        # focus chain) for no visible difference.
        if chips == [(chip.text(), chip.tone) for chip in self.chips]:
            return
        for chip in self.chips:
            self._layout.removeWidget(chip)
            chip.deleteLater()
        self.chips = []
        for index, (text, tone) in enumerate(chips):
            # Parented at construction, not on insertion: ``StatusChip`` shows itself, and a
            # parentless widget that is shown becomes a *top-level window*. Qt then creates a
            # native window for it and destroys it again one line later, which measured ~45 ms
            # per chip. On a real display, it also flickers.
            chip = StatusChip(text, tone, self)
            self._layout.insertWidget(index, chip)
            self.chips.append(chip)
        self.setVisible(bool(chips))

    def texts(self) -> list[str]:
        """What the row currently says. The tests read this rather than pixel positions."""
        return [chip.text() for chip in self.chips]


#: The glyph in a banner's coloured disc, per tone. Plain ASCII and U+2713, both of which the
#: platform UI font and the offscreen fallback have. A warning triangle or an emoji
#: would render from Segoe UI Emoji on Windows, in full colour, at the wrong size.
BANNER_ICONS: dict[str, str] = {
    "error": "!",
    "warning": "!",
    "info": "i",
    "success": "✓",
    "neutral": "i",
}


class Banner(QFrame):
    """One in-page message: a coloured disc, a headline, and an optional second line.

    In-page and not a ``QMessageBox``, always. A modal stops the application to say something
    the page had room for, forces an acknowledgement the analyst cannot re-read, and steals
    the keyboard. On a screen whose job is "here is what is wrong with this run", the
    message *is* the content.

    Hidden when it has nothing to say, so a page can own three banners permanently and show
    whichever apply.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("Banner")
        self._tone = ""
        row = QHBoxLayout(self)
        pad = theme.SPACING + 4
        row.setContentsMargins(pad, pad, pad, pad)
        row.setSpacing(pad)

        self.icon_label = QLabel("")
        self.icon_label.setObjectName("BannerIcon")
        self.icon_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.icon_label.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        row.addWidget(self.icon_label, 0, Qt.AlignmentFlag.AlignTop)

        text = QVBoxLayout()
        text.setContentsMargins(0, 0, 0, 0)
        text.setSpacing(2)
        self.title_label = QLabel("")
        self.title_label.setObjectName("BannerTitle")
        self.title_label.setWordWrap(True)
        text.addWidget(self.title_label)
        self.body_label = QLabel("")
        self.body_label.setObjectName("BannerBody")
        self.body_label.setWordWrap(True)
        self.body_label.setVisible(False)
        text.addWidget(self.body_label)
        row.addLayout(text, 1)

        #: Where a page can park one trailing action ("Show details", "Open folder").
        self.action_slot = QHBoxLayout()
        self.action_slot.setContentsMargins(0, 0, 0, 0)
        self.action_slot.setSpacing(theme.SPACING)
        row.addLayout(self.action_slot, 0)

        self.setVisible(False)

    @property
    def tone(self) -> str:
        return self._tone

    def show_message(self, title: str, body: str = "", tone: str = "info") -> None:
        self.title_label.setText(title)
        self.body_label.setText(body)
        self.body_label.setVisible(bool(body))
        self.icon_label.setText(BANNER_ICONS.get(tone, BANNER_ICONS["neutral"]))
        if tone != self._tone:
            self._tone = tone
            for widget in (self, self.icon_label, self.title_label):
                widget.setProperty("tone", tone)
            theme.repolish(self, self.icon_label, self.title_label)
        self.setVisible(True)

    def clear(self) -> None:
        self.title_label.setText("")
        self.body_label.setText("")
        self.setVisible(False)

    def add_action(self, widget: QWidget) -> None:
        self.action_slot.addWidget(widget, 0, Qt.AlignmentFlag.AlignTop)


def clear_layout(layout: QLayout) -> None:
    """Delete every child widget of ``layout``. Used by the lists pages rebuild in place."""
    while layout.count():
        item = layout.takeAt(0)
        widget = item.widget()
        if widget is not None:
            widget.setParent(None)
            widget.deleteLater()
        child = item.layout()
        if child is not None:
            clear_layout(child)


def meta_label(text: str) -> QLabel:
    """The small bold key half of a key/value pair inside a card."""
    label = QLabel(text)
    label.setObjectName("MetaLabel")
    return label


def meta_value(text: str = "") -> QLabel:
    """The value half of a key/value pair inside a card."""
    label = QLabel(text)
    label.setObjectName("MetaValue")
    label.setWordWrap(True)
    label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
    return label


def path_label(text: str = "") -> QLabel:
    """A selectable, wrapping well for a filesystem path.

    Selectable because the single most common thing to do with a path on screen is paste it
    somewhere else, and wrapping rather than elided because a truncated path is not a path.
    """
    label = QLabel(text)
    label.setObjectName("PathValue")
    label.setWordWrap(True)
    label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
    return label


class Card(QFrame):
    """A white 1-px-bordered panel. The only container the design system has.

    Exposes :attr:`body` (a vertical layout, already padded) rather than letting callers
    lay out on the frame itself, so card padding is identical everywhere.
    """

    def __init__(self, title: str = "", parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("Card")
        outer = QVBoxLayout(self)
        outer.setContentsMargins(theme.SPACING * 2, theme.SPACING * 2, theme.SPACING * 2, theme.SPACING * 2)
        outer.setSpacing(theme.SPACING)
        self.title_label = QLabel(title)
        self.title_label.setObjectName("CardTitle")
        self.title_label.setVisible(bool(title))
        outer.addWidget(self.title_label)
        self.body = QVBoxLayout()
        self.body.setContentsMargins(0, 0, 0, 0)
        self.body.setSpacing(theme.SPACING)
        outer.addLayout(self.body)


#: The reading width of an empty state's paragraph. It is one sentence
#: of explanation, centred, and a 900-px line of it would be unreadable.
EMPTY_BODY_WIDTH = 480


class WrappedParagraph(QLabel):
    """A word-wrapped label that reports the height its own text actually needs.

    Qt's word-wrap height is only correct when the layout hands the label the width the
    layout itself chose. A ``QVBoxLayout`` carrying an alignment flag does not: it gives each
    child its ``sizeHint()`` instead, and ``QLabel``'s wrapped ``sizeHint`` is computed for a
    width nobody ends up using. The paragraph is then given three lines of room, wraps to
    five, and, because its text is centred inside a rectangle too short for it, paints
    straight over the headline above and the button below.

    Pinning the width and answering ``sizeHint`` with ``heightForWidth`` of that same width
    makes the two agree, for whatever text is set later.
    """

    def __init__(self, text: str = "", width: int = EMPTY_BODY_WIDTH, parent: QWidget | None = None) -> None:
        super().__init__(text, parent)
        self.setWordWrap(True)
        self._wrap_width = width
        self.setFixedWidth(width)

    def setText(self, text: str) -> None:  # noqa: N802 (Qt's spelling)
        """Set the text, and be visible exactly when there is some.

        A paragraph that has something to say and is hidden anyway is a silent screen, and
        every page that writes into an empty state would otherwise have to remember to unhide
        it. Two of the six did not, and their body copy never appeared at all.
        """
        super().setText(text)
        self.setVisible(bool(text))

    def heightForWidth(self, width: int) -> int:  # noqa: N802, ARG002 (Qt's spelling)
        """The height for the *pinned* width, whatever width the layout is asking about.

        This is the half that actually fixes it. ``QLabel.heightForWidth`` answers for the
        width it is handed and ignores the widget's own maximum, so a layout offering the
        full 842-px column is told the paragraph needs three lines. The label then renders
        it at its fixed 480 px, in five. Answering for the width the text will really
        be laid out at is the only way the two can agree.
        """
        return super().heightForWidth(self._wrap_width)

    def sizeHint(self) -> QSize:  # noqa: N802
        return QSize(self._wrap_width, max(self.heightForWidth(self._wrap_width), 0))

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        return self.sizeHint()


class EmptyState(QFrame):
    """What a page shows when there is nothing on it yet.

    Not an icon and a shrug: a headline naming what is missing, and a body
    naming the action that fixes it. The product spec's rule is *no dead ends: every
    empty state explains what to do next*, and a dashed border reads as "a space waiting to
    be filled" rather than as a failure.
    """

    def __init__(
        self,
        headline: str,
        body: str = "",
        action: QWidget | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("EmptyState")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(theme.SPACING * 5, theme.SPACING * 6, theme.SPACING * 5, theme.SPACING * 6)
        layout.setSpacing(theme.SPACING)
        layout.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.headline_label = QLabel(headline)
        self.headline_label.setObjectName("EmptyStateHeadline")
        self.headline_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.headline_label)

        self.body_label = WrappedParagraph(body)
        self.body_label.setObjectName("EmptyStateBody")
        self.body_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.body_label.setVisible(bool(body))
        layout.addWidget(self.body_label, alignment=Qt.AlignmentFlag.AlignHCenter)

        self.action = action
        if action is not None:
            layout.addSpacing(theme.SPACING)
            layout.addWidget(action, alignment=Qt.AlignmentFlag.AlignHCenter)


def _button(text: str, variant: str) -> QPushButton:
    button = QPushButton(text)
    button.setProperty("variant", variant)
    button.setCursor(Qt.CursorShape.PointingHandCursor)
    return button


def primary_button(text: str) -> QPushButton:
    """The one action a screen most wants the user to take. At most one per screen."""
    return _button(text, "primary")


def secondary_button(text: str) -> QPushButton:
    """Everything else with a border: Back, Browse, Open folder."""
    return _button(text, "secondary")


def ghost_button(text: str) -> QPushButton:
    """A borderless, accent-coloured action for tertiary choices (Clear, Change...)."""
    return _button(text, "ghost")


def set_variant(button: QPushButton, variant: str) -> None:
    """Re-rank a button at runtime, keeping "at most one primary per screen" true.

    The Process screen is the case that needs it: before a run, the hero is *Run analysis*.
    After one, the hero is *Go to Review* and running again is the secondary choice. Two
    filled accent buttons side by side ask the analyst to decide which one the application
    thinks they should press, which is the application's job.
    """
    if button.property("variant") == variant:
        return
    button.setProperty("variant", variant)
    theme.repolish(button)


def hint(text: str) -> QLabel:
    """Small secondary-coloured explanatory text under a control."""
    label = QLabel(text)
    label.setObjectName("Hint")
    label.setWordWrap(True)
    return label


def divider() -> QFrame:
    line = QFrame()
    line.setObjectName("SidebarDivider")
    line.setFrameShape(QFrame.Shape.HLine)
    return line


def scroll_wrap(widget: QWidget) -> QScrollArea:
    """``widget`` in a frameless, vertically-scrolling area styled by the theme.

    A page that scrolls sideways has a layout bug: the column
    is already capped at ``theme.CONTENT_MAX_WIDTH``, and a horizontal scrollbar appearing
    under a form is the classic sign of one. Making the axis unavailable turns that bug into
    a visible squeeze instead of a hidden one.
    """
    area = QScrollArea()
    area.setObjectName("PageScroll")
    area.setWidgetResizable(True)
    area.setFrameShape(QFrame.Shape.NoFrame)
    area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
    area.setWidget(widget)
    return area


# --------------------------------------------------------------------------------------
# The page scaffold
# --------------------------------------------------------------------------------------


class WorkflowPage(QWidget):
    """Base class for the six workflow pages."""

    #: Emitted when the page's header chip changed. The window re-reads
    #: :meth:`page_status` in response. Pages emit it, never the window.
    status_changed = Signal()

    #: ``int``: "please show this step". A page never navigates by itself: it asks, and
    #: :meth:`MainWindow.goto_step` still applies the workflow gate, so a request for a step
    #: that is somehow not unlocked is refused rather than half-honoured. This is how Select
    #: hands off to Assign and Process hands off to Review. Both are *advances the user
    #: expects*, which is why they happen at all. See :meth:`advance_to`.
    navigate_requested = Signal(int)

    #: Which of the six steps this page is. Set by every subclass.
    step: int = 0
    #: Page title, shown in the window's header strip. Defaults to the canonical copy in
    #: :mod:`quant209.app.state` so the sidebar, the header and the page can never disagree
    #: about what a step is called.
    title: str = ""
    #: ``0`` disables the cap (full-bleed pages: chromatograms, heatmaps).
    max_content_width: int = theme.CONTENT_MAX_WIDTH
    #: Opt-in vertical scrolling for pages whose content is taller than the window.
    #:
    #: Off by default: the first four pages each answer one question and
    #: are *designed* to fit, and a scroll area on a page that fits hides content below a
    #: fold nobody knows is there. It is switched on only where the content genuinely grows
    #: with the data: Export, whose second card and success paths appear after a run. A
    #: page that fills its own height instead (Review's tab widget, whose tables and plots
    #: scroll internally) must leave this off, or the elastic child collapses to its minimum
    #: inside the scroll area's unbounded height.
    scrollable: bool = False

    def __init__(self, state: AppState, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.state = state
        if not self.title:
            self.title = STEP_TITLES[self.step]
        self.setObjectName(f"Page{self.__class__.__name__}")

        outer = QVBoxLayout(self)
        outer.setContentsMargins(
            theme.CONTENT_MARGIN, theme.CONTENT_MARGIN, theme.CONTENT_MARGIN, theme.CONTENT_MARGIN
        )
        outer.setSpacing(0)

        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(0)
        self._container = QWidget()
        if self.max_content_width:
            self._container.setMaximumWidth(self.max_content_width)
            # Expanding, not the default Preferred: without this the column is only as wide
            # as its widest child wants to be, so a page with one narrow card would sit in a
            # 540-px gutter instead of filling the reading column up to the cap.
            self._container.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        #: Where subclasses put their widgets.
        self.content = QVBoxLayout(self._container)
        self.content.setContentsMargins(0, 0, 0, 0)
        self.content.setSpacing(theme.SPACING * 2)
        if self.max_content_width:
            # Stretch factor 1 on the column and 0 on the spacer: the column takes the space
            # first and stops at its cap, and only the leftover goes to the spacer. That
            # leftover on the right is what keeps the capped column left-aligned rather than
            # centred. A centred column drifts away from the sidebar as the window widens,
            # and the eye reads the growing gap as a broken layout.
            row.addWidget(self._container, 1)
            row.addStretch(0)
        else:
            row.addWidget(self._container, 1)

        #: The scroll area, when :attr:`scrollable` is set. ``None`` otherwise.
        self.scroll_area: QScrollArea | None = None
        if self.scrollable:
            body = QWidget()
            body.setObjectName("PageScrollBody")
            body.setLayout(row)
            self.scroll_area = scroll_wrap(body)
            outer.addWidget(self.scroll_area, 1)
        else:
            outer.addLayout(row, 1)

        self.build()

    # ---- subclass hooks ---------------------------------------------------------------

    def build(self) -> None:
        """Populate :attr:`content`. Called once, at the end of ``__init__``."""

    def on_activated(self) -> None:
        """The page just became visible. Do expensive imports and first loads here."""

    def refresh(self) -> None:
        """Application state changed while this page is visible."""

    def page_status(self) -> tuple[str, str] | None:
        """``(text, tone)`` for the header chip, or ``None`` for no chip."""
        return None

    # ---- navigation -------------------------------------------------------------------

    def advance_to(self, step: int, delay_ms: int = 0) -> None:
        """Ask the window to move to ``step``, optionally after ``delay_ms``.

        The delay exists for one reason: an auto-advance that happens in the same instant as
        the thing that caused it reads as a glitch rather than as a step. A beat long enough
        to see the success banner land, then the move. It is also what keeps the request out
        of the middle of a modal file dialog's or a worker signal's own call stack, where
        rebuilding the page stack underneath the emitter is asking for trouble.
        """
        if delay_ms <= 0:
            self.navigate_requested.emit(step)
            return
        # ``self`` as the context object, not a bare lambda: Qt then drops the pending call
        # if this page is destroyed before the timer fires. Without it, closing the window
        # during the delay leaves a timer holding a deleted widget, and the emit raises
        # "Signal source has been deleted" out of the event loop, where nothing can catch it.
        QTimer.singleShot(delay_ms, self, lambda: self.navigate_requested.emit(step))
