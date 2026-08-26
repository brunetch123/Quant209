"""The Quant 209 design system: one palette, one stylesheet, applied once at startup.

Why a single hand-written QSS string
-------------------------------------
Default widget chrome showing through is the single largest contributor to an inconsistent
look in any Qt or web app: stock gray dialog backgrounds, native scrollbars, native table
headers, and the platform focus rectangle, mixed with whatever the app itself styles. So
this module does two things together and neither is optional:

* ``QApplication.setStyle("Fusion")``. Fusion is the one Qt style that behaves identically
  on every Windows version and honours every stylesheet property used below. The native
  Windows style silently ignores a good half of them (it draws through the OS theme engine),
  which is exactly how an app ends up half-styled.
* one stylesheet covering **every widget class the app will ever show**, including tables,
  combo boxes, tabs and progress bars. Styling a widget the first time it appears on screen
  is how inconsistency creeps in. The whole system is written down here once.

The palette
------------
:data:`PALETTE` is the single source of colour truth and is exported for pyqtgraph, which
cannot read a Qt stylesheet and would otherwise drift into a second, slightly different set
of colours. :data:`CHART_SERIES` is its categorical companion, ordered for colourblind-safe
pairwise contrast.

Sizing constants (:data:`SPACING`, :data:`CONTENT_MARGIN`, :data:`CONTENT_MAX_WIDTH`,
:data:`SIDEBAR_WIDTH`) live here rather than in the widgets so the 8-px grid is enforced from
one place: every margin, padding and gap in the app is a multiple of :data:`SPACING`.

Font sizes are in **points**, not pixels, so the app respects the user's display scaling. A
lab machine at 125% is the common case, and a px-based UI shrinks to unreadable there.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtGui import QColor, QFont, QFontDatabase, QPalette
from PySide6.QtWidgets import QApplication

#: Where the handful of SVG marks the stylesheet needs live. See :data:`_ASSET_URLS`.
ASSETS = Path(__file__).resolve().parent / "assets"

# --------------------------------------------------------------------------------------
# Palette
# --------------------------------------------------------------------------------------

#: Every colour in the application. Read by the stylesheet below and, later, by the
#: pyqtgraph charts. It is a plain dict of hex strings rather than ``QColor``s because both
#: consumers want the string form.
PALETTE: dict[str, str] = {
    # surfaces ----------------------------------------------------------------------
    "background": "#FAFBFC",  # the page behind the cards
    "surface": "#FFFFFF",  # cards, sidebar, inputs, tables
    "surface_alt": "#F4F6F7",  # hover fills, inert wells
    "surface_sunken": "#EDEFF1",  # progress-bar troughs, scrollbar gutters
    "border": "#E3E6E8",  # the 1 px card/divider border
    "border_strong": "#D4D8DB",  # input and secondary-button borders
    # text --------------------------------------------------------------------------
    "text": "#1A1D1F",
    "text_secondary": "#6A7178",
    "text_disabled": "#A9AFB4",
    # accent ------------------------------------------------------------------------
    "accent": "#0F6E6B",
    "accent_hover": "#0C5956",
    "accent_pressed": "#094643",
    "accent_contrast": "#FFFFFF",
    "accent_soft": "#E7F1F0",  # selection tint, active step row, accent chips
    "accent_disabled": "#E9EEEE",
    # status ------------------------------------------------------------------------
    "success": "#1E7F4F",
    "success_soft": "#E6F3EC",
    "warning": "#B7791F",
    "warning_soft": "#FBF1DF",
    "error": "#B42318",
    "error_soft": "#FBE9E7",
    "info": "#1F5FA8",
    "info_soft": "#E8F0FA",
}

#: Categorical series colours for charts. Ordered so that any two adjacent entries are
#: distinguishable under deuteranopia and protanopia. The accent leads so a single-series
#: chart matches the rest of the UI.
CHART_SERIES: tuple[str, ...] = (
    "#0F6E6B",  # teal (accent)
    "#B7791F",  # amber
    "#1F5FA8",  # blue
    "#B42318",  # red
    "#5B4B8A",  # violet
    "#3F7D20",  # green
    "#8A6A3D",  # brown
    "#4A5459",  # slate
)

# --------------------------------------------------------------------------------------
# Type and metrics
# --------------------------------------------------------------------------------------

# The typeface: the platform's own UI font, at the app's size
# ------------------------------------------------------
# Segoe UI Variable is the intended typeface, with Segoe UI as the fallback. Qt already
# provides the second one for free: ``QApplication.font()`` is Segoe UI on every Windows
# 10/11 machine, resolved from system metrics without touching the font database. Asking for
# the first one by name is expensive in a way that is easy to miss:
#
#   naming *any* font family      +600 ms   (populates the Windows font database)
#   changing only the point size   +40 ms
#
# Measured on the target host with ``tools/measure_startup.py``. ``QFontDatabase.families()``
# alone costs 365 ms. That is a fifth of a 3000 ms cold-start budget spent on the difference
# between Segoe UI and its variable-font refresh. At a 10.5 pt body size, that difference is
# not visible side by side. The cold-start budget exists to keep startup fast, and the
# typeface upgrade is not worth a fifth of it.
#
# So no family is named anywhere, not here and not in the stylesheet, and the app inherits
# the platform UI font. Size, weight and colour are set explicitly. **Do not "fix" this by
# adding a font-family declaration back.** It will cost half a second of startup and change
# almost nothing on screen. If Segoe UI Variable ever becomes worth it, set it *after*
# ``MainWindow.first_painted`` and check for a reflow flash.
FONT_SIZE_PT = 10.5  # body
FONT_SIZE_SMALL_PT = 9  # chips, footnotes, table headers
FONT_SIZE_TITLE_PT = 16  # page titles
FONT_SIZE_WORDMARK_PT = 15

#: The 8-px grid. Every margin/padding/gap in the app is a multiple of this.
SPACING = 8
#: Content inset inside the content area, on all four sides.
CONTENT_MARGIN = 24
#: Content is left-aligned and stops growing here, so a page does not sprawl into unreadable
#: 2000-px line lengths on a wide monitor.
CONTENT_MAX_WIDTH = 980
#: The sidebar width is fixed rather than proportional, so the stepper never reflows.
SIDEBAR_WIDTH = 230

#: Width of the mark above the wordmark. Its height follows from the drawing's own
#: proportions, so replacing the artwork does not need a second number edited here.
MARK_WIDTH = 132

#: The window's minimum size. Below this the six-step sidebar plus a full-width results
#: table stops being usable.
MIN_WINDOW_WIDTH = 1200
MIN_WINDOW_HEIGHT = 760


def _q(name: str) -> QColor:
    return QColor(PALETTE[name])


# --------------------------------------------------------------------------------------
# Stylesheet assets
# --------------------------------------------------------------------------------------


def _asset_url(name: str) -> str:
    """A QSS ``url(...)`` argument for an asset file, with the separators QSS expects.

    Forward slashes on Windows too: Qt's CSS parser treats a backslash as an escape, so a
    native path silently resolves to nothing and the rule paints an empty pixmap.
    """
    return (ASSETS / name).as_posix()


#: Substituted into the stylesheet alongside :data:`PALETTE`.
#:
#: Why the combo-box arrow is an SVG file and not a CSS border triangle
#: ---------------------------------------------------------------------
#: The usual trick, a zero-sized ``::down-arrow`` with two transparent side borders and one
#: coloured top border, which every CSS engine mitres into a triangle, **does not work in
#: Qt**. Qt draws each border edge as its own rectangle, so the rule renders as a solid
#: 8x5 grey block: a blob, not an arrow. (Verified by pixel-checking a rendered ``QComboBox``.
#: ``tests/test_app_pages.py::test_the_combo_box_arrow_renders_as_a_chevron`` keeps it
#: verified.) Nor can the rule simply be deleted and the native arrow left to show through:
#: once a ``QComboBox`` has any box rule at all, ``QStyleSheetStyle`` takes over the whole
#: complex control and never asks the base style to draw the arrow. Removing the rule
#: leaves no arrow whatsoever, and a ``QProxyStyle`` override of ``PE_IndicatorArrowDown`` is
#: never called either.
#:
#: So the arrow is a 10x6 SVG chevron, which Qt's image loader renders crisply at any DPI and
#: which is still text in the repository. It is two files rather than one because a disabled
#: combo box needs the muted stroke. Recolouring is not something QSS can do to an image.
_ASSET_URLS: dict[str, str] = {
    "chevron_down": _asset_url("chevron-down.svg"),
    "chevron_down_disabled": _asset_url("chevron-down-disabled.svg"),
}


# --------------------------------------------------------------------------------------
# Stylesheet
# --------------------------------------------------------------------------------------


def stylesheet() -> str:
    """The whole application stylesheet, with :data:`PALETTE` substituted in.

    Written as one ``str.format``-free template with ``%(name)s`` placeholders: QSS is full
    of braces, so ``str.format`` would need every one of them doubled and the result would
    be unreadable and unmaintainable.

    Two conventions worth knowing before editing:

    * **Focus rings are drawn as a 2 px border with the padding reduced by 1 px** in the same
      rule. QSS has no ``outline`` that Qt honours, so a focus ring has to be a border. If
      the padding is not compensated, focusing a widget nudges its contents by a pixel and
      the whole form twitches. Every ``:focus`` rule below carries its padding partner.
    * **State is carried on Qt dynamic properties** (``[state="active"]``, ``[variant=
      "primary"]``, ``[tone="warning"]``) rather than on object names, so one widget class
      can move between states at runtime. Qt does not re-evaluate a property selector by
      itself. The widget must be unpolished and re-polished. ``repolish`` below is the
      helper that does it.
    """
    return _QSS % {**PALETTE, **_ASSET_URLS}


_QSS = """
/* ---------------------------------------------------------------- base ---- */
QWidget {
    color: %(text)s;
}
QMainWindow, QDialog, QMessageBox {
    background-color: %(background)s;
}
QWidget#ContentArea, QStackedWidget, QStackedWidget > QWidget {
    background-color: %(background)s;
}
QLabel {
    background: transparent;
}
QLabel:disabled {
    color: %(text_disabled)s;
}
QToolTip {
    background-color: %(text)s;
    color: %(surface)s;
    border: none;
    border-radius: 4px;
    padding: 6px 8px;
    font-size: 9pt;
}

/* ------------------------------------------------------------- sidebar ---- */
QWidget#Sidebar {
    background-color: %(surface)s;
    border-right: 1px solid %(border)s;
}
QLabel#Mark {
    padding: 0 0 2px 0;
}
QLabel#Wordmark {
    font-size: 15pt;
    font-weight: 600;
    color: %(accent)s;
}
QLabel#WordmarkTag {
    font-size: 9pt;
    color: %(text_secondary)s;
}
QLabel#SidebarSectionLabel {
    font-size: 9pt;
    font-weight: 600;
    color: %(text_secondary)s;
}
QLabel#SidebarFooter {
    font-size: 9pt;
    color: %(text_secondary)s;
}
QLabel#RunFolderLabel {
    font-size: 9pt;
    color: %(text)s;
    background-color: %(surface_alt)s;
    border: 1px solid %(border)s;
    border-radius: 6px;
    padding: 6px 8px;
}
QFrame#SidebarDivider {
    background-color: %(border)s;
    max-height: 1px;
    min-height: 1px;
    border: none;
}

/* ---------------------------------------------------------- step rows ---- */
/* min-height belongs here, not in StepRow.__init__: QStyleSheetStyle recomputes a styled
   widget's minimum size from its own box model when it polishes, which silently discards a
   programmatic setMinimumHeight. 40px = the 24px number circle plus the row layout's 8px
   above and below; anything less clips the circle. */
QPushButton#StepRow {
    background-color: transparent;
    border: 2px solid transparent;
    border-radius: 6px;
    padding: 0px;
    min-height: 40px;
    text-align: left;
}
QPushButton#StepRow:enabled:hover {
    background-color: %(surface_alt)s;
}
QPushButton#StepRow:enabled:pressed {
    background-color: %(accent_soft)s;
}
QPushButton#StepRow[state="active"] {
    background-color: %(accent_soft)s;
}
QPushButton#StepRow:focus {
    border: 2px solid %(accent)s;
}
QLabel#StepNumber {
    min-width: 24px;
    max-width: 24px;
    min-height: 24px;
    max-height: 24px;
    border-radius: 12px;
    font-size: 9pt;
    font-weight: 600;
    background-color: %(surface_sunken)s;
    color: %(text_secondary)s;
}
QLabel#StepNumber[state="active"] {
    background-color: %(accent)s;
    color: %(accent_contrast)s;
}
QLabel#StepNumber[state="complete"] {
    background-color: %(success)s;
    color: %(accent_contrast)s;
}
QLabel#StepNumber[state="locked"] {
    /* surface_alt, not background: at #FAFBFC on a white sidebar the circle disappears
       entirely and the locked rows read as loose text rather than as steps not yet
       reached. Muted, but still a shape. */
    background-color: %(surface_alt)s;
    color: %(text_disabled)s;
}
QLabel#StepLabel {
    color: %(text)s;
}
QLabel#StepLabel[state="active"] {
    font-weight: 600;
    color: %(accent)s;
}
QLabel#StepLabel[state="locked"] {
    color: %(text_disabled)s;
}

/* ---------------------------------------------------------- page head ---- */
QLabel#PageTitle {
    font-size: 16pt;
    font-weight: 600;
    color: %(text)s;
}
QLabel#PageSubtitle {
    color: %(text_secondary)s;
}
QFrame#HeaderRule {
    background-color: %(border)s;
    max-height: 1px;
    min-height: 1px;
    border: none;
}

/* -------------------------------------------------------------- cards ---- */
QFrame#Card {
    background-color: %(surface)s;
    border: 1px solid %(border)s;
    border-radius: 6px;
}
QFrame#EmptyState {
    background-color: %(surface)s;
    border: 1px dashed %(border_strong)s;
    border-radius: 6px;
}
QLabel#EmptyStateHeadline {
    font-size: 12pt;
    font-weight: 600;
    color: %(text)s;
}
QLabel#EmptyStateBody {
    color: %(text_secondary)s;
}
QLabel#CardTitle {
    font-weight: 600;
    color: %(text)s;
}
QLabel#Hint {
    color: %(text_secondary)s;
    font-size: 9pt;
}

/* --------------------------------------------------------------- chips ---- */
QLabel#StatusChip {
    border-radius: 10px;
    padding: 3px 10px;
    font-size: 9pt;
    font-weight: 600;
    background-color: %(surface_sunken)s;
    color: %(text_secondary)s;
}
QLabel#StatusChip[tone="accent"] {
    background-color: %(accent_soft)s;
    color: %(accent)s;
}
QLabel#StatusChip[tone="success"] {
    background-color: %(success_soft)s;
    color: %(success)s;
}
QLabel#StatusChip[tone="warning"] {
    background-color: %(warning_soft)s;
    color: %(warning)s;
}
QLabel#StatusChip[tone="error"] {
    background-color: %(error_soft)s;
    color: %(error)s;
}
QLabel#StatusChip[tone="info"] {
    background-color: %(info_soft)s;
    color: %(info)s;
}

/* ------------------------------------------------------------ banners ---- */
/* An in-page banner, never a modal. A message box stops the application to say something
   the page could have said in place, and the analyst then has nowhere to read it back. */
QFrame#Banner {
    border-radius: 6px;
    border: 1px solid %(border)s;
    background-color: %(surface_alt)s;
}
QFrame#Banner[tone="error"] {
    border-color: %(error)s;
    background-color: %(error_soft)s;
}
QFrame#Banner[tone="warning"] {
    border-color: %(warning)s;
    background-color: %(warning_soft)s;
}
QFrame#Banner[tone="info"] {
    border-color: %(info)s;
    background-color: %(info_soft)s;
}
QFrame#Banner[tone="success"] {
    border-color: %(success)s;
    background-color: %(success_soft)s;
}
QLabel#BannerIcon {
    min-width: 20px;
    max-width: 20px;
    min-height: 20px;
    max-height: 20px;
    border-radius: 10px;
    font-size: 9pt;
    font-weight: 700;
    background-color: %(text_secondary)s;
    color: %(surface)s;
}
QLabel#BannerIcon[tone="error"] { background-color: %(error)s; }
QLabel#BannerIcon[tone="warning"] { background-color: %(warning)s; }
QLabel#BannerIcon[tone="info"] { background-color: %(info)s; }
QLabel#BannerIcon[tone="success"] { background-color: %(success)s; }
QLabel#BannerTitle {
    font-weight: 600;
    color: %(text)s;
}
QLabel#BannerTitle[tone="error"] { color: %(error)s; }
QLabel#BannerTitle[tone="warning"] { color: %(warning)s; }
QLabel#BannerTitle[tone="info"] { color: %(info)s; }
QLabel#BannerTitle[tone="success"] { color: %(success)s; }
QLabel#BannerBody {
    color: %(text_secondary)s;
}

/* ---------------------------------------------------- card inner parts ---- */
QLabel#MetaLabel {
    font-size: 9pt;
    font-weight: 600;
    color: %(text_secondary)s;
}
QLabel#MetaValue {
    color: %(text)s;
}
QLabel#PathValue {
    color: %(text)s;
    background-color: %(surface_alt)s;
    border: 1px solid %(border)s;
    border-radius: 6px;
    padding: 7px 10px;
}
QLabel#PlanRow {
    color: %(text)s;
    padding: 2px 0px;
}
QLabel#PlanRow[tone="warning"] { color: %(warning)s; }
QLabel#PlanRow[tone="error"] { color: %(error)s; }
QLabel#PlanRow[tone="success"] { color: %(success)s; }
QLabel#RunningLabel {
    color: %(text)s;
    font-weight: 600;
}

/* --------------------------------------------------------------- logs ---- */
QListWidget#RunLog {
    background-color: %(surface_alt)s;
    border: 1px solid %(border)s;
    border-radius: 6px;
    padding: 4px;
}
QListWidget#RunLog::item {
    padding: 1px 6px;
    border: none;
    color: %(text_secondary)s;
}
QListWidget#RunLog::item:selected {
    background: transparent;
    color: %(text_secondary)s;
}
QPlainTextEdit#Traceback {
    background-color: %(surface_alt)s;
    border: 1px solid %(border)s;
    border-radius: 6px;
    color: %(text_secondary)s;
}

/* ------------------------------------------------------------ buttons ---- */
QPushButton {
    background-color: %(surface)s;
    color: %(text)s;
    border: 1px solid %(border_strong)s;
    border-radius: 6px;
    padding: 7px 16px;
    min-height: 18px;
}
QPushButton:hover {
    background-color: %(surface_alt)s;
    border-color: %(text_disabled)s;
}
QPushButton:pressed {
    background-color: %(surface_sunken)s;
}
QPushButton:focus {
    border: 2px solid %(accent)s;
    padding: 6px 15px;
}
QPushButton:disabled {
    background-color: %(background)s;
    color: %(text_disabled)s;
    border-color: %(border)s;
}
QPushButton[variant="primary"] {
    background-color: %(accent)s;
    color: %(accent_contrast)s;
    border: 1px solid %(accent)s;
    font-weight: 600;
}
QPushButton[variant="primary"]:hover {
    background-color: %(accent_hover)s;
    border-color: %(accent_hover)s;
}
QPushButton[variant="primary"]:pressed {
    background-color: %(accent_pressed)s;
    border-color: %(accent_pressed)s;
}
QPushButton[variant="primary"]:focus {
    border: 2px solid %(accent_pressed)s;
    padding: 6px 15px;
}
QPushButton[variant="primary"]:disabled {
    background-color: %(accent_disabled)s;
    border-color: %(accent_disabled)s;
    color: %(text_disabled)s;
}
QPushButton[variant="ghost"] {
    background-color: transparent;
    border: 1px solid transparent;
    color: %(accent)s;
    font-weight: 600;
}
QPushButton[variant="ghost"]:hover {
    background-color: %(accent_soft)s;
}
QPushButton[variant="ghost"]:disabled {
    color: %(text_disabled)s;
    background-color: transparent;
}

/* ------------------------------------------------------------- inputs ---- */
QLineEdit, QPlainTextEdit, QTextEdit, QSpinBox, QDoubleSpinBox, QComboBox {
    background-color: %(surface)s;
    color: %(text)s;
    border: 1px solid %(border_strong)s;
    border-radius: 6px;
    padding: 6px 10px;
    min-height: 18px;
    selection-background-color: %(accent)s;
    selection-color: %(accent_contrast)s;
}
QLineEdit:focus, QPlainTextEdit:focus, QTextEdit:focus,
QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {
    border: 2px solid %(accent)s;
    padding: 5px 9px;
}
QLineEdit:disabled, QPlainTextEdit:disabled, QTextEdit:disabled,
QSpinBox:disabled, QDoubleSpinBox:disabled, QComboBox:disabled {
    background-color: %(background)s;
    color: %(text_disabled)s;
    border-color: %(border)s;
}
QLineEdit[state="error"] {
    border: 1px solid %(error)s;
}
QComboBox::drop-down {
    border: none;
    width: 22px;
}
QComboBox::down-arrow {
    /* An SVG, not a CSS border triangle -- see _ASSET_URLS for why the trick fails in Qt. */
    image: url("%(chevron_down)s");
    width: 10px;
    height: 6px;
}
QComboBox::down-arrow:disabled {
    image: url("%(chevron_down_disabled)s");
}
QComboBox QAbstractItemView {
    background-color: %(surface)s;
    border: 1px solid %(border)s;
    border-radius: 6px;
    padding: 4px;
    outline: none;
    selection-background-color: %(accent_soft)s;
    selection-color: %(text)s;
}
QCheckBox, QRadioButton {
    spacing: 8px;
    background: transparent;
}
QCheckBox::indicator, QRadioButton::indicator {
    width: 16px;
    height: 16px;
    background-color: %(surface)s;
    border: 1px solid %(border_strong)s;
}
QCheckBox::indicator {
    border-radius: 3px;
}
QRadioButton::indicator {
    border-radius: 9px;
}
QCheckBox::indicator:hover, QRadioButton::indicator:hover {
    border-color: %(accent)s;
}
QCheckBox::indicator:checked, QRadioButton::indicator:checked {
    background-color: %(accent)s;
    border-color: %(accent)s;
}
QCheckBox::indicator:disabled, QRadioButton::indicator:disabled {
    background-color: %(background)s;
    border-color: %(border)s;
}

/* ------------------------------------------------------------- tables ---- */
QTableView, QTreeView, QListView {
    background-color: %(surface)s;
    alternate-background-color: %(surface)s;
    border: 1px solid %(border)s;
    border-radius: 6px;
    gridline-color: transparent;
    outline: none;
    selection-background-color: %(accent_soft)s;
    selection-color: %(text)s;
}
QTableView::item, QTreeView::item, QListView::item {
    padding: 6px 10px;
    border: none;
    border-bottom: 1px solid %(background)s;
}
QTableView::item:hover, QTreeView::item:hover, QListView::item:hover {
    background-color: %(surface_alt)s;
}
QTableView::item:selected, QTreeView::item:selected, QListView::item:selected {
    background-color: %(accent_soft)s;
    color: %(text)s;
}
QHeaderView {
    background-color: %(surface)s;
}
QHeaderView::section {
    background-color: %(surface)s;
    color: %(text_secondary)s;
    font-size: 9pt;
    font-weight: 600;
    padding: 8px 10px;
    border: none;
    border-bottom: 1px solid %(border)s;
}
QHeaderView::section:hover {
    color: %(text)s;
}
QTableCornerButton::section {
    background-color: %(surface)s;
    border: none;
    border-bottom: 1px solid %(border)s;
}

/* --------------------------------------------------------------- tabs ---- */
QTabWidget::pane {
    border: 1px solid %(border)s;
    border-radius: 6px;
    background-color: %(surface)s;
    top: -1px;
}
QTabBar {
    qproperty-drawBase: 0;
}
QTabBar::tab {
    background: transparent;
    color: %(text_secondary)s;
    border: none;
    border-bottom: 2px solid transparent;
    padding: 8px 16px;
    margin-right: 8px;
}
QTabBar::tab:hover {
    color: %(text)s;
}
QTabBar::tab:selected {
    color: %(accent)s;
    font-weight: 600;
    border-bottom: 2px solid %(accent)s;
}

/* ------------------------------------------------------ review page ---- */
/* A tab's own ground is the page background, not the pane's white: the Review tabs are
   built out of white cards and white tables, and a white card on a white pane has only its
   1-px border to say it is a card at all. */
QWidget#ReviewTab {
    background-color: %(background)s;
}

/* The dashboard's three genuinely new widget classes. Everything else on Review reuses
   Card / StatusChip / Banner / QTableView above -- a fifth kind of panel on the one screen
   an analyst spends the most time on is how a design system stops being one. */

/* A stat tile: one number, one caption. Same frame as Card so the two sit on the same
   grid; the number carries the tone, because a red "3" reads before any word does. */
QFrame#StatCard {
    background-color: %(surface)s;
    border: 1px solid %(border)s;
    border-radius: 6px;
}
QLabel#StatValue {
    font-size: 17pt;
    font-weight: 600;
    color: %(text)s;
}
QLabel#StatValue[tone="success"] { color: %(success)s; }
QLabel#StatValue[tone="warning"] { color: %(warning)s; }
QLabel#StatValue[tone="error"] { color: %(error)s; }
QLabel#StatValue[tone="accent"] { color: %(accent)s; }
QLabel#StatCaption {
    font-size: 9pt;
    font-weight: 600;
    color: %(text_secondary)s;
}

/* The QA/QC filter chips. A checkable button, so it needs a :checked rule -- without one
   Qt paints a pressed-in checked button exactly like an unpressed one under a stylesheet,
   and a filter row where nothing shows which filter is on is worse than no filter row. */
QPushButton[variant="filter"] {
    background-color: %(surface)s;
    border: 1px solid %(border_strong)s;
    border-radius: 12px;
    padding: 4px 12px;
    min-height: 16px;
    font-size: 9pt;
    font-weight: 600;
    color: %(text_secondary)s;
}
QPushButton[variant="filter"]:hover {
    background-color: %(surface_alt)s;
    color: %(text)s;
}
QPushButton[variant="filter"]:checked {
    background-color: %(accent_soft)s;
    border-color: %(accent)s;
    color: %(accent)s;
}
QPushButton[variant="filter"]:focus {
    border: 2px solid %(accent)s;
    padding: 3px 11px;
}

/* The recovery heatmap's column headers. The ordinary 10-px side padding is right for a
   header that labels a column of prose; here it labels a 55-px cell holding a percentage, and the
   twenty pixels it takes are the difference between a readable file id and an elided one -- a
   sample id that is not merely hard to read but *wrong*. The cells themselves are painted by
   RecoveryCellDelegate and have no padding to match. */
QHeaderView#HeatmapHeader::section {
    padding: 8px 3px;
}

/* The chromatogram browser's file rail. A plain QListView otherwise, but it sits directly
   against the plot, so it keeps the card border rather than floating. */
QListWidget#FileRail {
    border: 1px solid %(border)s;
    border-radius: 6px;
    padding: 2px;
}
QListWidget#FileRail::item {
    padding: 6px 10px;
    border-radius: 4px;
    border: none;
}

/* ---------------------------------------------------------- progress ---- */
QProgressBar {
    background-color: %(surface_sunken)s;
    border: none;
    border-radius: 4px;
    min-height: 8px;
    max-height: 8px;
    text-align: center;
    color: %(text_secondary)s;
}
QProgressBar::chunk {
    background-color: %(accent)s;
    border-radius: 4px;
}

/* ------------------------------------------------------- scroll areas ---- */
QScrollArea {
    background: transparent;
    border: none;
}
QScrollArea > QWidget > QWidget {
    background: transparent;
}
QScrollBar:vertical {
    background: transparent;
    width: 10px;
    margin: 0px;
}
QScrollBar::handle:vertical {
    background-color: %(border_strong)s;
    border-radius: 5px;
    min-height: 32px;
}
QScrollBar::handle:vertical:hover {
    background-color: %(text_disabled)s;
}
QScrollBar:horizontal {
    background: transparent;
    height: 10px;
    margin: 0px;
}
QScrollBar::handle:horizontal {
    background-color: %(border_strong)s;
    border-radius: 5px;
    min-width: 32px;
}
QScrollBar::handle:horizontal:hover {
    background-color: %(text_disabled)s;
}
QScrollBar::add-line, QScrollBar::sub-line {
    width: 0px;
    height: 0px;
    border: none;
    background: none;
}
QScrollBar::add-page, QScrollBar::sub-page {
    background: none;
}

/* ------------------------------------------------------------- menus ---- */
QMenuBar {
    background-color: %(surface)s;
    border-bottom: 1px solid %(border)s;
}
QMenuBar::item {
    padding: 6px 12px;
    background: transparent;
}
QMenuBar::item:selected {
    background-color: %(accent_soft)s;
    color: %(accent)s;
}
QMenu {
    background-color: %(surface)s;
    border: 1px solid %(border)s;
    border-radius: 6px;
    padding: 4px;
}
QMenu::item {
    padding: 6px 24px 6px 12px;
    border-radius: 4px;
}
QMenu::item:selected {
    background-color: %(accent_soft)s;
    color: %(text)s;
}
QMenu::separator {
    height: 1px;
    background-color: %(border)s;
    margin: 4px 8px;
}

/* ---------------------------------------------------------- splitter ---- */
QSplitter::handle {
    background-color: %(border)s;
}
QSplitter::handle:horizontal {
    width: 1px;
}
QSplitter::handle:vertical {
    height: 1px;
}
"""


# --------------------------------------------------------------------------------------
# Application wiring
# --------------------------------------------------------------------------------------


def base_font(app: QApplication) -> QFont:
    """The platform's UI font at :data:`FONT_SIZE_PT`. See the typeface note above.

    Derived from ``app.font()`` rather than constructed from a family name, which is what
    keeps the font database out of the startup path. The default is 9 pt, which is a
    point-and-a-half tighter than this design wants. The size is the part worth changing.
    """
    font = QFont(app.font())
    font.setPointSizeF(FONT_SIZE_PT)
    return font


def mono_font(point_size: float = FONT_SIZE_SMALL_PT) -> QFont:
    """The platform's fixed-pitch UI font, for the run log and tracebacks.

    From ``QFontDatabase.systemFont(FixedFont)`` rather than a family name: the platform theme
    answers that from system parameters, where naming a family ("Consolas", "Courier New")
    populates the whole Windows font database and costs ~600 ms, the same trap documented
    above :data:`FONT_SIZE_PT`. It would be paid as a freeze the first time the Process
    page opened rather than at startup, which is worse, not better.
    """
    font = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
    font.setPointSizeF(point_size)
    return font


def qt_palette() -> QPalette:
    """A ``QPalette`` matching :data:`PALETTE`.

    The stylesheet covers every widget the app draws, but not everything Qt draws *for*
    it: a native file dialog's item view, a tooltip's drop shadow, the text cursor
    colour. Those read the palette, so the palette has to agree with the stylesheet or
    the stock gray reappears in exactly the places the user notices last and trusts least.
    """
    p = QPalette()
    p.setColor(QPalette.ColorRole.Window, _q("background"))
    p.setColor(QPalette.ColorRole.WindowText, _q("text"))
    p.setColor(QPalette.ColorRole.Base, _q("surface"))
    p.setColor(QPalette.ColorRole.AlternateBase, _q("surface_alt"))
    p.setColor(QPalette.ColorRole.Text, _q("text"))
    p.setColor(QPalette.ColorRole.Button, _q("surface"))
    p.setColor(QPalette.ColorRole.ButtonText, _q("text"))
    p.setColor(QPalette.ColorRole.Highlight, _q("accent"))
    p.setColor(QPalette.ColorRole.HighlightedText, _q("accent_contrast"))
    p.setColor(QPalette.ColorRole.ToolTipBase, _q("text"))
    p.setColor(QPalette.ColorRole.ToolTipText, _q("surface"))
    p.setColor(QPalette.ColorRole.PlaceholderText, _q("text_secondary"))
    p.setColor(QPalette.ColorRole.Link, _q("accent"))
    for group in (QPalette.ColorGroup.Disabled,):
        p.setColor(group, QPalette.ColorRole.WindowText, _q("text_disabled"))
        p.setColor(group, QPalette.ColorRole.Text, _q("text_disabled"))
        p.setColor(group, QPalette.ColorRole.ButtonText, _q("text_disabled"))
    return p


def apply_theme(app: QApplication) -> None:
    """Apply the design system to ``app``: Fusion base style, palette, font, stylesheet.

    Idempotent: calling it twice (the tests do) leaves the application in the same state.
    """
    app.setStyle("Fusion")
    app.setPalette(qt_palette())
    app.setFont(base_font(app))  # after setStyle: Fusion may carry its own default font
    app.setStyleSheet(stylesheet())


def repolish(*widgets) -> None:
    """Re-evaluate stylesheet rules for widgets whose dynamic properties just changed.

    Qt matches ``[state="active"]``-style selectors when a widget is *polished*, and does
    not watch the property afterwards. Anything that calls ``setProperty`` for styling has
    to follow it with this, or the widget keeps the look it had when it was created: the
    classic "my QSS state selector does nothing" bug.
    """
    for widget in widgets:
        style = widget.style()
        style.unpolish(widget)
        style.polish(widget)
        widget.update()
