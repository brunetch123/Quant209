"""Step 2, Assign: confirm the file role detected from each filename.

What this screen is for
------------------------
Every downstream number depends on which files are the calibration injections, which is the
reference, and which are blanks. Those are guessed from the filename convention
(``core.sequence.detect_type``), and the guess is right almost always. But "almost always"
is not a basis for a quantification, so the guess is shown, labelled as a guess, and made
editable in one click.

Hence the two type columns. **Detected type** is what the filename said and never changes.
**Type** is what the run will actually be processed with. While they agree, the detected
cell carries a small ``auto`` badge. The moment the analyst overrides one, the badge goes and
the two columns visibly disagree, which is the whole point. An application that silently
overwrote the detection would leave nothing on screen to show that a human had intervened,
and this run's audit trail would be one row shorter than the truth.

Why a model/view table rather than a widget table
--------------------------------------------------
``QTableWidget`` stores the data in the cells, which makes the widget the source of truth for
something :class:`~quant209.app.state.AppState` already owns. Two copies of the file list
drift the first time anything else changes one. Something else does: the Convert step
writes mzML, which changes this table's Status column without touching a single type. A
:class:`QAbstractTableModel` reading ``state.entries`` has one copy, and refreshing is a
signal rather than a rebuild.

The role dropdown is a **persistent** editor
---------------------------------------------
Every row's Type cell holds a real, open combo box rather than a piece of text that turns
into one when double-clicked. A cell that looks like a label but is secretly editable is
unclear how to use: nothing on screen says the type can be changed, so nobody changes it.
Thirty combo boxes cost nothing at this scale and the affordance is unmistakable.

The banners
------------
Three run-level facts are worth interrupting for, at three different weights, and all three
are about what will happen *later*. No calibration injection means nothing can be quantified
at all (error). No reference injection means no recoveries (warning). No blanks means no LOQ
and therefore no censoring (info). They are stated here, on the screen that can fix them,
rather than as a failure three steps later.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from PySide6.QtCore import QAbstractTableModel, QModelIndex, QRectF, QSize, Qt
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QComboBox,
    QHBoxLayout,
    QHeaderView,
    QStyle,
    QStyledItemDelegate,
    QStyleOptionViewItem,
    QTableView,
)

from quant209.app import theme
from quant209.app.state import STEP_ASSIGN, STEP_CONVERT, STEP_PROCESS
from quant209.app.workflow.base import (
    ROLE_LABELS,
    ROLE_ORDER,
    Banner,
    ChipRow,
    EmptyState,
    WorkflowPage,
    meta_label,
    primary_button,
    role_chips,
)

# --------------------------------------------------------------------------------------
# Table shape
# --------------------------------------------------------------------------------------

COL_FILE = 0
COL_DETECTED = 1
COL_TYPE = 2
COL_SIZE = 3
COL_STATUS = 4
COLUMN_COUNT = 5

HEADERS: tuple[str, ...] = ("File", "Detected type", "Type", "Size", "Status")

#: ``index.data(AUTO_ROLE)``: True while the row still carries the detected type. Read by
#: :class:`DetectedDelegate` to decide whether to draw the ``auto`` badge.
AUTO_ROLE = int(Qt.ItemDataRole.UserRole) + 1

#: The badge itself. Lower case and small: it is a footnote on the value, not a second value.
AUTO_BADGE = "auto"

#: Status column vocabulary. "Ready" and "Converted" are different words because one file
#: arrived as mzML and one was converted here, and an analyst chasing a bad number wants to
#: know which.
STATUS_READY = "Ready (mzML)"
STATUS_CONVERTED = "Converted"
STATUS_PENDING = "Needs conversion"
STATUS_NOT_READ = "Not read"


@dataclass
class _Row:
    """One table row, computed once per refresh rather than per repaint.

    The Size and Status columns both ask the filesystem, and ``data()`` is called for every
    visible cell on every repaint. A scroll would otherwise ``stat`` a run's worth of
    directories tens of times a second.
    """

    path: Path
    detected: str
    file_type: str
    auto: bool
    size_text: str
    status_text: str


def _format_size(n_bytes: int) -> str:
    if n_bytes <= 0:
        return "--"
    for unit, scale in (("GB", 1 << 30), ("MB", 1 << 20), ("KB", 1 << 10)):
        if n_bytes >= scale:
            return f"{n_bytes / scale:.1f} {unit}"
    return f"{n_bytes} B"


class EntryTableModel(QAbstractTableModel):
    """``state.entries`` as a five-column table. The state is the source of truth."""

    def __init__(self, state, parent=None) -> None:
        super().__init__(parent)
        self.state = state
        self._rows: list[_Row] = []
        #: path -> size in bytes. An Agilent ``.D`` is a directory of a dozen files and the
        #: answer cannot change while the run is open, so it is measured once.
        self._sizes: dict[Path, int] = {}

    # ---- population ---------------------------------------------------------------------

    def sync(self) -> bool:
        """Re-read ``state.entries``. Returns True if the row *set* changed (a full reset).

        A type override changes cells, not rows, and a reset would close every persistent
        combo box on the table, including the one the analyst is in the middle of using.
        So the two cases are distinguished by the file paths: same files, same rows.
        """
        rows = [self._row_for(entry) for entry in self.state.entries]
        structural = [r.path for r in rows] != [r.path for r in self._rows]
        if structural:
            self.beginResetModel()
            self._rows = rows
            self.endResetModel()
            return True
        self._rows = rows
        if rows:
            self.dataChanged.emit(
                self.index(0, 0), self.index(len(rows) - 1, COLUMN_COUNT - 1)
            )
        return False

    def _row_for(self, entry) -> _Row:
        from quant209.core.sequence import detect_type

        path = entry.path
        if path not in self._sizes:
            self._sizes[path] = self._measure(path)
        mzml = self.state.mzml_for(entry)
        if path.suffix.lower() == ".mzml":
            status = STATUS_READY if path.is_file() else STATUS_NOT_READ
        elif mzml is not None and mzml.is_file():
            status = STATUS_CONVERTED
        else:
            status = STATUS_PENDING
        return _Row(
            path=path,
            # The *detected* type is re-derived from the filename rather than remembered.
            # ``FileEntry`` keeps only the type in force, so once a row is overridden the
            # original is gone. This column exists to still show it.
            detected=detect_type(path),
            file_type=entry.file_type,
            auto=bool(entry.auto),
            size_text=_format_size(self._sizes[path]),
            status_text=status,
        )

    @staticmethod
    def _measure(path: Path) -> int:
        try:
            if path.is_dir():
                return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())
            return path.stat().st_size
        except OSError:  # a file removed mid-scan: report nothing rather than crash a repaint
            return 0

    def row_at(self, row: int) -> _Row:
        return self._rows[row]

    # ---- QAbstractTableModel ------------------------------------------------------------

    def rowCount(self, parent=QModelIndex()) -> int:  # noqa: N802, B008 (Qt's spelling)
        return 0 if parent.isValid() else len(self._rows)

    def columnCount(self, parent=QModelIndex()) -> int:  # noqa: N802, B008
        return 0 if parent.isValid() else COLUMN_COUNT

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):  # noqa: N802
        if orientation != Qt.Orientation.Horizontal or role != Qt.ItemDataRole.DisplayRole:
            return None
        return HEADERS[section]

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        row = self._rows[index.row()]
        column = index.column()
        if role == Qt.ItemDataRole.DisplayRole:
            if column == COL_FILE:
                return row.path.name
            if column == COL_DETECTED:
                return ROLE_LABELS.get(row.detected, row.detected)
            if column == COL_TYPE:
                return ROLE_LABELS.get(row.file_type, row.file_type)
            if column == COL_SIZE:
                return row.size_text
            if column == COL_STATUS:
                return row.status_text
        if role == Qt.ItemDataRole.EditRole and column == COL_TYPE:
            return row.file_type
        if role == AUTO_ROLE:
            return row.auto
        if role == Qt.ItemDataRole.ToolTipRole:
            return str(row.path)
        if role == Qt.ItemDataRole.ForegroundRole and column == COL_STATUS:
            if row.status_text == STATUS_PENDING:
                return QColor(theme.PALETTE["warning"])
            if row.status_text == STATUS_NOT_READ:
                return QColor(theme.PALETTE["error"])
            return QColor(theme.PALETTE["success"])
        return None

    def flags(self, index):
        base = super().flags(index)
        if index.isValid() and index.column() == COL_TYPE:
            return base | Qt.ItemFlag.ItemIsEditable
        return base

    def setData(self, index, value, role=Qt.ItemDataRole.EditRole) -> bool:  # noqa: N802
        if role != Qt.ItemDataRole.EditRole or index.column() != COL_TYPE:
            return False
        if value not in ROLE_ORDER:
            return False
        # Straight through to the state, which sets ``auto=False`` and emits. The page's
        # handler calls ``sync`` and every derived cell (and the badge) updates from there.
        self.state.set_entry_type(index.row(), value)
        return True


# --------------------------------------------------------------------------------------
# Delegates
# --------------------------------------------------------------------------------------


class DetectedDelegate(QStyledItemDelegate):
    """Draws the detected type, followed by an ``auto`` pill while nobody has overridden it.

    A second column of plain text would say the same words. The pill is what makes "this is
    the machine's guess" legible at a glance down thirty rows, and its *absence* is what
    marks the rows a human touched.
    """

    #: Matches ``QTableView::item``'s horizontal padding in the stylesheet, so the text in
    #: this column sits on the same left edge as every other column's.
    PADDING = 10
    GAP = 8

    def paint(self, painter: QPainter, option, index) -> None:
        opt = QStyleOptionViewItem(option)
        self.initStyleOption(opt, index)
        text = opt.text
        opt.text = ""  # the background/selection is the style's job, the text below is not
        widget = opt.widget
        style = widget.style() if widget is not None else QApplication.style()
        style.drawControl(QStyle.ControlElement.CE_ItemViewItem, opt, painter, widget)

        rect = opt.rect.adjusted(self.PADDING, 0, -self.PADDING, 0)
        painter.save()
        painter.setFont(opt.font)
        painter.setPen(QColor(theme.PALETTE["text"]))
        painter.drawText(
            rect,
            int(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter),
            text,
        )
        if index.data(AUTO_ROLE):
            metrics = QFontMetrics(opt.font)
            self._draw_badge(painter, rect, metrics.horizontalAdvance(text), opt.font)
        painter.restore()

    def _draw_badge(self, painter: QPainter, rect, text_width: int, font: QFont) -> None:
        badge_font = QFont(font)
        badge_font.setPointSizeF(theme.FONT_SIZE_SMALL_PT)
        badge_font.setBold(True)
        metrics = QFontMetrics(badge_font)
        width = metrics.horizontalAdvance(AUTO_BADGE) + 14
        height = metrics.height() + 2
        pill = QRectF(
            rect.left() + text_width + self.GAP,
            rect.center().y() - height / 2 + 1,
            width,
            height,
        )
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor(theme.PALETTE["surface_sunken"]))
        painter.drawRoundedRect(pill, height / 2, height / 2)
        painter.setPen(QColor(theme.PALETTE["text_secondary"]))
        painter.setFont(badge_font)
        painter.drawText(pill, int(Qt.AlignmentFlag.AlignCenter), AUTO_BADGE)

    def sizeHint(self, option, index) -> QSize:  # noqa: N802
        size = super().sizeHint(option, index)
        badge_font = QFont(option.font)
        badge_font.setPointSizeF(theme.FONT_SIZE_SMALL_PT)
        badge_font.setBold(True)
        extra = QFontMetrics(badge_font).horizontalAdvance(AUTO_BADGE) + 14 + self.GAP
        return QSize(size.width() + extra, size.height())


class RoleDelegate(QStyledItemDelegate):
    """A combo box of the five file roles, held permanently open on every row.

    ``activated`` rather than ``currentIndexChanged``: the first fires only when a *person*
    picks something, the second also fires when the view pushes new data into the editor
    after a refresh, which would turn every refresh into an edit, and an edit into a
    refresh, forever.
    """

    def createEditor(self, parent, option, index) -> QComboBox:  # noqa: N802
        combo = QComboBox(parent)
        for role in ROLE_ORDER:
            combo.addItem(ROLE_LABELS[role], role)
        combo.setCursor(Qt.CursorShape.PointingHandCursor)
        combo.activated.connect(lambda _i, editor=combo: self._commit(editor))
        return combo

    def _commit(self, editor: QComboBox) -> None:
        self.commitData.emit(editor)

    def setEditorData(self, editor: QComboBox, index) -> None:  # noqa: N802
        position = editor.findData(index.data(Qt.ItemDataRole.EditRole))
        if position >= 0 and position != editor.currentIndex():
            editor.setCurrentIndex(position)

    def setModelData(self, editor: QComboBox, model, index) -> None:  # noqa: N802
        model.setData(index, editor.currentData(), Qt.ItemDataRole.EditRole)


# --------------------------------------------------------------------------------------
# The page
# --------------------------------------------------------------------------------------


class AssignPage(WorkflowPage):
    step = STEP_ASSIGN

    #: Enough height for a run's worth of rows without the page itself scrolling.
    ROW_HEIGHT = 40
    TYPE_COLUMN_WIDTH = 170

    def build(self) -> None:
        self.banner_no_cal = Banner()
        self.banner_no_ref = Banner()
        self.banner_no_blank = Banner()
        for banner in (self.banner_no_cal, self.banner_no_ref, self.banner_no_blank):
            self.content.addWidget(banner)

        summary = QHBoxLayout()
        summary.setContentsMargins(0, 0, 0, 0)
        summary.setSpacing(theme.SPACING + 4)
        summary.addWidget(meta_label("FILE TYPES"))
        self.chips = ChipRow()
        summary.addWidget(self.chips, 1)
        self.content.addLayout(summary)

        self.model = EntryTableModel(self.state, self)
        self.table = QTableView()
        self.table.setModel(self.model)
        self.table.setObjectName("EntryTable")
        self.table.setAlternatingRowColors(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setShowGrid(False)
        self.table.verticalHeader().setVisible(False)
        self.table.verticalHeader().setDefaultSectionSize(self.ROW_HEIGHT)
        self.table.horizontalHeader().setHighlightSections(False)
        self.table.setWordWrap(False)
        self.table.setCornerButtonEnabled(False)
        self.detected_delegate = DetectedDelegate(self.table)
        self.role_delegate = RoleDelegate(self.table)
        self.table.setItemDelegateForColumn(COL_DETECTED, self.detected_delegate)
        self.table.setItemDelegateForColumn(COL_TYPE, self.role_delegate)

        header = self.table.horizontalHeader()
        # Left, not Fusion's centred default: every cell in this table is left-aligned, and a
        # centred caption floating over a wide left-aligned column is the single most visible
        # sign of a table nobody laid out on purpose.
        header.setDefaultAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
        )
        header.setSectionResizeMode(COL_FILE, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(COL_DETECTED, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(COL_TYPE, QHeaderView.ResizeMode.Fixed)
        header.setSectionResizeMode(COL_SIZE, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(COL_STATUS, QHeaderView.ResizeMode.ResizeToContents)
        self.table.setColumnWidth(COL_TYPE, self.TYPE_COLUMN_WIDTH)
        self.content.addWidget(self.table, 1)

        self.empty = EmptyState(
            "No files to assign",
            "Select a run folder first.",
        )
        self.content.addWidget(self.empty)

        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.addStretch(1)
        self.continue_button = primary_button("Continue")
        self.continue_button.clicked.connect(self._continue)
        row.addWidget(self.continue_button)
        self.content.addLayout(row)

        self.state.entries_changed.connect(self._on_entries_changed)
        self.refresh()

    # ---- behaviour ----------------------------------------------------------------------

    def _on_entries_changed(self, _entries) -> None:
        self.refresh()

    def _continue(self) -> None:
        self.advance_to(self._next_step())

    def _next_step(self) -> int:
        """Process when there is nothing left to convert, Convert when there is."""
        return STEP_PROCESS if self.state.mzml_ready() else STEP_CONVERT

    def open_role_editors(self) -> None:
        """Hold every row's role combo box open. See the module docstring."""
        for row in range(self.model.rowCount()):
            index = self.model.index(row, COL_TYPE)
            if not self.table.isPersistentEditorOpen(index):
                self.table.openPersistentEditor(index)

    # ---- rendering ----------------------------------------------------------------------

    def refresh(self) -> None:
        entries = self.state.entries
        self.model.sync()
        self.open_role_editors()

        has_files = bool(entries)
        self.table.setVisible(has_files)
        self.empty.setVisible(not has_files)
        self.continue_button.setVisible(has_files)
        self.chips.set_chips(role_chips(entries))
        self._refresh_banners(entries)
        self.continue_button.setText(
            "Continue to Process" if self.state.mzml_ready() else "Continue to Convert"
        )
        self.status_changed.emit()

    def _refresh_banners(self, entries) -> None:
        types = {entry.file_type for entry in entries}
        if not entries:
            for banner in (self.banner_no_cal, self.banner_no_ref, self.banner_no_blank):
                banner.clear()
            return

        if "cal" not in types:
            self.banner_no_cal.show_message(
                "A calibration injection is required",
                "The retention-time anchors and the response factors both come from one, so "
                "this run cannot be quantified until a file is assigned the calibration "
                "type.",
                tone="error",
            )
        else:
            self.banner_no_cal.clear()

        if "ref" not in types:
            self.banner_no_ref.show_message(
                "No reference injection: recoveries will be unavailable",
                "Surrogate recoveries are measured against a reference injection. The run "
                "will still be quantified, with uncorrected masses.",
                tone="warning",
            )
        else:
            self.banner_no_ref.clear()

        if "blank" not in types:
            self.banner_no_blank.show_message(
                "No blanks: LOQ censoring disabled",
                "The limit of quantification is the spread of this run's blanks, so without "
                "them no result can be censored as below the laboratory background.",
                tone="info",
            )
        else:
            self.banner_no_blank.clear()

    def page_status(self) -> tuple[str, str] | None:
        entries = self.state.entries
        if not entries:
            return None
        if not any(e.file_type == "cal" for e in entries):
            return ("No calibration injection", "error")
        overridden = sum(1 for e in entries if not e.auto)
        if overridden:
            return (f"{overridden} overridden", "accent")
        return ("Auto-detected", "neutral")
