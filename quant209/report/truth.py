"""Read analyst ground-truth from the lab's quantification workbooks.

Each analysis-run folder holds one or more analyst-filled Excel workbooks, named for the
samples they cover, plus an instrument sequence-log workbook that must be skipped. The
per-workbook layout:

    - One worksheet per sample, named by sample id. Row 15 is the
      header (mass transition | ... | Cal RT | RT | Area | IS | RRF | Mass (ng) |
      Corrected). Row 16+ is one row per analyte, column C = analyte label.
      Columns used: C=label (3), P=area (16), J=mass in ng (10). These are the
      analyst's manually integrated peak area and the workbook-calculated congener
      mass, before recovery correction. A handful of sheets (``SamplesList_Standards``,
      ``cal_info``, ``cal``, ``ref``, ``QuickView``, ``Summary``, ``Uncorrected_Mass``,
      ``SS_Corrected_Mass``, ``Profile``) are workbook scaffolding, not per-sample data,
      and are excluded.
    - A ``Summary`` sheet holding a "Recoveries" block: a header row (a "Recoveries"
      label cell, a "REF" column, then one column per sample) followed by exactly ten
      data rows, one per homolog in mono->deca order, giving each homolog's 13C
      surrogate recovery fraction for that sample. Row *labels* in this block are
      unreliable analyst free text (one workbook in the archive labels the tri-CB row
      "13C 31" instead of the canonical "13C 28"), so homolog identity is read from
      row position (matching ``core.standards.HOMOLOG_SURROGATE``'s mono->deca
      iteration order), not from the label text.
    - A ``cal`` sheet: the calibration injection this workbook quantified against. Row 5
      is the header (mass transition | congener | RT | Area | IS mass/area | mass | RRF |
      ... | Peak | RT | Area | Height). Row 6+ is one row per analyte. Columns used:
      C=label (3), E=area (5), H=RRF (8). These are the analyst's manually integrated
      *calibration* peak area and the RRF derived from it,
      ``RRF = (cal_mass / cal_area) / (cal_mass_IS / cal_area_IS)``, the formula
      ``core.quantify.build_calibration`` implements. Verified against the sheet: for
      congener 1 in one workbook, ``29.0606 / 115270.09`` over the d30
      response reproduces the sheet's 2.20945 exactly.
    - A ``ref`` sheet: the reference injection whose 13C surrogate masses define 100%
      recovery for that workbook's samples. Two blocks. A ten-row summary block at the top
      left of the sheet holds each surrogate's mass (ng) in that reference injection. It
      is a header row carrying a "Recoveries" label cell and a "Mass" column, then one row
      per homolog in mono->deca order. Below it, from row 15, a per-analyte table in the
      same shape as a sample sheet but at its own column offsets. Only the summary block is
      read: its "Mass" column is cell-for-cell the per-analyte table's "Mass (ng)" for the
      surrogate rows, and it is already keyed by homolog position. As in the ``Summary``
      block, row labels are analyst free text ("13C 31" for the tri-CB surrogate) and
      homolog identity comes from row position.

**The ``cal`` sheet's "Source:" cell is not the calibration injection's name.** Across the
archive it names a sample range belonging to some *other* workbook, often one from a
different project entirely. It is a breadcrumb naming the workbook this one was copied from,
not a data file in this run. The cal sheet is therefore keyed by the workbook that holds it
and nothing else. Matching a workbook to the injection it actually used is left to a caller
that can compare the areas, which ``tools/bench_quant.py`` does.

Orientation convention: every DataFrame below is indexed by analyte label or homolog on
axis 0 (rows). Axis 1 (columns) carries the sample id, or, for the ``cal``/``ref`` frames,
the workbook that holds the sheet. That orientation lets a caller look up one
analyte or homolog's values across every sample with a single ``.loc[label]``.

Why the ``cal``/``ref`` frames are keyed by *workbook* and not by injection: each workbook
holds exactly one ``cal`` sheet and one ``ref`` sheet, and neither names the data file it came
from. Inventing an injection id here would be an inference dressed up as ground truth. What
the workbook does say is which samples it covers, so ``workbook_samples`` is exposed alongside
and a caller that knows the run's own reference bracketing (``RunResult.reference_for``) can
resolve workbook -> reference injection from evidence rather than from a guess.

Masses here are *uncorrected* (workbook column J, "Mass (ng)"). They are the RRF-based
mass computed straight from the manually integrated area, before surrogate-recovery
correction. Recoveries are exposed separately so a caller can check each stage of the
calculation chain independently, from area to uncorrected mass to recovery correction,
rather than only the final corrected number. The analyst's own corrected mass (column K)
is reproducible as ``mass / recovery`` when ``recovery < 1``, following the EPA-1668
convention of correcting for losses and never deflating for apparent gains.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import openpyxl
import pandas as pd

from quant209.core.standards import HOMOLOG_SURROGATE

# Per-sample sheet geometry is *discovered*, not assumed. See ``_discover_layout``.
# Two layouts occur in this archive and they differ by one column and three rows:
#
#   modern (2022-10 onward)  header row 15: B=mass transition C=label ... G=Area
#                                           H=IS I=RRF J=Mass (ng) K=Corrected Mass
#   early (pre-2022-10)       header row 12: A=mass transition B=label ... F=Area
#                                           H=RRF I=Mass (ng) J=Corrected Mass
#
# Both also repeat the manually integrated area in a right-hand "Peak" block (modern
# column P). That block and the calculation block's "Area" are the same number, because P
# is the measurement and G is a formula referencing it, verified cell for cell against a
# real workbook. Reading the *first* "Area" at or after the transition column therefore
# works for both layouts and needs no per-layout special case.
_LAYOUT_SCAN_ROWS = 30  # generous: the header sits at row 12 or 15 in every workbook seen
_LAYOUT_MAX_COL = 20

# Sheets that are workbook scaffolding, not per-sample analyte data. Compared after
# ``_normalised`` so that "Sample List_Standard Info" (early workbooks) and
# "SamplesList_Standards" (modern) are one entry, and spacing variants cannot slip through.
_NON_SAMPLE = {
    "SamplesList_Standards", "Sample List_Standard Info", "cal_info", "cal", "ref",
    "QuickView", "Summary", "Uncorrected_Mass", "SS_Corrected_Mass", "Profile",
    "Aroclor Analysis", "Aroclor Analysis 1016, 1254", "Non Aroclor PCBs + PCB3",
}

def _normalised(sheet_name: str) -> str:
    """Sheet name reduced for comparison: case-folded, spaces and underscores removed.

    Analysts name the same sheet ``cal``, ``Cal 01`` and ``cal_01`` across workbooks.
    """
    return sheet_name.strip().casefold().replace(" ", "").replace("_", "")


def _is_cal_sheet(sheet_name: str) -> bool:
    """True for the calibration sheet under any of its observed names (``cal``, ``Cal 01``)."""
    return re.fullmatch(r"cal\d*", _normalised(sheet_name)) is not None


# mono..deca, the row order every Recoveries block uses.
_HOMOLOG_ORDER = tuple(HOMOLOG_SURROGATE)

#: ``_NON_SAMPLE`` reduced by ``_normalised`` for comparison.
_NON_SAMPLE_NORMALISED = frozenset()  # populated below, after _normalised exists

# Bounds for the two bounded iter_rows() scans _read_recoveries_block does over the
# Summary sheet: generous enough for every observed layout, still small enough that a
# single-pass read stays cheap. Explicit constants rather than ws.max_row/max_column,
# which read_only worksheets can report inflated for sheets with sparse formatting.
_SUMMARY_HEADER_SCAN_ROWS = 15
# Bounds for the per-sample inline Recoveries block (_read_inline_recoveries): it sits
# in the first rows, right of the analyte table (column N in every workbook seen).
_RECOVERY_SCAN_ROWS = 20
_RECOVERY_MAX_COL = 24
_SUMMARY_MAX_COL = 30

# ``cal`` sheet: column indices (1-based) and the first data row. The header sits at row 5 and
# is *located* rather than assumed (see ``_read_cal_sheet``), so a workbook that gained a row
# above the table is read correctly instead of silently one row out.
_CAL_COL_LABEL, _CAL_COL_AREA, _CAL_COL_RRF = 3, 5, 8
_CAL_MAX_COL = 8
_CAL_HEADER_SCAN_ROWS = 12

# ``ref`` sheet: bound for the scan that finds the surrogate summary block's header. The block
# sits in the first few rows. The per-analyte table below it starts at row 15 and is not read.
_REF_HEADER_SCAN_ROWS = 8
_REF_MAX_COL = 24


@dataclass(eq=False)
class RunTruth:
    """Analyst ground truth for one analysis-run folder.

    areas: manually integrated peak area per (analyte label, sample id).
    masses: analyst-calculated, pre-recovery-correction congener mass (ng) per
        (analyte label, sample id), with the same index and columns as ``areas``.
    recoveries: fractional (~0-1.5) 13C-surrogate recovery per (homolog, sample id).
    cal_areas: manually integrated *calibration*-injection peak area per (analyte label,
        workbook stem). These are the areas every RRF in that workbook was derived from.
    cal_rrfs: the analyst's relative response factor per (analyte label, workbook stem),
        same index/columns as ``cal_areas``.
    ref_surrogate_masses: 13C-surrogate mass (ng) in the reference injection, per (homolog,
        workbook stem). This is the denominator of every recovery in that workbook.
    workbook_samples: workbook stem -> the sample sheets it holds, so a caller can resolve
        which injection a ``cal``/``ref`` column belongs to (see the module docstring).

    A workbook whose ``cal``/``ref`` sheets were never filled in still contributes a column.
    This happens when an issue was superseded and re-quantified under a later filename. Its
    RRF and reference mass cells hold ``#DIV/0!`` and read as NaN. Its cal *area* cells hold
    a literal 0 and read as 0.0, because that is what the sheet says. The column is kept for
    the same reason an unquantified sample column is kept in ``masses``. "The analyst
    produced no number here" is itself a fact a caller needs to see. Callers screen such columns out on the RRF frame
    (``.notna().any()``) rather than on the area frame, whose zeros are not missing values.

    ``eq=False``: the generated ``__eq__`` compares fields with ``==``, which on
    ``pd.DataFrame`` fields returns another DataFrame rather than a bool and raises
    ("truth value of a DataFrame is ambiguous"). ``Chromatogram`` in ``core/mzml.py``
    carries the same defect and the same fix. Identity-based equality is what is wanted here.
    """

    areas: pd.DataFrame
    masses: pd.DataFrame
    recoveries: pd.DataFrame
    cal_areas: pd.DataFrame = field(default_factory=pd.DataFrame)
    cal_rrfs: pd.DataFrame = field(default_factory=pd.DataFrame)
    ref_surrogate_masses: pd.DataFrame = field(default_factory=pd.DataFrame)
    workbook_samples: dict[str, list[str]] = field(default_factory=dict)
    #: The calibration-standard mass (ng) each workbook's ``cal_info`` sheet enters per
    #: (analyte label, workbook stem). It is what every RRF in that workbook divides by, so a
    #: workbook that enters a different mass for a co-eluting group reports group masses on a
    #: different basis (see ``report.validation.rebase_group_masses``).
    cal_masses: pd.DataFrame = field(default_factory=pd.DataFrame)
    #: Sample id -> the workbook stem its sheet was finally read from, after re-issues.
    sample_book: dict[str, str] = field(default_factory=dict)


_NON_SAMPLE_NORMALISED = frozenset(_normalised(n) for n in _NON_SAMPLE)


def _num(value: object) -> float | None:
    """Coerce a workbook cell value to float.

    Excel error strings such as '#DIV/0!' and blank cells become None rather than raising.
    """
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


@dataclass(frozen=True)
class SheetLayout:
    """Where one per-sample sheet keeps the columns this reader needs (1-based)."""

    header_row: int
    label: int
    area: int
    mass: int


def _discover_layout(ws) -> SheetLayout | None:
    """Locate the analyte table on ``ws``, or None if this is not a per-sample sheet.

    The header row is the one carrying "mass transition". The analyte label sits in the
    column immediately right of it in both known layouts, and "Area" and "Mass (ng)" are
    found by their own header text rather than by a fixed offset. Discovering the geometry
    instead of hardcoding it is what lets one reader serve workbooks written years apart by
    different analysts. It also fails loudly, returning None so the sheet is skipped and
    reported, rather than silently reading the wrong column.

    Uses bounded ``iter_rows`` fetches rather than ``ws.cell()``: on a ``read_only``
    worksheet each ``.cell()`` re-parses the sheet XML from row 1, so a per-cell loop is
    quadratic over a workbook.
    """
    if _normalised(ws.title) in _NON_SAMPLE_NORMALISED:
        return None
    for offset, row in enumerate(
        ws.iter_rows(min_row=1, max_row=_LAYOUT_SCAN_ROWS, max_col=_LAYOUT_MAX_COL,
                     values_only=True)
    ):
        header = {
            i: v.strip().casefold()
            for i, v in enumerate(row, start=1)
            if isinstance(v, str) and v.strip()
        }
        transition = next((i for i, t in header.items() if t == "mass transition"), None)
        if transition is None:
            continue
        area = next((i for i, t in sorted(header.items()) if t == "area" and i > transition), None)
        mass = next((i for i, t in sorted(header.items()) if t.startswith("mass (ng)")), None)
        if area is None or mass is None:
            continue
        return SheetLayout(header_row=offset + 1, label=transition + 1, area=area, mass=mass)
    return None


def _read_sample_sheet(ws, layout: SheetLayout) -> dict[str, tuple[float | None, float | None]]:
    """{analyte label: (area, mass)} for every labelled row in one per-sample sheet.

    One ``iter_rows`` pass. See ``_discover_layout`` for why per-cell ``.cell()`` access
    must be avoided on a ``read_only`` worksheet.
    """
    out: dict[str, tuple[float | None, float | None]] = {}
    max_col = max(layout.label, layout.area, layout.mass)
    for row in ws.iter_rows(min_row=layout.header_row + 1, max_col=max_col, values_only=True):
        label = row[layout.label - 1]
        if label is None:
            continue
        out[str(label)] = (_num(row[layout.area - 1]), _num(row[layout.mass - 1]))
    return out


def _read_inline_recoveries(ws) -> dict[str, float]:
    """{homolog: recovery} from the "Recoveries" block on one per-sample sheet.

    Every layout in this archive repeats each sample's ten surrogate recoveries in a small
    block to the right of the analyte table ("Recoveries" label, then ten rows of
    ``label | recovery | mass``). Reading it per sheet is preferred to reading the
    ``Summary`` sheet's grid. The value is already attached to the sample whose sheet it is,
    so no sample-name column matching is needed. It also works for the early workbooks,
    which have no ``Summary`` sheet at all.

    Homolog identity comes from row *position* (mono..deca), not from the row label: the
    labels are analyst free text and disagree between workbooks. The tri-CB surrogate is
    written "13C 28" in some and "13C 31" in others.
    """
    block = _find_header_cell(ws, "recoveries", _RECOVERY_SCAN_ROWS, _RECOVERY_MAX_COL)
    if block is None:
        return {}
    # _find_header_cell reports the column 0-based. iter_rows wants it 1-based.
    header_row, label_col0 = block
    label_col = label_col0 + 1
    rows = list(
        ws.iter_rows(
            min_row=header_row + 1,
            max_row=header_row + len(_HOMOLOG_ORDER),
            min_col=label_col,
            max_col=label_col + 1,
            values_only=True,
        )
    )
    out: dict[str, float] = {}
    for homolog, row in zip(_HOMOLOG_ORDER, rows):
        value = _num(row[1]) if len(row) > 1 else None
        if value is not None:
            out[homolog] = value
    return out


def _read_recoveries_block(ws) -> dict[str, dict[str, float]]:
    """{homolog: {sample id: recovery}} from a Summary sheet's "Recoveries" block.

    Two bounded ``iter_rows`` passes. See ``_looks_like_sample_sheet`` for why per-cell
    ``.cell()`` access must be avoided on a ``read_only`` worksheet.
    """
    header_scan = ws.iter_rows(
        min_row=1, max_row=_SUMMARY_HEADER_SCAN_ROWS, max_col=6, values_only=True
    )
    header_row: int | None = None
    header_col: int | None = None  # 0-based index of the "Recoveries" label cell
    for i, row in enumerate(header_scan, start=1):
        for j, value in enumerate(row):
            if isinstance(value, str) and value.strip().lower() == "recoveries":
                header_row, header_col = i, j
                break
        if header_row is not None:
            break
    if header_row is None or header_col is None:
        return {}

    block_rows = list(
        ws.iter_rows(
            min_row=header_row,
            max_row=header_row + len(_HOMOLOG_ORDER),
            max_col=_SUMMARY_MAX_COL,
            values_only=True,
        )
    )
    header = block_rows[0]

    # Sample columns: the contiguous run of non-empty string cells right after the label
    # column, skipping a leading "REF" column. The reference injection is not a sample, and
    # its recovery is 1 by definition, since recovery is measured relative to it.
    samples: list[tuple[int, str]] = []  # (0-based column index, sample id)
    is_first = True
    idx = header_col + 1
    while idx < len(header):
        value = header[idx]
        if not isinstance(value, str) or not value.strip():
            break
        name = value.strip()
        if not (is_first and name.upper() == "REF"):
            samples.append((idx, name))
        is_first = False
        idx += 1

    out: dict[str, dict[str, float]] = {}
    for i, homolog in enumerate(_HOMOLOG_ORDER):
        data_row = block_rows[i + 1]
        row_values: dict[str, float] = {}
        for col_idx, sample in samples:
            value = _num(data_row[col_idx]) if col_idx < len(data_row) else None
            if value is not None:
                row_values[sample] = value
        out[homolog] = row_values
    return out


def _find_header_cell(ws, text: str, max_row: int, max_col: int) -> tuple[int, int] | None:
    """(1-based row, 0-based column) of the first cell equal to ``text``, case-insensitively.

    One bounded ``iter_rows`` pass. See ``_looks_like_sample_sheet`` for why per-cell
    ``.cell()`` access must be avoided on a ``read_only`` worksheet. Locating the header rather
    than hardcoding its row number is what keeps a workbook that gained a banner row above the
    table from being read one row out. On a table of numbers with no other structure, that
    would produce a plausible-looking frame of wrong values rather than an error.
    """
    scan = ws.iter_rows(min_row=1, max_row=max_row, max_col=max_col, values_only=True)
    wanted = text.strip().lower()
    for i, row in enumerate(scan, start=1):
        for j, value in enumerate(row):
            if isinstance(value, str) and value.strip().lower() == wanted:
                return i, j
    return None


def _read_cal_sheet(ws) -> dict[str, tuple[float | None, float | None]]:
    """{analyte label: (calibration area, RRF)} from one workbook's ``cal`` sheet.

    Returns ``{}`` when the header cannot be found, which is the honest answer for a sheet
    whose layout this reader does not recognise: an empty column in a benchmark is visible,
    whereas a column read at the wrong offsets is not.
    """
    header = _find_header_cell(ws, "RRF", _CAL_HEADER_SCAN_ROWS, _CAL_MAX_COL)
    if header is None:
        return {}
    header_row, header_col = header
    if header_col != _CAL_COL_RRF - 1:  # the layout moved, so refuse rather than mis-read
        return {}
    out: dict[str, tuple[float | None, float | None]] = {}
    for row in ws.iter_rows(min_row=header_row + 1, max_col=_CAL_MAX_COL, values_only=True):
        label = row[_CAL_COL_LABEL - 1]
        if label is None:
            continue
        out[str(label)] = (_num(row[_CAL_COL_AREA - 1]), _num(row[_CAL_COL_RRF - 1]))
    return out


def _read_ref_surrogates(ws) -> dict[str, float | None]:
    """{homolog: surrogate mass (ng)} from a ``ref`` sheet's ten-row summary block.

    The block's header row carries a "Recoveries" label cell and, further right, a "Mass"
    column. The ten rows under it are the homologs in mono->deca order. Homolog identity
    comes from row *position* and never from the row label, for the reason documented on
    ``_read_recoveries_block``. The labels are analyst free text, and one of them names the
    tri-CB surrogate "13C 31" rather than the canonical "13C 28".
    """
    header = _find_header_cell(ws, "Recoveries", _REF_HEADER_SCAN_ROWS, _REF_MAX_COL)
    if header is None:
        return {}
    header_row, header_col = header
    block_rows = list(
        ws.iter_rows(
            min_row=header_row,
            max_row=header_row + len(_HOMOLOG_ORDER),
            max_col=_REF_MAX_COL,
            values_only=True,
        )
    )
    # The mass column is located from the header row rather than assumed to be a fixed offset
    # from the "Recoveries" cell. The recovery column between them is a constant 1 and would
    # read as a plausible "mass" if the offset were ever wrong by one.
    mass_col = next(
        (
            j
            for j, value in enumerate(block_rows[0])
            if j > header_col and isinstance(value, str) and value.strip().lower() == "mass"
        ),
        None,
    )
    if mass_col is None:
        return {}
    # Every homolog is reported, ``None`` included: a block found but full of ``#DIV/0!`` (a
    # superseded workbook that was never re-quantified) is a different fact from a workbook
    # with no ``ref`` sheet at all, and only the former should occupy a column.
    return {
        homolog: (_num(block_rows[i + 1][mass_col]) if mass_col < len(block_rows[i + 1]) else None)
        for i, homolog in enumerate(_HOMOLOG_ORDER)
    }


def _find_quant_workbooks(run_folder: Path) -> list[Path]:
    """Analyst quantification workbooks in a run folder: every .xlsx except Excel's
    temp lock files (``~$*``) and the instrument sequence-log workbook."""
    return [
        p
        for p in sorted(run_folder.glob("*.xlsx"))
        if not p.name.startswith("~$") and "sequence" not in p.name.lower()
    ]


#: Words that mark a workbook as superseded or experimental (``OLD``, ``test...``, ``Copy``).
#: Such a workbook never overrides one without the mark.
_SUPERSEDED = re.compile(r"(?i)(?:^|[^a-z])(old|test[a-z0-9]*|copy)(?:[^a-z]|$)")
#: A YYYYMMDD date anywhere in a workbook's name, the analyst's own version stamp. Matched
#: without word boundaries so that a stamp run into a tag (``_tag20220808``) is still read.
_DATE_TOKEN = re.compile(r"20\d{2}(?:0[1-9]|1[0-2])(?:0[1-9]|[12]\d|3[01])")


def _version_key(path: Path) -> tuple[bool, str, str]:
    """Sort key putting a run's workbooks in the order their re-issues were made.

    Re-issues of a workbook cover the same samples, and the last one read wins. Ordering by
    filename alone made that the one whose name sorts last, which need not be the latest: a
    tag run into a name sorts after a later date stamp, because a letter sorts after a digit.
    Workbooks marked as superseded come first, so anything else overrides them. Then
    the latest date stamp in the name, or the file's own modification date when the name has
    none. Then the name, which keeps the old behaviour wherever the stamps tie.
    """
    stamps = _DATE_TOKEN.findall(path.stem)
    try:
        when = max(stamps) if stamps else datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y%m%d")
    except OSError:
        when = ""
    return (_SUPERSEDED.search(path.stem) is None, when, path.name)


def _read_cal_info(ws) -> dict[str, float | None]:
    """{analyte label: calibration-standard mass (ng)} from a workbook's ``cal_info`` sheet."""
    names = ("congener", "label", "analyte")
    for i, row in enumerate(ws.iter_rows(min_row=1, max_row=_CAL_HEADER_SCAN_ROWS,
                                          max_col=_CAL_MAX_COL, values_only=True)):
        cells = [str(c).strip().lower() if isinstance(c, str) else "" for c in row]
        if "mass" in cells and any(c in names for c in cells):
            label_col = next(j for j, c in enumerate(cells) if c in names)
            mass_col = cells.index("mass")
            out: dict[str, float | None] = {}
            for r in ws.iter_rows(min_row=i + 2, max_col=max(label_col, mass_col) + 1,
                                  values_only=True):
                if r[label_col] is not None:
                    out[str(r[label_col])] = _num(r[mass_col])
            return out
    return {}


def _warn_skipped(xlsx: Path, exc: Exception) -> None:
    """Warn that a workbook was skipped.

    A malformed or unreadable workbook should not abort the whole run, but it must not be
    dropped silently either. This is the ground-truth source every validation benchmark
    trusts. Shared by both catch sites in ``read_run_truth``, the load failure and the
    per-workbook parse failure, so the reason for the broad ``except Exception`` there is
    stated once rather than twice.
    """
    warnings.warn(f"skipping unreadable workbook {xlsx.name}: {exc!r}", stacklevel=3)


def read_run_truth(run_folder: str | Path) -> RunTruth:
    """Merge analyst ground truth from every quantification workbook in a run folder.

    Workbooks are read in the order their re-issues were made (``_version_key``), and a
    later one overrides an earlier one for the same
    (analyte/homolog, sample) cell, so a corrected re-issue wins over an earlier version
    covering the same samples. A malformed or unreadable workbook is skipped with a warning
    rather than aborting the whole run. See ``_warn_skipped``.
    """
    run_folder = Path(run_folder)
    area_records: dict[str, dict[str, float | None]] = {}
    mass_records: dict[str, dict[str, float | None]] = {}
    recovery_records: dict[str, dict[str, float]] = {}
    cal_area_records: dict[str, dict[str, float | None]] = {}
    cal_rrf_records: dict[str, dict[str, float | None]] = {}
    ref_records: dict[str, dict[str, float | None]] = {}
    workbook_samples: dict[str, list[str]] = {}
    cal_mass_records: dict[str, dict[str, float | None]] = {}
    sample_book: dict[str, str] = {}

    for xlsx in sorted(_find_quant_workbooks(run_folder), key=_version_key):
        try:
            wb = openpyxl.load_workbook(xlsx, data_only=True, read_only=True)
        except Exception as exc:  # noqa: BLE001 - see _warn_skipped
            _warn_skipped(xlsx, exc)
            continue
        book = xlsx.stem
        try:
            for ws in wb.worksheets:
                if _normalised(ws.title) == "calinfo":
                    for label, mass in _read_cal_info(ws).items():
                        cal_mass_records.setdefault(label, {})[book] = mass
                    continue
                if _is_cal_sheet(ws.title):
                    for label, (cal_area, rrf) in _read_cal_sheet(ws).items():
                        cal_area_records.setdefault(label, {})[book] = cal_area
                        cal_rrf_records.setdefault(label, {})[book] = rrf
                    continue
                layout = _discover_layout(ws)
                if layout is None:
                    continue
                sample = ws.title
                measured = _read_sample_sheet(ws, layout)
                # A sheet whose every area is blank or zero is an unfilled template or a
                # stale copy the analyst abandoned. Both occur in these workbooks: the
                # numbered "9" and "10" sheets, and re-issues that kept the previous batch's
                # sample names over #DIV/0! cells. Neither was quantified, so admitting one
                # would invent a column of zeros for the benchmark to score against.
                if not any(area for area, _ in measured.values()):
                    continue
                workbook_samples.setdefault(book, []).append(sample)
                sample_book[sample] = book
                for label, (area, mass) in measured.items():
                    area_records.setdefault(label, {})[sample] = area
                    mass_records.setdefault(label, {})[sample] = mass
                # Each sheet carries its own sample's recoveries. The Summary grid below
                # is consulted only for samples no sheet supplied.
                for homolog, value in _read_inline_recoveries(ws).items():
                    recovery_records.setdefault(homolog, {})[sample] = value
            if "Summary" in wb.sheetnames:
                for homolog, by_sample in _read_recoveries_block(wb["Summary"]).items():
                    for sample, value in by_sample.items():
                        recovery_records.setdefault(homolog, {}).setdefault(sample, value)
            if "ref" in wb.sheetnames:
                for homolog, mass in _read_ref_surrogates(wb["ref"]).items():
                    ref_records.setdefault(homolog, {})[book] = mass
        except Exception as exc:  # noqa: BLE001 - see _warn_skipped
            _warn_skipped(xlsx, exc)
            continue
        finally:
            wb.close()

    areas = pd.DataFrame(area_records).T.sort_index().sort_index(axis=1)
    masses = pd.DataFrame(mass_records).T.reindex(index=areas.index, columns=areas.columns)
    recoveries = pd.DataFrame(recovery_records).T.sort_index().sort_index(axis=1)
    cal_areas = pd.DataFrame(cal_area_records).T.sort_index().sort_index(axis=1)
    cal_rrfs = pd.DataFrame(cal_rrf_records).T.reindex(
        index=cal_areas.index, columns=cal_areas.columns
    )
    # mono->deca row order, matching ``recoveries``' own homolog ordering rather than the
    # alphabetical sort a bare DataFrame construction would impose.
    ref_masses = pd.DataFrame(ref_records).T.reindex(index=list(_HOMOLOG_ORDER))

    return RunTruth(
        areas=areas,
        masses=masses,
        recoveries=recoveries,
        cal_areas=cal_areas,
        cal_rrfs=cal_rrfs,
        ref_surrogate_masses=ref_masses.sort_index(axis=1),
        workbook_samples={k: sorted(v) for k, v in sorted(workbook_samples.items())},
        cal_masses=pd.DataFrame(cal_mass_records).T.sort_index().sort_index(axis=1),
        sample_book=sample_book,
    )
