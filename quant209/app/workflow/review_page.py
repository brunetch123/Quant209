"""Step 5, Review: the run verdict, the chromatograms, and every flag behind them.

What this screen is for
------------------------
Everything before this page is machinery. This is where a human decides whether a run's
numbers may leave the building.

The dashboard is built around one claim, that **every number on it can be traced to a
picture of the trace it came from**, and one reading order, left to right along the tab bar:

============  =============================================================================
Overview      the verdict, in one banner and five numbers. A triage screen: no charts, and
              a table that lists *only* what is wrong. If it is empty, stop reading.
Samples       the audit centrepiece. Any file, any transition, every analyte's integration
              drawn over the acquired trace (see :class:`~quant209.app.charts.
              ChromatogramView` for what each mark means).
Calibration   the response factors, per injection and combined, with the mis-integration
              guard's refusals stated rather than left as blank cells.
Recoveries    homolog x sample, as a heatmap with the 50-150% acceptance band drawn on it.
QA/QC         every flag event in the run, filterable, exportable.
============  =============================================================================

The import discipline, which is not optional
----------------------------------------------
**This module imports nothing heavy at module scope and must stay that way.** The workflow
package's page registry resolves module names lazily. The moment this file grows an
``import pandas`` or an ``import pyqtgraph`` at the top, every one of those costs lands on
whichever navigation first touches it. A stray ``from ... import`` anywhere on the startup
path would also drag it into the three-second budget. So:

* pyqtgraph arrives via :mod:`quant209.app.charts`, imported when a tab that draws is first
  *built*. An analyst who only reads the Overview never pays for it at all.
* ``core.mzml`` (numpy) arrives inside :class:`ChromatogramCache.load`.
* ``core.quantify`` (numpy) arrives inside the Calibration tab's builder.
* the science constants this page states on screen (the recovery band, the homolog order,
  the surrogate assignment) are **duplicated here as literals** rather than imported from
  ``core.qaqc``/``core.standards``, and ``tests/test_app_review.py`` pins each copy to its
  original so the two cannot drift.

Tabs build on first click
--------------------------
Only Overview is built when the page is activated. The other four are QWidget shells until
the analyst opens them, because building the Samples tab means constructing a plot, reading
an mzML, and laying out 185 rows. That work is invisible when it happens on a click and very
visible when it happens on a page transition.

Switching files stays warm
---------------------------
Re-plotting reads the file's chromatograms back from disk, which is the one expensive thing
this screen does (~1 s for a 30 MB mzML). :class:`ChromatogramCache` keeps the last four
files, matching how the screen is actually used: comparing a sample against its blank, its
reference, and the calibration injection. It holds the memory flat on a host with 2 GB. A
cold read says so on screen before it blocks, rather than freezing silently.
"""

from __future__ import annotations

import csv
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from PySide6.QtCore import QAbstractTableModel, QModelIndex, QRectF, QSize, Qt
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QComboBox,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QSizePolicy,
    QSplitter,
    QStyle,
    QStyledItemDelegate,
    QStyleOptionViewItem,
    QTableView,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from quant209.app import theme
from quant209.app.state import STEP_REVIEW
from quant209.app.workflow.base import (
    ROLE_LABELS,
    Banner,
    Card,
    EmptyState,
    WorkflowPage,
    ghost_button,
    hint,
    meta_label,
    scroll_wrap,
    secondary_button,
)

# --------------------------------------------------------------------------------------
# Science constants, duplicated rather than imported. See the module docstring.
# --------------------------------------------------------------------------------------

#: ``core.pipeline.STATUS_OK``.
STATUS_OK = "OK"

#: ``core.qaqc.QaqcParams.recovery_low`` / ``recovery_high``: the acceptance band, as
#: fractions. The heatmap draws it, the verdict counts breaches of it, and the QA/QC tab
#: names it, so it is stated once here.
RECOVERY_LOW = 0.50
RECOVERY_HIGH = 1.50

#: ``core.standards.HOMOLOG_SURROGATE``: which labelled surrogate tracks each homolog's
#: recovery, in chlorination order. The heatmap's row order *is* this dict's key order.
HOMOLOG_SURROGATE: dict[str, str] = {
    "mono": "13C 3",
    "di": "13C 15",
    "tri": "13C 28",
    "tetra": "13C 52",
    "penta": "13C 118",
    "hexa": "13C 153",
    "hepta": "13C 180",
    "octa": "13C 194",
    "nona": "13C 206",
    "deca": "13C 209",
}

#: The two ``core.integrate.KNOWN_FLAGS`` that say "this area is imprecise". Counted as
#: *flagged integrations*. ``absent`` is not, because it is a presence verdict rather than
#: a defect in an integration (an absent congener is usually the right answer). ``clustered``
#: is a statement about which baseline was drawn.
QUALITY_FLAGS: tuple[str, ...] = ("low_sn", "narrow")

#: The share of a file's integrated peaks that may carry a quality flag before the run is
#: worth a second look. A *fraction* rather than a count: a run is 185 analytes in every one
#: of twenty-odd files, so any absolute number is either always tripped or never tripped
#: depending only on how many files were injected.
FLAGGED_FRACTION_LIMIT = 0.05

#: How far one congener's response factor may vary *between* this run's calibration
#: injections before the Overview says so, as a relative standard deviation over the
#: injections that produced a factor at all.
#:
#: Tied to a number this project already uses elsewhere, rather than chosen independently:
#: the validation scores every mass against a ±25% band, so a response factor whose own
#: injections disagree by more than that band can carry a mass outside it on its own, with
#: nothing wrong anywhere else in the chain. It is also a failure the archive shows twice:
#: on one run a calibration anomaly accounts for 40% of the mass misses, and on another a
#: single calibration injection leaves 19 response factors outside the same band.
RRF_SPREAD_LIMIT = 0.25

#: How many files' chromatograms are held in memory at once. Four is the working set of the
#: question this screen answers: a sample, its blank, its reference, the calibration.
CACHE_FILES = 4

# --------------------------------------------------------------------------------------
# Tab and table shape
# --------------------------------------------------------------------------------------

TAB_OVERVIEW, TAB_SAMPLES, TAB_CALIBRATION, TAB_RECOVERIES, TAB_QAQC = range(5)
TAB_TITLES: tuple[str, ...] = ("Overview", "Samples", "Calibration", "Recoveries", "QA/QC")

#: The five file groups in the Samples rail, in the order an analyst looks for them.
RAIL_GROUPS: tuple[tuple[str, str], ...] = (
    ("sample_ids", "SAMPLES"),
    ("blank_ids", "BLANKS"),
    ("ref_ids", "REFERENCE"),
    ("cal_ids", "CALIBRATION"),
)

#: ``index.data(FLAGS_ROLE)`` -> ``list[str]``, painted as chips by :class:`FlagChipDelegate`.
FLAGS_ROLE = int(Qt.ItemDataRole.UserRole) + 1
#: ``index.data(VALUE_ROLE)`` -> the raw ``float`` behind a formatted cell (the heatmap).
VALUE_ROLE = int(Qt.ItemDataRole.UserRole) + 2
#: ``index.data(OUT_OF_BAND_ROLE)`` -> ``True`` for a recovery outside the acceptance band.
OUT_OF_BAND_ROLE = int(Qt.ItemDataRole.UserRole) + 3

#: Chip tone per flag kind. Amber for anything that weakens a number, blue for the three
#: that only describe how it was measured, grey for a censoring decision that worked as
#: designed.
FLAG_TONES: dict[str, str] = {
    "absent": "neutral",
    "low_sn": "warning",
    "narrow": "warning",
    "clustered": "info",
    "shoulder": "info",
    "forced": "info",
    "displaced": "info",
    "apportioned": "info",
    "unresolved": "info",
    "ghost": "info",
    "realigned": "info",
    "censored": "neutral",
    "recovery": "error",
    "missing_is": "error",
}

#: The QA/QC tab's filter chips, in order. ``All`` is not a kind. It clears the others.
FLAG_KINDS: tuple[str, ...] = (
    "recovery", "missing_is", "absent", "low_sn", "narrow", "clustered", "shoulder", "forced", "displaced",
    "apportioned", "unresolved", "ghost", "realigned", "censored",
)

#: How a flag kind is *written* on screen. The stored kind is unchanged: it is what the CSV
#: export writes and what every filter compares. But ``low_sn`` is an identifier, not a word,
#: and a screen that shows an analyst the engine's variable names is a screen that was never
#: read by one. Anything absent from this map falls back to its own kind.
FLAG_LABELS: dict[str, str] = {"low_sn": "low S/N", "missing_is": "missing IS"}

#: The same rule for ``core.integrate.KNOWN_METHODS``, the peak table's "Method" column.
METHOD_LABELS: dict[str, str] = {
    "valley_drop": "valley drop",
    "rt_window": "RT window",
    "tangent_skim": "tangent skim",
    "deconvolved": "deconvolved",
}


def flag_text(kind: str) -> str:
    """One flag kind, as it is written on screen."""
    return FLAG_LABELS.get(kind, kind)


ROW_HEIGHT = 30
#: Narrowest a heatmap column may become before the grid scrolls instead of shrinking.
#: Measured rather than chosen, and against the *header* rather than the cell: ``"107%"`` is
#: 34 px at the body size but the column is titled with a file id, and a five-character id at
#: the header's 9 pt semibold is 46. Below this the ids elide to a leading ellipsis, not
#: merely hard to read but a different sample.
HEATMAP_MIN_COLUMN = 54


def _is_number(value: Any) -> bool:
    """True for a finite real number. NaN fails its own equality test, which is the point."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return number == number and number not in (float("inf"), float("-inf"))


# --------------------------------------------------------------------------------------
# The run summary: everything the Overview tab states, computed once
# --------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------
# Severity: what "worst first" means on the Overview's table
# --------------------------------------------------------------------------------------
#
# The table sorts on :attr:`FailingItem.weight`, and the weights are **banded** rather than
# free so the order is a stated policy instead of an accident of how each class's number
# happens to scale. One band per class, a hundred apart, worst first:
#
#   ==========  ==================================================================
#   900 .. 600  the run's numbers may be wrong wholesale: a run-level guard failed
#   500 .. 400  a congener is unquantifiable, or quantified against an unstable factor
#          300  a whole file reported nothing
#   200 .. 100  one file x one homolog is in doubt
#     0 ..  99  an advisory count
#   ==========  ==================================================================
#
# Within a band an item adds at most 99, so a severe item never outranks a whole class above
# it. The ordering *between* classes is the policy, and the magnitude only breaks ties inside
# one.
W_FIT = 900.0
W_NO_PEAKS = 800.0
W_NO_REFERENCE = 700.0
W_NO_BLANKS = 600.0
W_RRF_REFUSED = 500.0
W_RRF_SPREAD = 400.0
W_NOTHING_ABOVE_LOQ = 300.0
W_RECOVERY_OUT = 200.0
W_RECOVERY_GAP = 100.0
W_FLAGGED = 0.0


def _within_band(fraction: float) -> float:
    """A 0..1 magnitude as a 0..99 tie-break, so it can never reach the next band up."""
    return min(max(fraction, 0.0), 1.0) * 99.0


#: Why a file has no usable internal standard, as the analyst reads it. Keys are the reasons in
#: ``core.pipeline.MISSING_IS_FLAG``.
MISSING_IS_REASONS: dict[str, str] = {
    "absent": "not in this injection",
    "not_found": "not found at its expected retention time",
    "rejected": "integrated, but inconsistent with the rest of the run",
}


def missing_standards(result) -> dict[str, dict[str, str]]:
    """File id -> {internal-standard label: reason} from the run's sample-level flags."""
    out: dict[str, dict[str, str]] = {}
    prefix = "missing_internal_standard:"
    for file_id, flags in (getattr(result, "flags", None) or {}).items():
        for flag in flags or ():
            if flag.startswith(prefix):
                label, _, why = flag[len(prefix):].partition("=")
                out.setdefault(str(file_id), {})[label] = why or "absent"
    return out


def missing_standard_text(missing: dict[str, str], whole_file: bool) -> str:
    """One sentence for a file's missing standards, ending in what it costs."""
    parts = [f"{label} {MISSING_IS_REASONS.get(why, why)}" for label, why in missing.items()]
    cost = "no mass reported for this file" if whole_file else (
        "no mass reported for the congeners " + " and ".join(missing) + " normalises")
    return "; ".join(parts) + ": " + cost


@dataclass(frozen=True)
class FailingItem:
    """One row of the Overview's "needs attention" table."""

    where: str
    what: str
    value: str
    tone: str = "warning"
    #: Sorts worst-first. Larger is worse. See the severity bands above.
    weight: float = 0.0


@dataclass
class RunSummary:
    """Every number the Overview tab shows, and the reasons behind the verdict."""

    n_samples: int = 0
    n_blanks: int = 0
    recovery_low: float | None = None
    recovery_high: float | None = None
    n_recoveries: int = 0
    n_out_of_band: int = 0
    flagged_integrations: int = 0
    total_peaks: int = 0
    censored_cells: int = 0
    censorable_cells: int = 0
    unusable: list[str] = field(default_factory=list)
    #: File id -> {internal-standard label: reason}, for every reported file lacking one.
    missing_standards: dict[str, dict[str, str]] = field(default_factory=dict)
    failing: list[FailingItem] = field(default_factory=list)
    fit_ok: bool = True
    fit_reason: str = ""

    @property
    def flagged_fraction(self) -> float:
        return (self.flagged_integrations / self.total_peaks) if self.total_peaks else 0.0

    @property
    def recovery_range_text(self) -> str:
        if self.recovery_low is None or self.recovery_high is None:
            return "n/a"
        return f"{self.recovery_low * 100:.0f}–{self.recovery_high * 100:.0f}%"


def iter_recoveries(frame) -> list[tuple[str, str, float]]:
    """``(homolog, file id, fraction)`` for every measurable cell of a recovery frame."""
    out: list[tuple[str, str, float]] = []
    if frame is None or getattr(frame, "empty", True):
        return out
    columns = [str(c) for c in frame.columns]
    for homolog in frame.index:
        row = frame.loc[homolog]
        for file_id in columns:
            value = row[file_id]
            if _is_number(value):
                out.append((str(homolog), file_id, float(value)))
    return out


@dataclass(frozen=True)
class Censoring:
    """What the limit of quantification did to this run, run-wide and per file."""

    #: Native cells set to 0 because they were at or below the limit.
    censored: int = 0
    #: Native cells that carried a mass and could therefore have been censored.
    censorable: int = 0
    #: File id -> native congeners still reporting a mass after censoring.
    kept: dict[str, int] = field(default_factory=dict)
    #: File id -> native congeners that carried a mass before censoring.
    measurable: dict[str, int] = field(default_factory=dict)


def censoring(result, reported: list[str]) -> Censoring:
    """The censoring decision over native analytes x reported files, run-wide and per file.

    Vectorised through the frames rather than iterated: the scope is 185 analytes by
    twenty-odd files, and a per-cell ``.loc`` over 4000 cells is a visible pause on a page
    that has to feel instant. The per-file counts come out of the same two frames in the same
    pass, because the Overview needs both the headline number and the one file that reported
    nothing.
    """
    analytes = getattr(result, "analytes", None) or []
    labels = [a.label for a in analytes if getattr(a, "role", "") == "native"]
    corrected = getattr(result, "masses_corrected", None)
    censored_frame = getattr(result, "masses_censored", None)
    if corrected is None or censored_frame is None or getattr(corrected, "empty", True):
        return Censoring()
    rows = [label for label in labels if label in corrected.index]
    columns = [fid for fid in reported if fid in corrected.columns]
    if not rows or not columns:
        return Censoring()
    before = corrected.loc[rows, columns]
    after = censored_frame.loc[rows, columns]
    zeroed = (after == 0) & (before > 0)
    return Censoring(
        censored=int(zeroed.to_numpy().sum()),
        censorable=int(before.notna().to_numpy().sum()),
        kept={fid: int((after[fid] > 0).sum()) for fid in columns},
        measurable={fid: int(before[fid].notna().sum()) for fid in columns},
    )


def recovery_gaps(frame, files: list[str]) -> dict[str, list[str]]:
    """``{file id: homologs with no recovery}``, for homologs this run measured elsewhere.

    Scoped to homologs the run measured *somewhere*. A run carries a labelled surrogate for
    each homolog it was spiked for and for no others, so an all-empty row is a statement
    about the method rather than about any file. Listing ten of those per file would bury the
    row that means something: a homolog measured in twenty-one files and missing in the
    twenty-second, whose masses in that file are therefore reported uncorrected while its
    neighbours' are not.
    """
    if frame is None or getattr(frame, "empty", True):
        return {}
    columns = [str(c) for c in frame.columns if str(c) in set(files)]
    gaps: dict[str, list[str]] = {}
    for homolog in frame.index:
        row = frame.loc[homolog]
        measured = [fid for fid in columns if _is_number(row[fid])]
        if not measured or len(measured) == len(columns):
            continue  # measured nowhere (not this run's homolog) or everywhere (no gap)
        for file_id in columns:
            if not _is_number(row[file_id]):
                gaps.setdefault(file_id, []).append(str(homolog))
    return gaps


def per_injection_models(result) -> dict[str, dict[str, float]]:
    """Each calibration injection's own RRF table, recomputed from its kept peaks.

    The pipeline averages the injections and keeps only the *names* of the ones that
    contributed (``CalibrationModel.injections``), which is the right thing for a result
    object. But the spread the average hid is what both the Calibration tab and the
    Overview's response-factor row are about. Since every injection's peaks are on the
    result, re-deriving each model costs a dictionary comprehension over 185 peaks and needs
    no change to the pipeline. It is the *same* function the pipeline called, so the numbers
    cannot disagree.

    ``core.quantify`` (and therefore numpy) is imported here rather than at module scope, and
    that is not merely style: see this module's import discipline. It costs nothing in
    practice, because this is only ever reached from a screen that has a ``RunResult`` on it,
    and a ``RunResult`` can only exist if ``core.pipeline`` already imported both.
    """
    from quant209.core.quantify import build_calibration

    analytes = getattr(result, "analytes", None) or []
    files = getattr(result, "files", None) or {}
    models: dict[str, dict[str, float]] = {}
    for file_id in getattr(result, "cal_ids", None) or []:
        record = files.get(file_id)
        peaks = list(getattr(record, "peaks", None) or [])
        if not peaks:
            continue
        models[str(file_id)] = dict(build_calibration(peaks, analytes, name=str(file_id)).rrf)
    return models


def rrf_spread(models: dict[str, dict[str, float]], limit: float = RRF_SPREAD_LIMIT) -> list[str]:
    """Analytes whose per-injection response factors disagree by more than ``limit`` (RSD).

    Relative, not absolute: response factors span orders of magnitude across the congener
    table, so a fixed spread would flag every large one and no small one. Analytes that only
    one injection produced a factor for are not judged. One measurement has no spread, and
    calling that agreement would be the wrong answer in the reassuring direction.

    ``statistics`` rather than numpy: this is a mean and a standard deviation over at most a
    handful of floats per analyte, and it keeps the arithmetic on the Overview's own path
    free of an array library.
    """
    from statistics import fmean, stdev

    names = list(models)
    if len(names) < 2:
        return []
    labels: list[str] = list(dict.fromkeys(label for model in models.values() for label in model))
    out: list[str] = []
    for label in labels:
        values = [models[name][label] for name in names if label in models[name]]
        if len(values) < 2:
            continue
        mean = fmean(values)
        if mean <= 0:
            continue
        if stdev(values) / mean > limit:
            out.append(label)
    return out


def summarise(result) -> RunSummary:
    """Reduce a ``RunResult`` to the Overview tab's five numbers and its failing items.

    The failing items are the triage list, and it is not one class of thing. A table whose
    ``What`` column reads ``flagged integrations`` for every row carries no information: it
    is really just a sorted count, and cannot tell an analyst what *kind* of trouble the run
    is in. Every class below is a different answer to "what is wrong", each is written so the
    ``What`` column names it, and each contributes no rows at all on a run that does not have
    it:

    ==========================  ===========================================================
    retention-time fit          the run-fit guard refused this run's calibration injection
    decoded to no chromatograms a file produced nothing integrable
    reference injection         no usable reference, so no recovery correction anywhere
    blanks                      none, so no limit of quantification and nothing censored
    response factors            congeners the calibration guard refused outright
    response-factor spread      congeners whose injections disagree (:data:`RRF_SPREAD_LIMIT`)
    no congener above the limit a *sample* in which nothing cleared the LOQ
    <homolog> surrogate recovery  outside the acceptance band, per file and homolog
    <homolog> recovery gap      measured elsewhere in the run, missing here
    flagged integrations        the advisory count, per file
    ==========================  ===========================================================

    Two scoping decisions matter here, because each is the difference between a useful row
    and a wall of them:

    * **"Nothing above the limit" is samples only.** A blank whose every congener censors to
      zero is a *clean blank*, the outcome the limit exists to produce. Reporting the run's
      seven blanks as seven problems would be the screen misreading its own QA/QC.
    * **It is "nothing", not "heavily".** Across one run's fifteen samples the censored share
      runs from 59% to 96%, so any fractional "heavily censored" threshold is a knob fitted to
      one run's chemistry. Whether a file reported *anything at all* is a property of the file,
      needs no threshold, and on that run picks out exactly the one sample that reported
      nothing.
    """
    summary = RunSummary()
    fit = getattr(result, "fit", None)
    summary.fit_ok = bool(getattr(fit, "fit_ok", True))
    summary.fit_reason = str(getattr(fit, "reason", "") or "")

    sample_ids = [str(s) for s in (getattr(result, "sample_ids", None) or [])]
    blank_ids = [str(s) for s in (getattr(result, "blank_ids", None) or [])]
    files = getattr(result, "files", None) or {}
    summary.n_samples, summary.n_blanks = len(sample_ids), len(blank_ids)
    reported = sample_ids + blank_ids

    if not summary.fit_ok:
        summary.failing.append(
            FailingItem(
                "this run",
                "retention-time fit",
                summary.fit_reason
                or "the retention-time table does not describe this run's calibration",
                tone="error",
                weight=W_FIT,
            )
        )

    # ---- recoveries: out of band, then measured elsewhere but missing here -------------
    for homolog, file_id, value in iter_recoveries(getattr(result, "recoveries", None)):
        summary.n_recoveries += 1
        summary.recovery_low = value if summary.recovery_low is None else min(summary.recovery_low, value)
        summary.recovery_high = value if summary.recovery_high is None else max(summary.recovery_high, value)
        if RECOVERY_LOW <= value <= RECOVERY_HIGH:
            continue
        summary.n_out_of_band += 1
        summary.failing.append(
            FailingItem(
                file_id,
                f"{homolog} surrogate recovery",
                f"{value * 100:.0f}%",
                tone="error" if value < 0.25 or value > 2.0 else "warning",
                weight=W_RECOVERY_OUT + _within_band(abs(value - 1.0)),
            )
        )

    for file_id, homologs in recovery_gaps(getattr(result, "recoveries", None), reported).items():
        for homolog in homologs:
            summary.failing.append(
                FailingItem(
                    file_id,
                    f"{homolog} recovery gap",
                    f"no {HOMOLOG_SURROGATE.get(homolog, homolog)} recovery, so "
                    f"{homolog} masses are reported uncorrected here",
                    weight=W_RECOVERY_GAP + _within_band(len(homologs) / len(HOMOLOG_SURROGATE)),
                )
            )

    # ---- integrations ------------------------------------------------------------------
    for file_id in reported:
        record = files.get(file_id)
        if record is None:
            continue
        peaks = list(getattr(record, "peaks", None) or [])
        if not peaks:
            summary.unusable.append(file_id)
            summary.failing.append(
                FailingItem(
                    file_id,
                    "decoded to no chromatograms",
                    "not quantified",
                    tone="error",
                    weight=W_NO_PEAKS,
                )
            )
            continue
        summary.total_peaks += len(peaks)
        flagged = sum(
            1 for peak in peaks if any(f in QUALITY_FLAGS for f in (peak.flags or ()))
        )
        summary.flagged_integrations += flagged
        if flagged / len(peaks) > FLAGGED_FRACTION_LIMIT:
            summary.failing.append(
                FailingItem(
                    file_id,
                    "flagged integrations",
                    f"{flagged} of {len(peaks)} peaks",
                    weight=W_FLAGGED + _within_band(flagged / len(peaks)),
                )
            )

    # ---- internal standards: no standard, no mass ----------------------------------------
    is_labels = [str(getattr(a, "label", "")) for a in (getattr(result, "analytes", None) or [])
                 if str(getattr(a, "role", "")).startswith("IS_")]
    for file_id, missing in missing_standards(result).items():
        if file_id not in reported:
            continue
        summary.missing_standards[file_id] = missing
        whole = bool(is_labels) and set(missing) >= set(is_labels)
        summary.failing.append(
            FailingItem(
                file_id,
                "internal standard missing",
                missing_standard_text(missing, whole),
                tone="error",
                weight=W_NO_PEAKS,
            )
        )

    # ---- the run's own losses -----------------------------------------------------------
    if not blank_ids:
        summary.failing.append(
            FailingItem(
                "this run",
                "blanks",
                "none, so there is no limit of quantification and nothing is censored",
                weight=W_NO_BLANKS,
            )
        )

    summary.failing.extend(_reference_items(result, reported))
    summary.failing.extend(_calibration_items(result))

    # ---- censoring: the headline number, and any sample that reported nothing -----------
    counts = censoring(result, reported)
    summary.censored_cells, summary.censorable_cells = counts.censored, counts.censorable
    for file_id in sample_ids:
        measurable = counts.measurable.get(file_id, 0)
        if measurable and counts.kept.get(file_id, 0) == 0:
            summary.failing.append(
                FailingItem(
                    file_id,
                    "no congener above the limit",
                    f"0 of {measurable} congeners cleared the limit of quantification",
                    weight=W_NOTHING_ABOVE_LOQ,
                )
            )

    summary.failing.sort(key=lambda item: item.weight, reverse=True)
    return summary


def _reference_items(result, reported: list[str]) -> list[FailingItem]:
    """The recovery denominator's own failures: no reference at all, or an unpaired file.

    ``reference_for`` absent entirely (rather than empty) means the result predates the
    pairing record and nothing can be claimed about it, so nothing is.
    """
    pairing = getattr(result, "reference_for", None)
    if pairing is None or not reported:
        return []
    if not pairing:
        return [
            FailingItem(
                "this run",
                "reference injection",
                "none usable, so every mass is reported uncorrected for recovery",
                tone="error",
                weight=W_NO_REFERENCE,
            )
        ]
    return [
        FailingItem(
            file_id,
            "reference pairing",
            "no bracketing reference injection, so masses are reported uncorrected",
            tone="error",
            weight=W_NO_REFERENCE,
        )
        for file_id in reported
        if file_id not in pairing
    ]


def _calibration_items(result) -> list[FailingItem]:
    """The two calibration classes: factors refused outright, and factors that disagree.

    Both are stated as one run-level row carrying a count, never one row per analyte. They
    are properties of the *run's calibration* rather than of any file, and on a 185-analyte
    table the per-analyte form would be twenty-five identical ``What`` cells that carry no
    information. The Calibration tab is where the analyst goes next. It already names every
    one of them, per injection, with the guard's own reason.
    """
    items: list[FailingItem] = []
    for name, reason in (getattr(result, "excluded_calibrations", None) or {}).items():
        items.append(
            FailingItem(
                str(name),
                "calibration injection excluded",
                f"{reason}, so it anchors no file and supplies no response factor",
                weight=W_RRF_REFUSED,
            )
        )
    calibration = getattr(result, "calibration", None)
    if calibration is None:
        return items
    invalid = dict(getattr(calibration, "invalid", None) or {})
    analytes = getattr(result, "analytes", None) or []
    if invalid:
        items.append(
            FailingItem(
                "this run",
                "response factors",
                f"{len(invalid)} of {len(analytes)} congeners have none, so they are "
                "not quantified in any file",
                weight=W_RRF_REFUSED,
            )
        )
    models = per_injection_models(result)
    unstable = rrf_spread(models)
    if unstable:
        items.append(
            FailingItem(
                "this run",
                "response-factor spread",
                f"{len(unstable)} congeners vary by more than {RRF_SPREAD_LIMIT:.0%} across "
                f"{len(models)} calibration injections",
                weight=W_RRF_SPREAD,
            )
        )
    return items


def verdict(summary: RunSummary) -> tuple[str, str, str]:
    """``(tone, headline, why)`` for the Overview banner.

    Three states and one rule each, so the colour is never a matter of taste:

    * **red**: the run-fit guard failed. Every number was still computed, and the reason is
      carried *verbatim* from ``core.calibrate.rt_fit_quality`` rather than paraphrased,
      because paraphrasing a scientific verdict is how one becomes wrong. Also red: a reported
      file with no usable internal standard, whose masses are therefore withheld.
    * **amber**: any surrogate recovery outside the band, or more than
      :data:`FLAGGED_FRACTION_LIMIT` of a file's integrations flagged, or no blanks, or a
      file that decoded to nothing.
    * **green**: none of the above. The sub-line still says what was checked. Stating "looks
      good" without the evidence would be a claim the analyst has no way to audit.
    """
    if not summary.fit_ok:
        return (
            "error",
            "Run failed the retention-time fit check",
            summary.fit_reason
            or "retention-time table does not fit this run's calibration injection",
        )
    if summary.missing_standards:
        files = sorted(summary.missing_standards)
        return (
            "error",
            "Internal standard missing",
            f"{len(files)} file(s) without a usable internal standard, so their masses are "
            "withheld: " + ", ".join(files),
        )

    reasons: list[str] = []
    if summary.n_out_of_band:
        reasons.append(
            f"{summary.n_out_of_band} of {summary.n_recoveries} recoveries outside "
            f"{RECOVERY_LOW:.0%}–{RECOVERY_HIGH:.0%}"
        )
    if summary.flagged_fraction > FLAGGED_FRACTION_LIMIT:
        reasons.append(
            f"{summary.flagged_integrations} of {summary.total_peaks} integrations flagged"
        )
    if not summary.n_blanks:
        reasons.append("no blanks in this run")
    if summary.unusable:
        reasons.append(
            f"{len(summary.unusable)} file(s) decoded to no chromatograms"
        )
    if reasons:
        return ("warning", "Review required", "   ·   ".join(reasons))

    parts = ["retention-time table fits"]
    if summary.n_recoveries:
        parts.append(
            f"{summary.n_recoveries} recoveries within "
            f"{RECOVERY_LOW:.0%}–{RECOVERY_HIGH:.0%}"
        )
    parts.append(
        f"{summary.flagged_integrations} of {summary.total_peaks} integrations flagged"
    )
    return ("success", "No problems found", "   ·   ".join(parts))


# --------------------------------------------------------------------------------------
# Small shared widgets
# --------------------------------------------------------------------------------------


class StatCard(QFrame):
    """One number and its caption. Five of them are the Overview's second glance."""

    def __init__(self, caption: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("StatCard")
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        layout = QVBoxLayout(self)
        pad = theme.SPACING + 6
        layout.setContentsMargins(pad, pad, pad, pad)
        layout.setSpacing(2)
        self.value_label = QLabel("—")
        self.value_label.setObjectName("StatValue")
        layout.addWidget(self.value_label)
        self.caption_label = QLabel(caption)
        self.caption_label.setObjectName("StatCaption")
        self.caption_label.setWordWrap(True)
        layout.addWidget(self.caption_label)
        self._tone = ""

    def set_value(self, text: str, tone: str = "") -> None:
        self.value_label.setText(text)
        if tone != self._tone:
            self._tone = tone
            self.value_label.setProperty("tone", tone or None)
            theme.repolish(self.value_label)


def review_table() -> QTableView:
    """A read-only table, configured the one way every table on this page is configured."""
    table = QTableView()
    table.setAlternatingRowColors(False)
    table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
    table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    table.setShowGrid(False)
    table.setWordWrap(False)
    table.setCornerButtonEnabled(False)
    table.verticalHeader().setVisible(False)
    table.verticalHeader().setDefaultSectionSize(ROW_HEIGHT)
    header = table.horizontalHeader()
    header.setHighlightSections(False)
    header.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
    return table


class ReviewTableModel(QAbstractTableModel):
    """A read-only table over a list of row objects, shared by all five tables here.

    Subclasses declare :attr:`HEADERS` and implement :meth:`display`. Everything else
    (counts, headers, tinting, tooltips, the flag-chip role) has one implementation, so
    five tables on one screen cannot render the same idea five slightly different ways.
    """

    HEADERS: tuple[str, ...] = ()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._rows: list[Any] = []

    # ---- population ---------------------------------------------------------------------

    def set_rows(self, rows: list[Any]) -> None:
        self.beginResetModel()
        self._rows = list(rows)
        self.endResetModel()

    def row_at(self, row: int) -> Any:
        return self._rows[row]

    @property
    def rows(self) -> list[Any]:
        return self._rows

    # ---- subclass hooks -------------------------------------------------------------------

    def display(self, row: Any, column: int) -> str:  # pragma: no cover: overridden
        raise NotImplementedError

    def foreground(self, row: Any, column: int) -> QColor | None:
        return None

    def background(self, row: Any, column: int) -> QColor | None:
        return None

    def flags_for(self, row: Any, column: int) -> list[str] | None:
        return None

    def tooltip(self, row: Any, column: int) -> str | None:
        return None

    def font(self, row: Any, column: int) -> QFont | None:
        """Per-cell type. Used only to weight the one column a table is *about*."""
        return None

    def alignment(self, column: int) -> Qt.AlignmentFlag:
        return Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter

    # ---- QAbstractTableModel ---------------------------------------------------------------

    def rowCount(self, parent=QModelIndex()) -> int:  # noqa: N802, B008 (Qt's spelling)
        return 0 if parent.isValid() else len(self._rows)

    def columnCount(self, parent=QModelIndex()) -> int:  # noqa: N802, B008
        return 0 if parent.isValid() else len(self.HEADERS)

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):  # noqa: N802
        if orientation != Qt.Orientation.Horizontal:
            return None
        if role == Qt.ItemDataRole.DisplayRole:
            return self.HEADERS[section]
        if role == Qt.ItemDataRole.TextAlignmentRole:
            # A header follows its column, always. "Area (counts·s)" left-aligned over
            # right-aligned figures leaves the caption floating three centimetres from the
            # numbers it names.
            return int(self.alignment(section))
        return None

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        row = self._rows[index.row()]
        column = index.column()
        if role == Qt.ItemDataRole.DisplayRole:
            return self.display(row, column)
        if role == Qt.ItemDataRole.ForegroundRole:
            return self.foreground(row, column)
        if role == Qt.ItemDataRole.BackgroundRole:
            return self.background(row, column)
        if role == Qt.ItemDataRole.ToolTipRole:
            return self.tooltip(row, column)
        if role == Qt.ItemDataRole.FontRole:
            return self.font(row, column)
        if role == Qt.ItemDataRole.TextAlignmentRole:
            return int(self.alignment(column))
        if role == FLAGS_ROLE:
            return self.flags_for(row, column)
        return None


class FlagChipDelegate(QStyledItemDelegate):
    """Paints ``index.data(FLAGS_ROLE)`` as a row of small pills.

    Chips rather than comma-separated words because a flags column is scanned down, not
    read across: the eye finds three amber pills in a column of grey ones instantly and
    finds ``low_sn, narrow`` in a wall of text not at all.
    """

    PADDING = 10
    GAP = 5
    PILL_PAD = 8

    def _chip_font(self, base: QFont) -> QFont:
        font = QFont(base)
        font.setPointSizeF(theme.FONT_SIZE_SMALL_PT)
        font.setBold(True)
        return font

    def paint(self, painter: QPainter, option, index) -> None:
        opt = QStyleOptionViewItem(option)
        self.initStyleOption(opt, index)
        opt.text = ""
        widget = opt.widget
        style = widget.style() if widget is not None else QApplication.style()
        style.drawControl(QStyle.ControlElement.CE_ItemViewItem, opt, painter, widget)

        flags = index.data(FLAGS_ROLE) or []
        if not flags:
            return
        font = self._chip_font(opt.font)
        metrics = QFontMetrics(font)
        height = metrics.height() + 2
        x = opt.rect.left() + self.PADDING
        top = opt.rect.center().y() - height / 2 + 1
        painter.save()
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setFont(font)
        for flag in flags:
            text = flag_text(flag)
            width = metrics.horizontalAdvance(text) + 2 * self.PILL_PAD
            if x + width > opt.rect.right() - 2:
                break
            tone = FLAG_TONES.get(flag, "neutral")
            pill = QRectF(x, top, width, height)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QColor(_tone_fill(tone)))
            painter.drawRoundedRect(pill, height / 2, height / 2)
            painter.setPen(QColor(_tone_text(tone)))
            painter.drawText(pill, int(Qt.AlignmentFlag.AlignCenter), text)
            x += width + self.GAP
        painter.restore()

    def sizeHint(self, option, index) -> QSize:  # noqa: N802
        flags = index.data(FLAGS_ROLE) or []
        metrics = QFontMetrics(self._chip_font(option.font))
        width = 2 * self.PADDING + sum(
            metrics.horizontalAdvance(flag_text(f)) + 2 * self.PILL_PAD + self.GAP for f in flags
        )
        return QSize(int(width), max(super().sizeHint(option, index).height(), ROW_HEIGHT))


class QualityMarkDelegate(QStyledItemDelegate):
    """Draws a short amber rule down the left edge of a row whose integration is flagged.

    Why a rule and not the row tint the model already supplies
    -----------------------------------------------------------
    :meth:`PeakTableModel.background` returns ``warning_soft`` for a flagged row and **it
    never reaches the screen**: the application stylesheet has a ``QTableView::item`` rule
    (for the row padding and the hairline separator), and once a stylesheet claims an item's
    box, ``QStyleSheetStyle`` paints the panel itself and ignores ``Qt::BackgroundRole``. The
    model role stays, because it answers "is this row flagged?" for the tests, but the
    *visible* mark is painted here, where nothing can swallow it.

    A rule is also the better mark. ``warning_soft`` on white is a 1.03:1 wash that vanishes
    on a lab monitor at any brightness. A 3-px bar at the start of the row is found by an eye
    scanning the analyte column, survives greyscale, and does not tint the numbers.
    """

    BAR_WIDTH = 3
    INSET = 4

    def paint(self, painter: QPainter, option, index) -> None:
        super().paint(painter, option, index)
        colour = index.data(Qt.ItemDataRole.BackgroundRole)
        if colour is None:
            return
        rect = option.rect
        painter.save()
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor(theme.PALETTE["warning"]))
        bar = QRectF(
            rect.left() + self.INSET,
            rect.top() + self.INSET + 1,
            self.BAR_WIDTH,
            max(rect.height() - 2 * (self.INSET + 1), 1),
        )
        painter.drawRoundedRect(bar, self.BAR_WIDTH / 2, self.BAR_WIDTH / 2)
        painter.restore()


def preview_placeholder(text: str, expanding: bool = True) -> QLabel:
    """The centred line a preview panel shows *instead of* an empty chart.

    Building the :class:`~quant209.app.charts.ChromatogramView` immediately and leaving it on
    screen with no data in it renders as a framed, gridded axis box running 0.1 to 0.9,
    indistinguishable from a chart that failed to load, and on two of the five tabs that would
    be the first thing an analyst sees. The plot stays hidden until it has a trace, and this
    stands in its place.

    ``expanding`` decides whether the emptiness fills the panel or hugs its own line, and the
    two callers want different answers. The Calibration tab's preview is a fixed-width
    *column* beside a full-height table, so an empty one that did not fill would leave a short
    card floating at the top of a tall gap. The Recoveries tab's preview is a full-width panel
    *below* a ten-row grid, where filling would turn one sentence into a 45%-of-screen sheet
    of white, reading as a screen that had not finished loading. See :class:`RecoveriesTab`.
    """
    label = hint(text)
    label.setWordWrap(True)
    label.setAlignment(Qt.AlignmentFlag.AlignCenter)
    label.setSizePolicy(
        QSizePolicy.Policy.Preferred,
        QSizePolicy.Policy.Expanding if expanding else QSizePolicy.Policy.Fixed,
    )
    return label


def filter_chip(text: str) -> QPushButton:
    """A checkable pill for the QA/QC filter row (``QPushButton[variant="filter"]``).

    Built here rather than in ``workflow.base``: it is the only checkable button in the
    application, and a shared vocabulary earns its place on the second page that wants it,
    not the first.
    """
    button = QPushButton(text)
    button.setProperty("variant", "filter")
    button.setCheckable(True)
    button.setCursor(Qt.CursorShape.PointingHandCursor)
    return button


def _tone_fill(tone: str) -> str:
    return theme.PALETTE.get(f"{tone}_soft", theme.PALETTE["surface_sunken"])


def _tone_text(tone: str) -> str:
    return theme.PALETTE.get(tone, theme.PALETTE["text_secondary"])


# --------------------------------------------------------------------------------------
# The chromatogram cache
# --------------------------------------------------------------------------------------


class ChromatogramCache:
    """The last :data:`CACHE_FILES` files' chromatograms, keyed by mzML path.

    An ``OrderedDict`` used as an LRU rather than ``functools.lru_cache``: the entries are
    tens of megabytes each and the page must be able to drop all of them when the run
    changes, which a decorator cache on a module-level function cannot be asked to do
    per-instance.
    """

    def __init__(self, size: int = CACHE_FILES) -> None:
        self.size = size
        self._items: OrderedDict[str, dict] = OrderedDict()

    def __len__(self) -> int:
        return len(self._items)

    def clear(self) -> None:
        self._items.clear()

    def cached(self, record) -> bool:
        path = record_path(record)
        return path is not None and str(path) in self._items

    def load(self, record, on_cold: Callable[[Path], None] | None = None) -> dict:
        """This file's ``{(q1, q3): Chromatogram}``, from memory when it is there.

        ``on_cold`` is called *before* the read with the path about to be opened, so the
        caller can say so on screen. A cold read blocks the UI thread for about a second.
        Blocking silently is what makes an application feel broken rather than busy.
        """
        path = record_path(record)
        if path is None:
            return {}
        key = str(path)
        hit = self._items.get(key)
        if hit is not None:
            self._items.move_to_end(key)
            return hit
        if on_cold is not None:
            on_cold(path)
        from quant209.core.mzml import read_srm_mzml

        chroms = read_srm_mzml(path) if path.exists() else {}
        self._items[key] = chroms
        while len(self._items) > self.size:
            self._items.popitem(last=False)
        return chroms


def record_path(record) -> Path | None:
    """The mzML a ``FileRecord`` was read from, under either of its two attribute names."""
    path = getattr(record, "mzml", None) or getattr(record, "mzml_path", None)
    return Path(path) if path is not None else None


# --------------------------------------------------------------------------------------
# Transitions
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Transition:
    """One canonical MRM channel, as the Samples tab's combo box offers it."""

    q1: float
    q3: float
    group: str
    labels: tuple[str, ...]

    @property
    def key(self) -> tuple[float, float]:
        return (self.q1, self.q3)

    @property
    def text(self) -> str:
        # The count is spelled out. A bare trailing "· 37" is a number an analyst has to
        # guess the meaning of, on the one control that decides what the screen is showing.
        count = len(self.labels)
        return (
            f"{self.group}  ·  {self.q1:.1f} → {self.q3:.1f}  ·  "
            f"{count} analyte{'' if count == 1 else 's'}"
        )


def canonical_transitions(analytes) -> list[Transition]:
    """The run's 21 canonical transitions, in acquisition (chlorination) order.

    Named by their *homolog* where they carry natives and by their standard's label where
    they do not, because "penta" and "13C 118" are what an analyst calls those two channels,
    not "325.9 → 255.9", which is what the instrument calls them.
    """
    grouped: dict[tuple[float, float], list] = {}
    for analyte in analytes:
        grouped.setdefault((round(analyte.q1, 1), round(analyte.q3, 1)), []).append(analyte)
    transitions: list[Transition] = []
    for (q1, q3), members in sorted(grouped.items()):
        members = sorted(members, key=lambda a: a.cal_rt)
        homologs = [h for h in dict.fromkeys(a.homolog for a in members) if h]
        group = "/".join(homologs) if homologs else members[0].label
        transitions.append(
            Transition(q1=q1, q3=q3, group=group, labels=tuple(a.label for a in members))
        )
    return transitions


# --------------------------------------------------------------------------------------
# Overview
# --------------------------------------------------------------------------------------


class FailingItemsModel(ReviewTableModel):
    #: ``Where``, not ``File``: the run-level classes (the retention-time fit, the missing
    #: blanks, both calibration rows) are properties of the run rather than of any one file,
    #: and they say so in this column. Heading it ``File`` put ``this run`` under it.
    HEADERS = ("Where", "What", "Value")

    def display(self, row: FailingItem, column: int) -> str:
        return (row.where, row.what, row.value)[column]

    def foreground(self, row: FailingItem, column: int) -> QColor | None:
        if column == 2:
            return QColor(theme.PALETTE.get(row.tone, theme.PALETTE["text"]))
        return None


class OverviewTab(QWidget):
    """The verdict, five numbers, and only what is wrong. No charts: this is triage."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        body = QWidget()
        layout = QVBoxLayout(body)
        pad = theme.SPACING * 2
        layout.setContentsMargins(pad, pad, pad, pad)
        layout.setSpacing(theme.SPACING * 2)

        self.banner = Banner()
        layout.addWidget(self.banner)

        stats = QHBoxLayout()
        stats.setContentsMargins(0, 0, 0, 0)
        stats.setSpacing(theme.SPACING + 4)
        self.cards: dict[str, StatCard] = {}
        for key, caption in (
            ("samples", "SAMPLES"),
            ("blanks", "BLANKS"),
            ("recovery", "SURROGATE RECOVERY RANGE"),
            ("flagged", "FLAGGED INTEGRATIONS"),
            ("censored", "CENSORED CELLS"),
        ):
            card = StatCard(caption)
            self.cards[key] = card
            stats.addWidget(card, 1)
        layout.addLayout(stats)

        card = Card("Needs attention")
        self.model = FailingItemsModel(self)
        self.table = review_table()
        self.table.setModel(self.model)
        header = self.table.horizontalHeader()
        # Only the *last* column stretches. With "What" stretching instead, every value ended
        # up pinned to the right edge of a wide monitor: "a blank · flagged integrations" on
        # the left and "41 of 185 peaks" 1300 px away, with nothing in between. Three snug
        # columns and trailing whitespace is how a table is read.
        for column in range(len(FailingItemsModel.HEADERS) - 1):
            header.setSectionResizeMode(column, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(
            len(FailingItemsModel.HEADERS) - 1, QHeaderView.ResizeMode.Stretch
        )
        self.table.setMinimumHeight(180)
        card.body.addWidget(self.table, 1)
        self.empty = EmptyState(
            "Nothing needs attention",
            "",
        )
        card.body.addWidget(self.empty)
        self.attention_card = card
        layout.addWidget(card, 1)
        layout.addStretch(0)

        outer.addWidget(scroll_wrap(body))
        self.summary = RunSummary()

    def set_result(self, result) -> None:
        self.summary = summarise(result) if result is not None else RunSummary()
        summary = self.summary
        tone, headline, why = verdict(summary)
        self.banner.show_message(headline, why, tone=tone)

        self.cards["samples"].set_value(str(summary.n_samples))
        self.cards["blanks"].set_value(
            str(summary.n_blanks), "" if summary.n_blanks else "warning"
        )
        self.cards["recovery"].set_value(
            summary.recovery_range_text, "warning" if summary.n_out_of_band else "success"
        )
        self.cards["flagged"].set_value(
            str(summary.flagged_integrations),
            "warning" if summary.flagged_fraction > FLAGGED_FRACTION_LIMIT else "",
        )
        self.cards["censored"].set_value(
            f"{summary.censored_cells} of {summary.censorable_cells}"
            if summary.censorable_cells
            else "—"
        )

        self.model.set_rows(summary.failing)
        has_items = bool(summary.failing)
        # The count belongs in the title: a table that scrolls says "22 items" once, where
        # counting the visible rows says "17" and is wrong.
        self.attention_card.title_label.setText(
            f"Needs attention  ·  {len(summary.failing)}" if has_items else "Needs attention"
        )
        self.table.setVisible(has_items)
        self.empty.setVisible(not has_items)


# --------------------------------------------------------------------------------------
# Samples: the chromatogram browser
# --------------------------------------------------------------------------------------


class PeakTableModel(ReviewTableModel):
    """The visible transition's peaks, one row each, in elution order."""

    HEADERS = ("Analyte", "RT (min)", "Area (counts·s)", "S/N", "Method", "Flags")

    def display(self, row, column: int) -> str:
        if column == 0:
            return str(row.label)
        if column == 1:
            return f"{row.rt_apex:.2f}"
        if column == 2:
            return f"{row.area:,.0f}"
        if column == 3:
            # Blank, not "n/a" and never a number: an analyte the presence test rejected has
            # its apex outside the integrated region, so its S/N is computed against an
            # extrapolated baseline and comes out large and meaningless (see
            # ``core.peaks.PeakResult.sn``). Printing it would invite exactly the inference
            # the whole detection contract exists to prevent.
            return "" if not _detected(row) else f"{row.sn:,.0f}"
        if column == 4:
            return METHOD_LABELS.get(str(row.method), str(row.method))
        return ""

    def flags_for(self, row, column: int) -> list[str] | None:
        return list(row.flags or ()) if column == 5 else None

    def background(self, row, column: int) -> QColor | None:
        if any(f in QUALITY_FLAGS for f in (row.flags or ())):
            return QColor(theme.PALETTE["warning_soft"])
        return None

    def alignment(self, column: int) -> Qt.AlignmentFlag:
        if column in (1, 2, 3):
            return Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        return super().alignment(column)

    def tooltip(self, row, column: int) -> str | None:
        if column == 3 and not _detected(row):
            return "no S/N: the presence test rejected this analyte, so its apex may lie "\
                   "outside the integrated region and its S/N is not a measurement"
        return None


def _detected(peak) -> bool:
    detected = getattr(peak, "detected", None)
    if detected is not None:
        return bool(detected)
    return "absent" not in (peak.flags or ())


class FileRailDelegate(QStyledItemDelegate):
    """A file's short id, with its role as a small right-aligned pill."""

    PADDING = 10
    ROLE_TONES = {
        "cal": "accent",
        "ref": "accent",
        "sample": "neutral",
        "blank": "info",
        "hex": "neutral",
    }

    def paint(self, painter: QPainter, option, index) -> None:
        opt = QStyleOptionViewItem(option)
        self.initStyleOption(opt, index)
        role = index.data(FLAGS_ROLE)
        if not role:  # a group header (bold, small, secondary), not a selectable row
            super().paint(painter, option, index)
            return
        text = opt.text
        opt.text = ""
        widget = opt.widget
        style = widget.style() if widget is not None else QApplication.style()
        style.drawControl(QStyle.ControlElement.CE_ItemViewItem, opt, painter, widget)

        badge_font = QFont(opt.font)
        badge_font.setPointSizeF(theme.FONT_SIZE_SMALL_PT)
        badge_font.setBold(True)
        metrics = QFontMetrics(badge_font)
        label = ROLE_LABELS.get(role, role)
        badge_width = metrics.horizontalAdvance(label) + 14
        height = metrics.height() + 2
        rect = opt.rect.adjusted(self.PADDING, 0, -self.PADDING, 0)

        painter.save()
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setFont(opt.font)
        painter.setPen(QColor(theme.PALETTE["text"]))
        painter.drawText(
            rect.adjusted(0, 0, -badge_width - 6, 0),
            int(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter),
            text,
        )
        tone = self.ROLE_TONES.get(role, "neutral")
        pill = QRectF(
            rect.right() - badge_width, rect.center().y() - height / 2 + 1, badge_width, height
        )
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor(_tone_fill(tone)))
        painter.drawRoundedRect(pill, height / 2, height / 2)
        painter.setPen(QColor(_tone_text(tone)))
        painter.setFont(badge_font)
        painter.drawText(pill, int(Qt.AlignmentFlag.AlignCenter), label)
        painter.restore()

    def sizeHint(self, option, index) -> QSize:  # noqa: N802
        size = super().sizeHint(option, index)
        return QSize(size.width(), max(size.height(), 28))


class SamplesTab(QWidget):
    """File rail, transition selector, the trace with every integration on it, peak table."""

    def __init__(self, charts, cache: ChromatogramCache, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.charts = charts
        self.cache = cache
        self._result = None
        self._transitions: list[Transition] = []

        outer = QVBoxLayout(self)
        pad = theme.SPACING * 2
        outer.setContentsMargins(pad, pad, pad, pad)
        outer.setSpacing(theme.SPACING)

        split = QSplitter(Qt.Orientation.Horizontal)
        split.setChildrenCollapsible(False)

        self.file_list = QListWidget()
        self.file_list.setObjectName("FileRail")
        self.file_list.setItemDelegate(FileRailDelegate(self.file_list))
        self.file_list.setMinimumWidth(190)
        self.file_list.setMaximumWidth(280)
        self.file_list.currentItemChanged.connect(lambda *_: self.replot())
        split.addWidget(self.file_list)

        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(theme.SPACING + 4, 0, 0, 0)
        right_layout.setSpacing(theme.SPACING)

        controls = QHBoxLayout()
        controls.setContentsMargins(0, 0, 0, 0)
        controls.setSpacing(theme.SPACING)
        controls.addWidget(meta_label("TRANSITION"))
        self.transition_combo = QComboBox()
        self.transition_combo.setMinimumWidth(320)
        self.transition_combo.currentIndexChanged.connect(lambda *_: self.replot())
        controls.addWidget(self.transition_combo)
        self.busy_hint = hint("")
        self.busy_hint.setVisible(False)
        controls.addWidget(self.busy_hint)
        # An explicit spacer, because the *hidden* busy hint cannot be the one that absorbs
        # the row. A layout ignores hidden widgets entirely, so with the hint invisible
        # (which is nearly always), no item carried a stretch, Qt shared the surplus among
        # the label and the combo instead, and "TRANSITION" ended up 340 px from the control
        # it labels.
        controls.addStretch(1)
        self.reset_button = ghost_button("Reset zoom")
        self.reset_button.clicked.connect(lambda: self.chart.reset_zoom())
        controls.addWidget(self.reset_button)
        right_layout.addLayout(controls)

        vertical = QSplitter(Qt.Orientation.Vertical)
        vertical.setChildrenCollapsible(False)
        self.chart = charts.ChromatogramView()
        self.chart.setMinimumHeight(240)
        vertical.addWidget(self.chart)

        self.model = PeakTableModel(self)
        self.table = review_table()
        self.table.setModel(self.model)
        self.table.setItemDelegateForColumn(0, QualityMarkDelegate(self.table))
        self.table.setItemDelegateForColumn(5, FlagChipDelegate(self.table))
        header = self.table.horizontalHeader()
        for column in range(5):
            header.setSectionResizeMode(column, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(5, QHeaderView.ResizeMode.Stretch)
        # Six rows of a 37-analyte channel, not four: the table is the audit trail beside the
        # picture, and a pane that shows a sixth of it makes the analyst scroll to read what
        # they are looking at.
        self.table.setMinimumHeight(6 * ROW_HEIGHT + 40)
        self.table.clicked.connect(self._on_row_clicked)
        vertical.addWidget(self.table)
        vertical.setStretchFactor(0, 5)
        vertical.setStretchFactor(1, 3)
        # Stretch factors decide how a *resize* is shared. The first layout comes from the
        # children's size hints, and a PlotWidget's hint is large enough to leave the table at
        # its bare minimum. These two numbers are scaled to whatever height the tab gets.
        vertical.setSizes([620, 340])
        right_layout.addWidget(vertical, 1)

        split.addWidget(right)
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        outer.addWidget(split, 1)

    # ---- population -----------------------------------------------------------------------

    def set_result(self, result) -> None:
        self._result = result
        self.file_list.blockSignals(True)
        self.transition_combo.blockSignals(True)
        try:
            self._fill_rail(result)
            self._fill_transitions(result)
        finally:
            self.file_list.blockSignals(False)
            self.transition_combo.blockSignals(False)
        self.replot()

    def _fill_rail(self, result) -> None:
        self.file_list.clear()
        if result is None:
            return
        first_file: QListWidgetItem | None = None
        for attribute, caption in RAIL_GROUPS:
            ids = [str(f) for f in (getattr(result, attribute, None) or [])]
            if not ids:
                continue
            header = QListWidgetItem(caption)
            font = QFont(self.file_list.font())
            font.setPointSizeF(theme.FONT_SIZE_SMALL_PT)
            font.setBold(True)
            header.setFont(font)
            header.setForeground(QColor(theme.PALETTE["text_secondary"]))
            # Enabled but not selectable: a disabled item is painted in the disabled palette
            # (grey on grey) and a *selectable* one would let the arrow keys land on a
            # heading that plots nothing.
            header.setFlags(Qt.ItemFlag.ItemIsEnabled)
            self.file_list.addItem(header)
            for file_id in ids:
                item = QListWidgetItem(file_id)
                record = (getattr(result, "files", None) or {}).get(file_id)
                item.setData(Qt.ItemDataRole.UserRole, file_id)
                item.setData(FLAGS_ROLE, getattr(record, "file_type", "sample"))
                self.file_list.addItem(item)
                if first_file is None:
                    first_file = item
        if first_file is not None:
            self.file_list.setCurrentItem(first_file)

    def _fill_transitions(self, result) -> None:
        self.transition_combo.clear()
        self._transitions = []
        if result is None:
            return
        analytes = getattr(result, "analytes", None)
        if not analytes:
            from quant209.core.config import load_canonical_peaks

            analytes = load_canonical_peaks()
        self._transitions = canonical_transitions(analytes)
        for transition in self._transitions:
            self.transition_combo.addItem(transition.text)
        if self._transitions:
            self.transition_combo.setCurrentIndex(0)

    # ---- plotting --------------------------------------------------------------------------

    def current_file_id(self) -> str | None:
        item = self.file_list.currentItem()
        return item.data(Qt.ItemDataRole.UserRole) if item is not None else None

    def current_transition(self) -> Transition | None:
        index = self.transition_combo.currentIndex()
        if 0 <= index < len(self._transitions):
            return self._transitions[index]
        return None

    def replot(self) -> None:
        """Draw the selected (file, transition). The one operation that must stay warm."""
        file_id, transition = self.current_file_id(), self.current_transition()
        self.model.set_rows([])
        if self._result is None or file_id is None or transition is None:
            self.chart.clear_plot()
            return
        record = (getattr(self._result, "files", None) or {}).get(file_id)
        if record is None:
            self.chart.clear_plot()
            return

        wanted = set(transition.labels)
        peaks = [p for p in (getattr(record, "peaks", None) or []) if p.label in wanted]
        peaks.sort(key=lambda p: p.rt_apex)
        chroms = self.cache.load(record, on_cold=self._announce_read)
        self.busy_hint.setVisible(False)
        chrom = self._chromatogram(chroms, transition, peaks)
        if chrom is None:
            self.chart.clear_plot()
            self.busy_hint.setText(
                f"{file_id} carries no {transition.group} transition ({transition.q1:.1f} → "
                f"{transition.q3:.1f})"
            )
            self.busy_hint.setVisible(True)
            return
        self.chart.show_trace(chrom.times, chrom.intensities, peaks)
        self.model.set_rows(peaks)

    @staticmethod
    def _chromatogram(chroms: dict, transition: Transition, peaks: list):
        """This transition's acquired trace, by the peaks' own recorded key when there is one.

        ``PeakResult.src_q1``/``src_q3`` name the transition actually integrated, which for a
        labelled standard is ~0.1 Da from the canonical table's chemically-derived key, so
        the canonical key alone would find no chromatogram for ten of the twenty-one
        channels. The tolerant lookup is the fallback for a file with no peaks at all here.
        """
        for peak in peaks:
            key = (getattr(peak, "src_q1", -1.0), getattr(peak, "src_q3", -1.0))
            chrom = chroms.get(key)
            if chrom is not None:
                return chrom
        from quant209.core.mzml import get_transition

        return get_transition(chroms, transition.q1, transition.q3)

    def _announce_read(self, path: Path) -> None:
        self.busy_hint.setText(f"Reading {path.name}…")
        self.busy_hint.setVisible(True)
        # repaint(), not update(): the read that follows blocks this thread, so a queued
        # paint would be delivered *after* it and the hint would flash up as the wait ended.
        self.busy_hint.repaint()

    def _on_row_clicked(self, index) -> None:
        if not index.isValid():
            return
        self.chart.zoom_to_peak(self.model.row_at(index.row()))


# --------------------------------------------------------------------------------------
# Calibration
# --------------------------------------------------------------------------------------


@dataclass
class CalibrationRow:
    label: str
    homolog: str
    transition: str
    per_injection: tuple[float | None, ...]
    combined: float | None
    flag: str


class CalibrationTableModel(ReviewTableModel):
    """One row per analyte: the response factor it will quantify with, then the evidence.

    Column order is answer-first: ``Combined RRF`` immediately after the analyte's identity,
    the per-injection values after it, the guard's refusal last. Per-injection-then-combined
    is the order the numbers are *computed* in, and it put the one figure that multiplies
    every mass in this run at column eight of nine: at the product's minimum window size it
    was off the right edge of the table, and the tab's whole subject had to be scrolled to.
    """

    def __init__(self, injections: tuple[str, ...] = (), parent: QWidget | None = None) -> None:
        # super() first, always: PySide refuses attribute assignment on a QObject whose base
        # __init__ has not run, and the error it raises names neither the attribute nor here.
        super().__init__(parent)
        self.injections = tuple(injections)

    @property
    def HEADERS(self) -> tuple[str, ...]:  # noqa: N802 (matches the base class attribute)
        return ("Analyte", "Homolog", "Transition", "Combined RRF", *self.injections, "Flag")

    def display(self, row: CalibrationRow, column: int) -> str:
        n = len(self.injections)
        if column == 0:
            return row.label
        if column == 1:
            return row.homolog or "labeled"
        if column == 2:
            return row.transition
        if column == 3:
            return _rrf(row.combined)
        if column < 4 + n:
            return _rrf(row.per_injection[column - 4])
        return row.flag

    def flags_for(self, row: CalibrationRow, column: int) -> list[str] | None:
        return None

    def foreground(self, row: CalibrationRow, column: int) -> QColor | None:
        if column == len(self.HEADERS) - 1 and row.flag:
            return QColor(theme.PALETTE["warning"])
        if row.combined is None and column == 0:
            return QColor(theme.PALETTE["text_secondary"])
        return None

    def alignment(self, column: int) -> Qt.AlignmentFlag:
        if 3 <= column <= 3 + len(self.injections):
            return Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        return super().alignment(column)

    def font(self, row: CalibrationRow, column: int) -> QFont | None:
        """Semibold on the combined RRF: five columns of near-identical four-digit numbers,
        and only one of them multiplies every mass in the run."""
        if column != 3:
            return None
        font = QFont()
        font.setWeight(QFont.Weight.DemiBold)
        return font

    def tooltip(self, row: CalibrationRow, column: int) -> str | None:
        return row.flag or None


def _rrf(value: float | None) -> str:
    return "—" if value is None else f"{value:.4g}"


def calibration_rows(result, models: dict[str, dict]) -> list[CalibrationRow]:
    """Join the canonical analyte table to the combined and per-injection response factors.

    ``models`` is ``{injection name: {label: rrf}}``. An analyte the guard refused carries
    the guard's own reason. One that survived in some injections and not others says so,
    because "averaged over three of four injections" is a different claim from "averaged
    over four" and only one of them is on the report.
    """
    calibration = getattr(result, "calibration", None)
    combined = dict(getattr(calibration, "rrf", None) or {})
    invalid = dict(getattr(calibration, "invalid", None) or {})
    names = list(models)
    rows: list[CalibrationRow] = []
    for analyte in getattr(result, "analytes", None) or []:
        per = tuple(models[name].get(analyte.label) for name in names)
        value = combined.get(analyte.label)
        if value is None:
            flag = invalid.get(analyte.label, "no response factor")
        else:
            dropped = [name for name, rrf in zip(names, per) if rrf is None]
            flag = f"dropped in {', '.join(dropped)}" if dropped else ""
        rows.append(
            CalibrationRow(
                label=analyte.label,
                homolog=analyte.homolog,
                transition=f"{analyte.q1:.1f} → {analyte.q3:.1f}",
                per_injection=per,
                combined=value,
                flag=flag,
            )
        )
    return rows


class CalibrationTab(QWidget):
    """The response-factor table, and the calibration peak behind whichever row is selected."""

    def __init__(self, charts, cache: ChromatogramCache, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.charts = charts
        self.cache = cache
        self._result = None
        self._all_rows: list[CalibrationRow] = []

        outer = QVBoxLayout(self)
        pad = theme.SPACING * 2
        outer.setContentsMargins(pad, pad, pad, pad)
        outer.setSpacing(theme.SPACING)

        top = QHBoxLayout()
        top.setContentsMargins(0, 0, 0, 0)
        top.setSpacing(theme.SPACING)
        self.search = QLineEdit()
        self.search.setPlaceholderText("Filter by congener, homolog or transition…")
        self.search.setClearButtonEnabled(True)
        # A minimum as well as a maximum: without one the field takes its size hint, the
        # placeholder elides to "Filter by conge…", and the only instruction on the tab is
        # the thing that got cut off.
        self.search.setMinimumWidth(300)
        self.search.setMaximumWidth(360)
        self.search.textChanged.connect(self._apply_filter)
        top.addWidget(self.search)
        self.summary_label = hint("")
        top.addWidget(self.summary_label, 1)
        outer.addLayout(top)

        split = QSplitter(Qt.Orientation.Horizontal)
        split.setChildrenCollapsible(False)
        self.model = CalibrationTableModel(parent=self)
        self.table = review_table()
        self.table.setModel(self.model)
        self.table.clicked.connect(self._on_row_clicked)
        split.addWidget(self.table)

        preview = Card("Calibration peak")
        self.preview_caption = hint("")
        self.preview_caption.setWordWrap(True)
        self.preview_caption.setVisible(False)
        preview.body.addWidget(self.preview_caption)
        self.chart = charts.ChromatogramView(compact=True)
        self.chart.setMinimumHeight(200)
        self.chart.setVisible(False)  # see preview_placeholder
        preview.body.addWidget(self.chart, 1)
        self.placeholder = preview_placeholder(
            "Select a row to see the calibration injection's peak for that analyte."
        )
        preview.body.addWidget(self.placeholder, 1)
        # 360, not 460: the columns that matter on this tab are Combined RRF and Flag, and at
        # 460 they fell off the right edge of the table and had to be scrolled to.
        preview.setMinimumWidth(300)
        preview.setMaximumWidth(360)
        split.addWidget(preview)
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 0)
        outer.addWidget(split, 1)

    def set_result(self, result) -> None:
        self._result = result
        if result is None:
            self.model.set_rows([])
            return
        models = per_injection_models(result)
        self.model.injections = tuple(models)
        self._all_rows = calibration_rows(result, models)
        flagged = sum(1 for row in self._all_rows if row.flag)
        self.summary_label.setText(
            f"RRFs from {len(models)} injection{'' if len(models) == 1 else 's'}  ·  "
            f"{flagged} analyte{'' if flagged == 1 else 's'} flagged"
        )
        self._apply_filter(self.search.text())
        header = self.table.horizontalHeader()
        for column in range(self.model.columnCount()):
            header.setSectionResizeMode(column, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(
            self.model.columnCount() - 1, QHeaderView.ResizeMode.Stretch
        )

    def _apply_filter(self, text: str) -> None:
        needle = (text or "").strip().lower()
        if not needle:
            self.model.set_rows(self._all_rows)
            return
        self.model.set_rows(
            [
                row
                for row in self._all_rows
                if needle in row.label.lower()
                or needle in row.homolog.lower()
                or needle in row.transition.lower()
            ]
        )

    def _on_row_clicked(self, index) -> None:
        if not index.isValid() or self._result is None:
            return
        row = self.model.row_at(index.row())
        cal_ids = getattr(self._result, "cal_ids", None) or []
        files = getattr(self._result, "files", None) or {}
        record = files.get(cal_ids[0]) if cal_ids else None
        if record is None:
            return
        peaks = [p for p in (getattr(record, "peaks", None) or []) if p.label == row.label]
        if not peaks:
            self._show_message(
                f"{row.label} was not placed in {cal_ids[0]}: {row.flag or 'no peak'}"
            )
            return
        peak = peaks[0]
        chroms = self.cache.load(record)
        key = (getattr(peak, "src_q1", -1.0), getattr(peak, "src_q3", -1.0))
        chrom = chroms.get(key)
        if chrom is None:
            from quant209.core.mzml import get_transition

            chrom = get_transition(chroms, peak.q1, peak.q3)
        if chrom is None:
            self._show_message(f"{cal_ids[0]} carries no {row.transition} transition")
            return
        # Every analyte on the channel is drawn, then the view is taken to the selected one:
        # a calibration peak is only judgeable next to the neighbours it was split from.
        neighbours = [
            p
            for p in (getattr(record, "peaks", None) or [])
            if (getattr(p, "src_q1", -1.0), getattr(p, "src_q3", -1.0)) == key
        ] or [peak]
        self._show_chart(
            f"{row.label} in {cal_ids[0]}  ·  RT {peak.rt_apex:.2f} min  ·  "
            f"area {peak.area:,.0f}"
        )
        self.chart.show_trace(chrom.times, chrom.intensities, neighbours)
        self.chart.zoom_to_peak(peak, animate=False)

    # ---- the preview panel's two states ------------------------------------------------

    def _show_chart(self, caption: str) -> None:
        self.preview_caption.setText(caption)
        self.preview_caption.setVisible(True)
        self.placeholder.setVisible(False)
        self.chart.setVisible(True)

    def _show_message(self, message: str) -> None:
        """No peak to draw: say why, where the picture would have been."""
        self.chart.clear_plot()
        self.chart.setVisible(False)
        self.preview_caption.setVisible(False)
        self.preview_caption.setText(message)
        self.placeholder.setText(message)
        self.placeholder.setVisible(True)


# --------------------------------------------------------------------------------------
# Recoveries
# --------------------------------------------------------------------------------------


class RecoveryMatrixModel(QAbstractTableModel):
    """``result.recoveries`` as a homolog x file grid. Values stay fractions. Text is %."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.homologs: list[str] = []
        self.files: list[str] = []
        self._values: dict[tuple[str, str], float] = {}

    def set_frame(self, frame) -> None:
        self.beginResetModel()
        self._values = {}
        if frame is None or getattr(frame, "empty", True):
            self.homologs, self.files = [], []
        else:
            self.files = [str(c) for c in frame.columns]
            # Rows in chlorination order, and only the ones this run actually carries.
            index = [str(i) for i in frame.index]
            self.homologs = [h for h in HOMOLOG_SURROGATE if h in index] or index
            for homolog, file_id, value in iter_recoveries(frame):
                self._values[(homolog, file_id)] = value
        self.endResetModel()

    def value_at(self, row: int, column: int) -> float | None:
        return self._values.get((self.homologs[row], self.files[column]))

    def rowCount(self, parent=QModelIndex()) -> int:  # noqa: N802, B008
        return 0 if parent.isValid() else len(self.homologs)

    def columnCount(self, parent=QModelIndex()) -> int:  # noqa: N802, B008
        return 0 if parent.isValid() else len(self.files)

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):  # noqa: N802
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if orientation == Qt.Orientation.Horizontal:
            return self.files[section]
        return self.homologs[section]

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        value = self.value_at(index.row(), index.column())
        if role == Qt.ItemDataRole.DisplayRole:
            return "—" if value is None else f"{value * 100:.0f}%"
        if role == VALUE_ROLE:
            return value
        if role == OUT_OF_BAND_ROLE:
            return value is not None and not (RECOVERY_LOW <= value <= RECOVERY_HIGH)
        if role == Qt.ItemDataRole.ToolTipRole and value is not None:
            return (
                f"{self.homologs[index.row()]} · {self.files[index.column()]}: "
                f"{value * 100:.1f}% (surrogate {HOMOLOG_SURROGATE.get(self.homologs[index.row()], '')})"
            )
        return None


class RecoveryCellDelegate(QStyledItemDelegate):
    """Paints one heatmap cell: diverging fill, band border, dark text on every cell.

    A ``QStyledItemDelegate`` rather than a pyqtgraph image. The numbers *are* the content
    here: an analyst reads "42%" off the grid and writes it down, and an image plot renders
    text as pixels that blur at any scaling. This keeps the type crisp and the row/column
    headers styled like every other table in the application.
    """

    def __init__(self, color_fn: Callable[[float], QColor], parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.color_fn = color_fn

    def paint(self, painter: QPainter, option, index) -> None:
        value = index.data(VALUE_ROLE)
        rect = option.rect
        painter.save()
        if value is None:
            painter.fillRect(rect, QColor(theme.PALETTE["surface"]))
            painter.setPen(QColor(theme.PALETTE["text_disabled"]))
        else:
            painter.fillRect(rect.adjusted(1, 1, -1, -1), self.color_fn(float(value)))
            if index.data(OUT_OF_BAND_ROLE):
                # A border as well as a darker fill: colour alone is not an accessible
                # encoding, and "outside the acceptance band" is the one fact on this grid
                # that must survive a colourblind reader and a black-and-white printout.
                edge = theme.PALETTE["warning"] if value < RECOVERY_LOW else theme.PALETTE["accent"]
                painter.setPen(QColor(edge))
                painter.drawRect(rect.adjusted(1, 1, -2, -2))
            painter.setPen(QColor(theme.PALETTE["text"]))
        if option.state & QStyle.StateFlag.State_Selected:
            # Selection is drawn as a ring rather than as a fill: filling would replace the
            # one thing the cell is for, its colour, with the selection colour.
            painter.save()
            painter.setPen(QColor(theme.PALETTE["accent_pressed"]))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(rect.adjusted(1, 1, -2, -2))
            painter.restore()
        painter.drawText(rect, int(Qt.AlignmentFlag.AlignCenter), index.data() or "")
        painter.restore()

    def sizeHint(self, option, index) -> QSize:  # noqa: N802
        return QSize(64, ROW_HEIGHT)


class RecoveryLegend(QWidget):
    """The diverging ramp, with the acceptance band bracketed on it."""

    STOPS = (0.0, 0.5, 1.0, 1.5, 2.0)
    BAR_HEIGHT = 14
    GAP = 4

    def __init__(self, color_fn: Callable[[float], QColor], parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.color_fn = color_fn
        # Measured, not a round number. A hard-coded 46 px clipped the descenders off the
        # caption's last line. A legend that cannot be read is decoration.
        metrics = QFontMetrics(self._caption_font())
        self.setFixedHeight(self.BAR_HEIGHT + self.GAP + 2 * metrics.height() + self.GAP)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def _caption_font(self) -> QFont:
        font = QFont(self.font())
        font.setPointSizeF(theme.FONT_SIZE_SMALL_PT)
        return font

    def paintEvent(self, event) -> None:  # noqa: N802 (Qt's spelling)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        width = max(self.width() - 2, 1)
        bar_height = self.BAR_HEIGHT
        for x in range(width):
            fraction = 2.0 * x / width
            painter.setPen(self.color_fn(fraction))
            painter.drawLine(x + 1, 0, x + 1, bar_height)
        painter.setPen(QColor(theme.PALETTE["border_strong"]))
        painter.drawRect(1, 0, width, bar_height)

        band_lo, band_hi = int(width * RECOVERY_LOW / 2.0), int(width * RECOVERY_HIGH / 2.0)
        painter.setPen(QColor(theme.PALETTE["text"]))
        painter.drawLine(band_lo + 1, 0, band_lo + 1, bar_height)
        painter.drawLine(band_hi + 1, 0, band_hi + 1, bar_height)

        font = self._caption_font()
        painter.setFont(font)
        painter.setPen(QColor(theme.PALETTE["text_secondary"]))
        metrics = QFontMetrics(font)
        for stop in self.STOPS:
            text = f"{stop * 100:.0f}%"
            x = int(width * stop / 2.0) + 1
            offset = 0 if stop == 0 else (metrics.horizontalAdvance(text) if stop == 2.0 else metrics.horizontalAdvance(text) // 2)
            painter.drawText(x - offset, bar_height + self.GAP + metrics.ascent(), text)
        painter.drawText(
            1,
            bar_height + self.GAP + metrics.height() + metrics.ascent(),
            f"acceptance band {RECOVERY_LOW:.0%}–{RECOVERY_HIGH:.0%}; cells outside it carry a border",
        )


class HeatmapView(QTableView):
    """A grid that asks the layout for exactly its rows' height and no more.

    The recovery grid is ten rows tall whatever the run, so letting it expand leaves half the
    tab as blank white below the cells. Answering with a real
    :meth:`sizeHint` (rather than calling ``setFixedHeight`` from a ``resizeEvent``, which
    resizes a child *during* a layout pass and leaves the widgets below it painted over the
    grid's last two rows) lets the layout do the arithmetic once, in the right order.
    """

    def _content_height(self) -> int:
        model = self.model()
        rows = model.rowCount() if model is not None else 0
        header = self.horizontalHeader()
        # The scrollbar's height is always reserved. Whether one is needed depends on the
        # viewport width, which is not known while the layout is asking for this hint, and a
        # dozen pixels of ground inside the frame is invisible where a clipped bottom row is
        # not.
        return (
            max(header.height(), header.sizeHint().height())
            + rows * ROW_HEIGHT
            + self.horizontalScrollBar().sizeHint().height()
            + 2 * self.frameWidth()
        )

    def sizeHint(self) -> QSize:  # noqa: N802 (Qt's spelling)
        return QSize(super().sizeHint().width(), self._content_height())

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        return QSize(super().minimumSizeHint().width(), self._content_height())


class RecoveriesTab(QWidget):
    """The homolog x sample heatmap, its legend, and the surrogate peak behind one cell.

    Stacked, not side by side
    --------------------------
    The grid is ten rows tall and as many columns wide as the run has files, which makes it a
    *banner*: it wants every pixel of width and almost none of the height. Beside a preview
    card it got two thirds of the width (ten of twenty-two samples, the rest behind a
    horizontal scrollbar) and left half the tab as blank white below ten rows of cells, with
    the legend stranded at the bottom of that emptiness, five hundred pixels from the colours
    it explains. Stacking puts the whole run on screen at once, the legend under the grid it
    belongs to, and the preview in the space the grid does not want.

    The stack sits in a scroll area because grid, legend and preview have a combined minimum
    near 700 px while the tab is about 580 px tall at the product's minimum window height, so
    *something* has to give. Left to Qt, what gave was the grid: four homologs were cut off
    the bottom, silently, with the vertical scrollbar switched off, the one failure mode a QA
    screen must not have. Scrolling the stack keeps every row reachable and puts the fold
    under the preview, where it belongs.

    The preview is a strip until it has something to show
    ------------------------------------------------------
    Adding the preview card at stretch 1 with an expanding placeholder inside it would mean
    that before any cell was clicked (which is *every* first view of this tab), it would be a
    45%-of-screen empty white card carrying one grey sentence. Half a screen showing nothing
    does not read as "waiting for you". It reads as a screen that failed to finish.

    So the card and the spacer below it trade one stretch factor between them
    (:meth:`_set_preview_expanded`): with no selection the card is at stretch 0 and hugs its
    title and its one line, the spacer takes the surplus, and the emptiness is a strip under
    the legend rather than a sheet. On the first click the chart is shown, the card takes the
    stretch and the spacer gives it up, and the trace gets every pixel the grid does not want.
    The alternative (giving the freed height to the grid) was rejected because the grid is
    a banner by construction (:class:`HeatmapView`): ten rows whatever the run, and stretching
    them would make the cells taller without making the heatmap say more.
    """

    def __init__(self, charts, cache: ChromatogramCache, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.charts = charts
        self.cache = cache
        self._result = None

        frame = QVBoxLayout(self)
        frame.setContentsMargins(0, 0, 0, 0)
        frame.setSpacing(0)
        body = QWidget()
        outer = QVBoxLayout(body)
        pad = theme.SPACING * 2
        outer.setContentsMargins(pad, pad, pad, pad)
        outer.setSpacing(theme.SPACING)
        frame.addWidget(scroll_wrap(body))

        self.model = RecoveryMatrixModel(self)
        self.table = HeatmapView()
        self.table.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        self.table.setModel(self.model)
        self.table.setItemDelegate(RecoveryCellDelegate(charts.recovery_color, self.table))
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setShowGrid(False)
        self.table.setCornerButtonEnabled(False)
        self.table.verticalHeader().setDefaultSectionSize(ROW_HEIGHT)
        self.table.verticalHeader().setHighlightSections(False)
        header = self.table.horizontalHeader()
        header.setObjectName("HeatmapHeader")  # tighter side padding (see theme._QSS)
        header.setHighlightSections(False)
        header.setDefaultAlignment(Qt.AlignmentFlag.AlignCenter)
        # Stretch with a floor, not plain Stretch: at 22 samples an equal division puts
        # "100%" into 38 px and the grid becomes a wall of ellipses. Below the floor the
        # header simply grows past the viewport and the view scrolls. A heatmap whose
        # cells are too small to read is not a heatmap.
        header.setMinimumSectionSize(HEATMAP_MIN_COLUMN)
        header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.table.clicked.connect(self._on_cell_clicked)
        outer.addWidget(self.table, 0)
        self.legend = RecoveryLegend(charts.recovery_color)
        outer.addWidget(self.legend, 0)
        self.empty = EmptyState(
            "No recoveries in this run",
            "No reference injection produced usable masses.",
        )
        outer.addWidget(self.empty)

        preview = Card("Surrogate peak")
        self.preview_caption = hint("")
        self.preview_caption.setWordWrap(True)
        self.preview_caption.setVisible(False)
        preview.body.addWidget(self.preview_caption)
        self.chart = charts.ChromatogramView(compact=True)
        # 160, not the 200 the two narrow previews use: this one sits *below* a ten-row grid
        # and a legend, and at the product's minimum window height the three together have to
        # fit without the card's bottom edge running off the screen.
        self.chart.setMinimumHeight(160)
        self.chart.setVisible(False)  # see preview_placeholder
        preview.body.addWidget(self.chart, 1)
        self.placeholder = preview_placeholder(
            "Select a cell to see the surrogate peak it was measured from.", expanding=False
        )
        preview.body.addWidget(self.placeholder, 0)
        self.preview_card = preview
        outer.addWidget(preview, 0)
        # The spacer the preview trades its stretch with. Held by index rather than by
        # reference because ``QBoxLayout`` re-weights an item by position, not by identity.
        self._preview_row = outer.indexOf(preview)
        outer.addStretch(1)
        self._filler_row = outer.count() - 1
        self._outer = outer
        self._set_preview_expanded(False)

    def set_result(self, result) -> None:
        self._result = result
        self.model.set_frame(getattr(result, "recoveries", None) if result is not None else None)
        has_data = self.model.rowCount() > 0 and self.model.columnCount() > 0
        self.table.setVisible(has_data)
        self.legend.setVisible(has_data)
        self.preview_card.setVisible(has_data)
        self.empty.setVisible(not has_data)
        # The row count decides the grid's height hint, so the layout has to be told to ask
        # again: a model reset alone does not invalidate a parent layout.
        self.table.updateGeometry()

    def _on_cell_clicked(self, index) -> None:
        if not index.isValid() or self._result is None:
            return
        homolog = self.model.homologs[index.row()]
        file_id = self.model.files[index.column()]
        surrogate = HOMOLOG_SURROGATE.get(homolog, "")
        record = (getattr(self._result, "files", None) or {}).get(file_id)
        peaks = [p for p in (getattr(record, "peaks", None) or []) if p.label == surrogate]
        value = self.model.value_at(index.row(), index.column())
        if record is None or not peaks:
            self._show_message(f"{surrogate} was not placed in {file_id}")
            return
        peak = peaks[0]
        chroms = self.cache.load(record)
        chrom = chroms.get((getattr(peak, "src_q1", -1.0), getattr(peak, "src_q3", -1.0)))
        if chrom is None:
            from quant209.core.mzml import get_transition

            chrom = get_transition(chroms, peak.q1, peak.q3)
        if chrom is None:
            self._show_message(f"{file_id} carries no {surrogate} transition")
            return
        self._show_chart(
            f"{surrogate} in {file_id}  ·  {homolog} recovery "
            f"{'—' if value is None else f'{value * 100:.0f}%'}  ·  area {peak.area:,.0f}"
        )
        self.chart.show_trace(chrom.times, chrom.intensities, [peak])
        self.chart.zoom_to_peak(peak, animate=False)

    # ---- the preview panel's two states ------------------------------------------------

    def _set_preview_expanded(self, expanded: bool) -> None:
        """Hand the surplus height to the preview card, or to the spacer under it.

        See the class docstring: this is the whole of the "strip until it has something to
        show" behaviour, and it is one call rather than a resize because the layout should do
        the arithmetic once, in the right order.
        """
        self._outer.setStretch(self._preview_row, 1 if expanded else 0)
        self._outer.setStretch(self._filler_row, 0 if expanded else 1)
        self.preview_card.updateGeometry()

    def _show_chart(self, caption: str) -> None:
        self.preview_caption.setText(caption)
        self.preview_caption.setVisible(True)
        self.placeholder.setVisible(False)
        self.chart.setVisible(True)
        self._set_preview_expanded(True)

    def _show_message(self, message: str) -> None:
        """No peak to draw: say why, in the strip, rather than opening an empty panel."""
        self.chart.clear_plot()
        self.chart.setVisible(False)
        self.preview_caption.setVisible(False)
        self.preview_caption.setText(message)
        self.placeholder.setText(message)
        self.placeholder.setVisible(True)
        self._set_preview_expanded(False)


# --------------------------------------------------------------------------------------
# QA/QC
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class FlagEvent:
    """One flag, on one analyte, in one file."""

    file_id: str
    analyte: str
    kind: str
    context: str


class FlagEventModel(ReviewTableModel):
    HEADERS = ("File", "Analyte", "Flag", "Value")

    def display(self, row: FlagEvent, column: int) -> str:
        return (row.file_id, row.analyte, "", row.context)[column]

    def flags_for(self, row: FlagEvent, column: int) -> list[str] | None:
        return [row.kind] if column == 2 else None


def flag_events(result) -> list[FlagEvent]:
    """Every flag event in the run: recoveries, peak flags, and censoring decisions.

    Everything, not just the samples: a calibration injection's flagged peaks are what the
    Calibration tab's refusals are made of, and an audit that could not see them would be
    asking the reader to take the response factors on trust.
    """
    events: list[FlagEvent] = []
    for homolog, file_id, value in iter_recoveries(getattr(result, "recoveries", None)):
        if not (RECOVERY_LOW <= value <= RECOVERY_HIGH):
            events.append(
                FlagEvent(
                    file_id,
                    HOMOLOG_SURROGATE.get(homolog, homolog),
                    "recovery",
                    f"{homolog} recovery {value * 100:.0f}%",
                )
            )

    for file_id, record in (getattr(result, "files", None) or {}).items():
        for peak in getattr(record, "peaks", None) or []:
            for flag in peak.flags or ():
                events.append(
                    FlagEvent(str(file_id), str(peak.label), flag, _flag_context(flag, peak))
                )

    for file_id, missing in missing_standards(result).items():
        for label, why in missing.items():
            events.append(
                FlagEvent(file_id, label, "missing_is",
                          f"{MISSING_IS_REASONS.get(why, why)}: no mass for what it normalises")
            )

    events.extend(_censored_events(result))
    return events


def _flag_context(flag: str, peak) -> str:
    if flag == "absent":
        return f"no resolved peak at RT {peak.rt_apex:.2f}"
    if flag == "low_sn":
        return f"S/N {peak.sn:,.0f}"
    if flag == "narrow":
        return f"{max(peak.i_end - peak.i_start, 0)} samples between boundaries"
    if flag == "clustered":
        return "area measured above a cluster baseline"
    if flag == "forced":
        return "not found by either detection test, integrated by a rule for this congener"
    if flag == "displaced":
        return "peak drifted outside its retention-time window; assigned to this congener by Aroclor abundance"
    if flag == "realigned":
        return ("searched at a retention time corrected for the drift of this file's "
                "labelled standards")
    if flag == "ghost":
        return ("a more chlorinated congener's in-source fragment (ghost peak) on this channel was "
                "kept out of this integration")
    if flag == "apportioned":
        return "merged with a neighbour of comparable abundance; area is its share by a two-peak fit"
    if flag == "unresolved":
        return "not separated from a neighbour in this run's calibration; area is its share of their peak by a two-peak fit"
    return ""


def _censored_events(result) -> list[FlagEvent]:
    """One event per native cell the LOQ censoring set to zero."""
    analytes = getattr(result, "analytes", None) or []
    labels = [a.label for a in analytes if getattr(a, "role", "") == "native"]
    corrected = getattr(result, "masses_corrected", None)
    censored = getattr(result, "masses_censored", None)
    loq = getattr(result, "loq", None)
    if corrected is None or censored is None or getattr(corrected, "empty", True):
        return []
    reported = [str(f) for f in (getattr(result, "sample_ids", None) or [])] + [
        str(f) for f in (getattr(result, "blank_ids", None) or [])
    ]
    rows = [label for label in labels if label in corrected.index]
    columns = [fid for fid in reported if fid in corrected.columns]
    if not rows or not columns:
        return []
    before = corrected.loc[rows, columns]
    after = censored.loc[rows, columns]
    mask = (after == 0) & (before > 0)
    events: list[FlagEvent] = []
    for label in rows:
        row_mask = mask.loc[label]
        for file_id in columns:
            if not bool(row_mask[file_id]):
                continue
            limit = None
            if loq is not None and not getattr(loq, "empty", True) and label in loq.index:
                candidate = loq[label]
                limit = float(candidate) if _is_number(candidate) else None
            measured = float(before.loc[label, file_id])
            events.append(
                FlagEvent(
                    file_id,
                    label,
                    "censored",
                    f"{measured:.3g} ng at or below the {limit:.3g} ng limit"
                    if limit is not None
                    else f"{measured:.3g} ng at or below the limit",
                )
            )
    return events


class QaqcTab(QWidget):
    """Every flag event in the run, filtered by kind and by file, exportable as CSV."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._events: list[FlagEvent] = []

        outer = QVBoxLayout(self)
        pad = theme.SPACING * 2
        outer.setContentsMargins(pad, pad, pad, pad)
        outer.setSpacing(theme.SPACING)

        filters = QHBoxLayout()
        filters.setContentsMargins(0, 0, 0, 0)
        filters.setSpacing(theme.SPACING // 2 + 2)
        self.kind_buttons: dict[str, Any] = {}
        for kind in ("All", *FLAG_KINDS):
            button = filter_chip(flag_text(kind))
            button.setChecked(kind == "All")
            button.clicked.connect(lambda _c, k=kind: self._on_kind_clicked(k))
            self.kind_buttons[kind] = button
            filters.addWidget(button)
        filters.addSpacing(theme.SPACING)
        self.file_combo = QComboBox()
        self.file_combo.setMinimumWidth(160)
        self.file_combo.currentIndexChanged.connect(lambda *_: self._apply_filter())
        filters.addWidget(self.file_combo)
        filters.addStretch(1)
        self.count_label = hint("")
        # No wrapping and no elision: a stretch before it means the row can always give it
        # the width it asks for, and without this it wrapped to "7661 of 7661 flag / events"
        # and sat against the export button as though the two were one broken control.
        self.count_label.setWordWrap(False)
        self.count_label.setSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Preferred)
        filters.addWidget(self.count_label)
        filters.addSpacing(theme.SPACING)
        self.export_button = secondary_button("Export CSV")
        self.export_button.clicked.connect(self.export_csv)
        filters.addWidget(self.export_button)
        outer.addLayout(filters)

        self.model = FlagEventModel(self)
        self.table = review_table()
        self.table.setModel(self.model)
        self.table.setItemDelegateForColumn(2, FlagChipDelegate(self.table))
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        outer.addWidget(self.table, 1)

        self.empty = EmptyState(
            "No flags match this filter",
            "",
        )
        outer.addWidget(self.empty)
        self.empty.setVisible(False)

    def set_result(self, result) -> None:
        self._events = flag_events(result) if result is not None else []
        self.file_combo.blockSignals(True)
        self.file_combo.clear()
        self.file_combo.addItem("All files", "")
        for file_id in dict.fromkeys(event.file_id for event in self._events):
            self.file_combo.addItem(file_id, file_id)
        self.file_combo.blockSignals(False)
        self._apply_filter()

    # ---- filtering ------------------------------------------------------------------------

    def selected_kinds(self) -> set[str]:
        """The flag kinds currently shown. Empty set means every kind ("All")."""
        if self.kind_buttons["All"].isChecked():
            return set()
        return {kind for kind in FLAG_KINDS if self.kind_buttons[kind].isChecked()}

    def _on_kind_clicked(self, kind: str) -> None:
        if kind == "All":
            for name, button in self.kind_buttons.items():
                button.setChecked(name == "All")
        else:
            # Any specific chip turns All off. Turning the last one off returns to All,
            # because an empty filter row showing an empty table is a dead end.
            self.kind_buttons["All"].setChecked(False)
            if not any(self.kind_buttons[k].isChecked() for k in FLAG_KINDS):
                self.kind_buttons["All"].setChecked(True)
        self._apply_filter()

    def filtered_events(self) -> list[FlagEvent]:
        kinds = self.selected_kinds()
        wanted_file = self.file_combo.currentData() or ""
        return [
            event
            for event in self._events
            if (not kinds or event.kind in kinds)
            and (not wanted_file or event.file_id == wanted_file)
        ]

    def _apply_filter(self) -> None:
        rows = self.filtered_events()
        self.model.set_rows(rows)
        self.count_label.setText(
            f"{len(rows)} of {len(self._events)} flag event{'' if len(self._events) == 1 else 's'}"
        )
        self.table.setVisible(bool(rows))
        self.empty.setVisible(not rows and bool(self._events))

    # ---- export ----------------------------------------------------------------------------

    def export_csv(self, path: str | Path | None = None) -> Path | None:
        """Write the *filtered* rows (what is on screen) to a CSV."""
        if path is None:
            chosen, _filter = QFileDialog.getSaveFileName(
                self, "Export flag events", "qaqc_flags.csv", "CSV files (*.csv)"
            )
            if not chosen:
                return None
            path = chosen
        path = Path(path)
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["file", "analyte", "flag", "value"])
            for event in self.filtered_events():
                writer.writerow([event.file_id, event.analyte, event.kind, event.context])
        return path


# --------------------------------------------------------------------------------------
# The page
# --------------------------------------------------------------------------------------


class ReviewPage(WorkflowPage):
    step = STEP_REVIEW
    #: Full-bleed: a chromatogram and a twenty-column heatmap are the two things in this
    #: application that genuinely want every pixel of a wide monitor.
    max_content_width = 0

    def build(self) -> None:
        self._charts = None
        #: The result the built tabs were filled from, compared by identity. ``id()`` would
        #: not do: CPython reuses the id of a freed object, so a run replaced by a second one
        #: could silently keep the first one's tabs.
        self._shown_result: Any = None
        #: The result :meth:`summary` last computed for, and what it computed. Same identity
        #: rule as ``_shown_result``.
        self._summary_for: Any = None
        self._summary: RunSummary | None = None
        self._tabs: dict[int, QWidget] = {}
        self.cache = ChromatogramCache()

        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)
        self._hosts: list[QWidget] = []
        for title in TAB_TITLES:
            host = QWidget()
            host.setObjectName("ReviewTab")
            layout = QVBoxLayout(host)
            layout.setContentsMargins(0, 0, 0, 0)
            layout.setSpacing(0)
            self._hosts.append(host)
            self.tabs.addTab(host, title)
        self.tabs.currentChanged.connect(self._on_tab_changed)
        self.content.addWidget(self.tabs, 1)

        self.empty = EmptyState(
            "No results yet",
            "Process the run first.",
        )
        self.content.addWidget(self.empty)
        self.refresh()

    # ---- lazy tab construction ------------------------------------------------------------

    def charts(self):
        """:mod:`quant209.app.charts`, imported the first time a drawing tab is built."""
        if self._charts is None:
            from quant209.app import charts

            self._charts = charts
        return self._charts

    def _build_tab(self, index: int) -> QWidget:
        if index == TAB_OVERVIEW:
            return OverviewTab()
        if index == TAB_QAQC:
            return QaqcTab()
        charts = self.charts()
        if index == TAB_SAMPLES:
            return SamplesTab(charts, self.cache)
        if index == TAB_CALIBRATION:
            return CalibrationTab(charts, self.cache)
        return RecoveriesTab(charts, self.cache)

    def ensure_tab(self, index: int) -> QWidget | None:
        """Build tab ``index`` if it has not been built, and fill it from the result."""
        result = self.state.result
        if result is None:
            return None
        widget = self._tabs.get(index)
        if widget is None:
            widget = self._build_tab(index)
            self._hosts[index].layout().addWidget(widget)
            self._tabs[index] = widget
            widget.set_result(result)
        return widget

    def _discard_tabs(self) -> None:
        for widget in self._tabs.values():
            widget.setParent(None)
            widget.deleteLater()
        self._tabs.clear()
        self.cache.clear()

    def _on_tab_changed(self, index: int) -> None:
        self.ensure_tab(index)

    # ---- lifecycle ---------------------------------------------------------------------------

    def on_activated(self) -> None:
        self.ensure_tab(self.tabs.currentIndex())

    def refresh(self) -> None:
        result = self.state.result
        has_result = result is not None
        self.tabs.setVisible(has_result)
        self.empty.setVisible(not has_result)
        if result is not self._shown_result:
            self._shown_result = result
            self._discard_tabs()
            if has_result:
                self.ensure_tab(self.tabs.currentIndex())
        self.status_changed.emit()

    def summary(self) -> RunSummary | None:
        """The Overview's numbers, whether or not that tab has been built.

        Cached on the result's *identity*, the same rule ``_shown_result`` uses and for the
        same reason (CPython reuses a freed object's ``id()``, so ``id()`` would not do).
        Every navigation to this page asks twice (once through ``refresh`` and once through
        the window's header chip) and :func:`summarise` costs ~15 ms on a 22-file run, about
        half of it re-deriving each calibration injection's response factors. Paying that once
        per result rather than once per question keeps a step change instant, and a
        ``RunResult`` is replaced rather than mutated, so there is nothing to go stale.
        """
        result = self.state.result
        if result is None:
            return None
        if result is not self._summary_for:
            self._summary = summarise(result)
            self._summary_for = result
        return self._summary

    def page_status(self) -> tuple[str, str] | None:
        summary = self.summary()
        if summary is None:
            return None
        tone, headline, _why = verdict(summary)
        return ({"error": "Fit failed", "warning": "Review required", "success": "OK"}[tone], tone)
