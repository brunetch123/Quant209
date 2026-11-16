"""The repeatable validation report generator, plus the joining and metric logic that
``tools/bench_areas.py``, ``tools/bench_quant.py`` and ``tools/master_bench.py`` share.

This module has two jobs, deliberately kept in one place rather than split, because the
second cannot be trusted without the first:

1. **The reusable joining and metric functions** (:func:`join_areas`, :func:`join_masses`,
   :func:`join_recoveries`, :func:`join_rrfs`, :func:`metrics`, and the row dataclasses and
   RPD/detectable-subset definitions they share). Every figure quoted in ``docs/METHODS.md``
   and ``docs/validation/`` comes from these functions, as does every number the application
   and the CLI show, so there is exactly one implementation of what "agreement" means.

2. **generate_validation_report**, a compact report for one run. It writes five CSVs an
   analyst or a future regression test can diff, three vector figures, and a ``summary.md``
   whose headline numbers are recomputed from those same CSVs rather than hand-typed, with a
   verdict against the acceptance criteria.

``RPD = (auto - manual) / ((auto + manual) / 2)`` on every comparison here except recovery
(already a fraction, so scored as absolute error) and the log-log regression figure, which
plots the values themselves rather than their difference.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: this module may run with no display, and imports pyplot below
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from quant209.core.config import Analyte  # noqa: E402
from quant209.core.peaks import PeakResult  # noqa: E402
from quant209.core.pipeline import RunResult  # noqa: E402
from quant209.core.standards import HOMOLOG_SURROGATE  # noqa: E402
from quant209.report.truth import RunTruth  # noqa: E402

# --------------------------------------------------------------------------------------
# Shared constants
# --------------------------------------------------------------------------------------

#: Agreement band every RPD-based table counts rows inside.
BAND = 0.25

#: A detectable row within this multiple of the threshold is still a trace peak: the
#: analyst's own manual boundaries there are worth a few tens of percent either way, so a
#: disagreement is not evidence about the engine. Used only to *classify* misses, never to
#: exclude them from the headline metrics.
NEAR_LOQ_MULTIPLE = 5.0

#: mono..deca, for ordering per-homolog tables/figures chromatographically.
HOMOLOG_ORDER = (
    "mono", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona", "deca",
)

#: Manual area (counts·s) below which a row is treated as noise, not a measurement.
DEFAULT_MIN_AREA = 1000.0

#: Analyst mass (ng) below which a row is treated as noise, not a measurement.
DEFAULT_MIN_MASS = 1.0

#: Acceptance criteria (see docs/METHODS.md section 6.1).
GATE_MEDIAN_RPD = 0.10
GATE_WITHIN_FRAC = 0.80
GATE_RECOVERY_MAE = 0.05


# --------------------------------------------------------------------------------------
# Joined rows
# --------------------------------------------------------------------------------------


@dataclass
class AreaRow:
    """One joined (sample, analyte) peak-area comparison."""

    sample: str
    file_type: str
    label: str
    homolog: str
    role: str
    manual: float
    auto: float
    flags: str
    method: str
    sn: float
    rt_apex: float

    @property
    def rpd(self) -> float:
        """Symmetric relative difference, and 0 when both areas are zero."""
        mean = (self.auto + self.manual) / 2.0
        return (self.auto - self.manual) / mean if mean else 0.0


@dataclass
class MassRow:
    """One joined (sample, analyte) mass comparison, in nanograms."""

    sample: str
    file_type: str
    label: str
    homolog: str
    role: str
    manual: float
    auto: float
    #: False when the pipeline produced no mass at all for this analyte in this file, either
    #: for want of a response factor or because the internal standard was unusable. ``auto``
    #: is then 0.0 by convention.
    quantified: bool
    flags: str
    method: str
    rt_apex: float
    detected: bool

    @property
    def rpd(self) -> float:
        mean = (self.auto + self.manual) / 2.0
        return (self.auto - self.manual) / mean if mean else 0.0


#: Same-file check (see :func:`verify_pairing`). A sheet disagrees with its own file when the
#: median |log10(ours / analyst)| over the labelled standards exceeds this (4.7 %). On the
#: archive, 99 % of sheets agree with their own file to within 0.0054 (1.2 %).
PAIRING_FLOOR = 0.02
#: A crossing needs the partner file to match at least this many times more closely.
PAIRING_MARGIN = 4.0
#: Fewest labelled standards measured on both sides for the fingerprint to say anything.
PAIRING_MIN_STANDARDS = 4


#: How close a workbook's group calibration mass must sit to 1/N of the certified summed mass,
#: N being the group's congener count, to be read as one congener's mass entered for the
#: whole group (see :func:`rebase_group_masses`).
GROUP_REBASE_TOLERANCE = 0.10


def rebase_group_masses(truth, analytes) -> tuple[pd.DataFrame, list[str]]:
    """The analyst's group masses, on the summed calibration basis Quant 209 reports on.

    Returns the masses and the notes to record. A co-eluting group's peak holds every congener
    in it, so the calibration standard's mass for it is the sum of theirs. A workbook can
    instead enter one member's mass for a group, which puts its group masses on a per-congener
    basis, a half or a third of the summed figure, and a comparison between the two would
    measure the convention rather than the integration.

    Where a workbook's ``cal_info`` mass for a group sits within
    :data:`GROUP_REBASE_TOLERANCE` of 1/N of the certified mass, the masses of the samples read
    from that workbook are multiplied by the ratio. The evidence is the workbook's own
    calibration masses; no scored value is consulted.
    """
    masses = truth.masses.copy()
    cal = getattr(truth, "cal_masses", None)
    if cal is None or cal.empty or masses.empty:
        return masses, []
    certified = {a.label: float(a.cal_mass) for a in analytes if a.role == "native"}
    factors: dict[tuple[str, str], float] = {}
    for label in cal.index:
        if "+" not in label or label not in certified:
            continue
        members = label.count("+") + 1
        for book, entered in cal.loc[label].items():
            if entered is None or not np.isfinite(entered) or entered <= 0:
                continue
            ratio = certified[label] / float(entered)
            if abs(ratio - members) / members <= GROUP_REBASE_TOLERANCE:
                factors[(str(book), label)] = ratio
    if not factors:
        return masses, []
    used: set[str] = set()
    for sample, book in (getattr(truth, "sample_book", None) or {}).items():
        if sample not in masses.columns:
            continue
        for (b, label), ratio in factors.items():
            if b == book and label in masses.index:
                masses.loc[label, sample] = masses.loc[label, sample] * ratio
                used.add(b)
    if not used:  # every sample came from a workbook already on the summed basis
        return masses, []
    groups = sorted({label for b, label in factors if b in used})
    return masses, [f"analyst masses of {len(groups)} co-eluting groups re-based to the summed "
                    f"calibration mass in {', '.join(sorted(used))}"]


@dataclass
class PairingRow:
    """How one analyst sheet was paired with a data file, and the evidence for it."""

    sheet: str
    file: str               # the data file this sheet is scored against
    own_distance: float     # fingerprint distance to the file sharing the sheet's name
    best_other: str         # the closest other file in the run
    best_distance: float
    verdict: str            # "same file" | "re-paired" | "unverified" | "no fingerprint"


def _standards_distance(theirs: pd.Series, ours: pd.Series) -> float:
    both = [(t, o) for t, o in zip(theirs.to_numpy(float), ours.to_numpy(float))
            if np.isfinite(t) and np.isfinite(o) and t > 0 and o > 0]
    if len(both) < PAIRING_MIN_STANDARDS:
        return float("nan")
    return float(np.median([abs(np.log10(o / t)) for t, o in both]))


def verify_pairing(
    result: RunResult, truth_areas: pd.DataFrame
) -> tuple[dict[str, str], list[PairingRow]]:
    """Check that every analyst sheet is compared with the data file it came from.

    A sheet and a file are matched by name, and a name can be wrong: a vial swapped at the
    bench, or a column pasted under the wrong heading. The evidence used here is the labelled
    standards (13C surrogates and internal standards), spiked into every file. For the same
    file our integration and the analyst's agree to a median 0.4 %, and different files
    differ by several percent. The native congeners, the very values being scored, are never
    consulted, so the check cannot choose a pairing for its agreement.

    A sheet is re-paired only when it and another sheet are each a decisive match for the
    other's file (:data:`PAIRING_MARGIN`), and each disagrees with its own
    (:data:`PAIRING_FLOOR`): a mutual crossing, which is what a swap looks like. A sheet that
    disagrees with its own file without such a partner is reported as ``unverified`` and left
    where it is, because the likelier explanation is an integration failure on one side.

    Returns ``(aliases, rows)``: ``aliases`` maps a data file id to the analyst sheet it is
    scored against, for :func:`join_areas`, :func:`join_masses` and :func:`join_recoveries`.
    """
    labelled = [a.label for a in result.analytes if a.role != "native"]
    ours = result.areas
    files = [f for f in result.reported_ids if f in ours.columns]
    sheets = [c for c in truth_areas.columns if c in files and truth_areas[c].notna().any()]
    labels = [l for l in labelled if l in ours.index and l in truth_areas.index]
    if not labels or not sheets:
        return {}, []
    theirs = truth_areas.loc[labels]
    mine = ours.loc[labels]
    dist = {(s, f): _standards_distance(theirs[s], mine[f]) for s in sheets for f in files}
    own = {s: dist[(s, s)] for s in sheets}
    best: dict[str, tuple[float, str]] = {}
    for s in sheets:
        others = [(d, f) for (ss, f), d in dist.items() if ss == s and f != s and np.isfinite(d)]
        best[s] = min(others) if others else (float("nan"), "")

    def crossed(s: str) -> bool:
        d_own, (d_best, partner) = own[s], best[s]
        return (np.isfinite(d_own) and np.isfinite(d_best) and d_own >= PAIRING_FLOOR
                and d_best * PAIRING_MARGIN <= d_own and partner in own)

    aliases: dict[str, str] = {}
    rows: list[PairingRow] = []
    for s in sheets:
        d_own, (d_best, partner) = own[s], best[s]
        mutual = crossed(s) and crossed(partner) and best[partner][1] == s
        if mutual:
            aliases[partner] = s
            verdict, file = "re-paired", partner
        elif not np.isfinite(d_own):
            verdict, file = "no fingerprint", s
        elif d_own < PAIRING_FLOOR:
            verdict, file = "same file", s
        else:
            verdict, file = "unverified", s
        rows.append(PairingRow(s, file, d_own, partner, d_best, verdict))
    return aliases, rows


@dataclass
class RecoveryRow:
    """One joined (sample, homolog) surrogate-recovery comparison, as a fraction."""

    sample: str
    homolog: str
    manual: float
    auto: float

    @property
    def error(self) -> float:
        return self.auto - self.manual


@dataclass
class RrfRow:
    """One analyte's response factor, measured against the analyst's, with its evidence."""

    label: str
    homolog: str
    ours: float
    analyst: float
    our_cal_area: float
    analyst_cal_area: float
    flags: str
    method: str

    @property
    def rpd(self) -> float:
        mean = (self.ours + self.analyst) / 2.0
        return (self.ours - self.analyst) / mean if mean else 0.0

    @property
    def area_rpd(self) -> float:
        mean = (self.our_cal_area + self.analyst_cal_area) / 2.0
        return (self.our_cal_area - self.analyst_cal_area) / mean if mean else 0.0


@dataclass
class RefRow:
    """One (reference injection, homolog) surrogate-mass comparison, in nanograms."""

    ref_id: str
    workbook: str
    homolog: str
    manual: float
    auto: float

    @property
    def rpd(self) -> float:
        mean = (self.auto + self.manual) / 2.0
        return (self.auto - self.manual) / mean if mean else 0.0


@dataclass
class InjectionMatch:
    """Which calibration injection an analyst workbook's ``cal`` sheet actually holds."""

    workbook: str
    cal_id: str
    scores: dict[str, float]  # cal file id -> median |RPD| of native calibration areas

    @property
    def score(self) -> float:
        return self.scores.get(self.cal_id, float("nan"))

    @property
    def runner_up(self) -> float:
        others = sorted(v for k, v in self.scores.items() if k != self.cal_id)
        return others[0] if others else float("nan")


# --------------------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------------------


@dataclass
class Metrics:
    """Agreement statistics for one subset of joined area or mass rows.

    One dataclass for both: ``n_unquantified`` is meaningless for areas (every integrated
    peak has an area) and stays 0 there by construction, since :class:`AreaRow` carries no
    ``quantified`` attribute for :func:`metrics` to inspect.
    """

    n_joined: int
    n_detectable: int
    median_abs_rpd: float
    within_band: float
    total_ratio: float
    median_rpd: float  # signed, which separates a bias from symmetric scatter
    n_unquantified: int = 0  # detectable rows with no computed mass at all (masses only)


def metrics(rows: list[AreaRow] | list[MassRow], threshold: float) -> Metrics:
    """RPD-based agreement statistics over ``rows`` whose *manual* value clears ``threshold``."""
    detectable = [r for r in rows if r.manual >= threshold]
    if not detectable:
        nan = float("nan")
        return Metrics(len(rows), 0, nan, nan, nan, nan, 0)
    rpds = np.array([r.rpd for r in detectable])
    manual = np.array([r.manual for r in detectable])
    auto = np.array([r.auto for r in detectable])
    return Metrics(
        n_joined=len(rows),
        n_detectable=len(detectable),
        median_abs_rpd=float(np.median(np.abs(rpds))),
        within_band=float(np.mean(np.abs(rpds) <= BAND)),
        total_ratio=float(auto.sum() / manual.sum()) if manual.sum() else float("nan"),
        median_rpd=float(np.median(rpds)),
        n_unquantified=sum(1 for r in detectable if not getattr(r, "quantified", True)),
    )


@dataclass
class RecoveryMetrics:
    n: int
    mae: float
    median_error: float  # signed
    max_abs_error: float
    within_005: float


def recovery_metrics(rows: list[RecoveryRow]) -> RecoveryMetrics:
    if not rows:
        nan = float("nan")
        return RecoveryMetrics(0, nan, nan, nan, nan)
    errors = np.array([r.error for r in rows])
    return RecoveryMetrics(
        n=len(rows),
        mae=float(np.mean(np.abs(errors))),
        median_error=float(np.median(errors)),
        max_abs_error=float(np.max(np.abs(errors))),
        within_005=float(np.mean(np.abs(errors) <= 0.05)),
    )


@dataclass
class RrfMetrics:
    n: int
    median_abs_rpd: float
    within_band: float
    n_off: int  # outside ±25%
    n_high: int  # ...of which read high (a calibration peak that integrated short)
    n_low: int
    median_abs_area_rpd: float  # the same comparison on the calibration areas behind the RRFs


def rrf_metrics(rows: list[RrfRow]) -> RrfMetrics:
    if not rows:
        nan = float("nan")
        return RrfMetrics(0, nan, nan, 0, 0, 0, nan)
    rpds = np.array([r.rpd for r in rows])
    areas = np.array([r.area_rpd for r in rows if np.isfinite(r.area_rpd)])
    return RrfMetrics(
        n=len(rows),
        median_abs_rpd=float(np.median(np.abs(rpds))),
        within_band=float(np.mean(np.abs(rpds) <= BAND)),
        n_off=int((np.abs(rpds) > BAND).sum()),
        n_high=int((rpds > BAND).sum()),
        n_low=int((rpds < -BAND).sum()),
        median_abs_area_rpd=float(np.median(np.abs(areas))) if areas.size else float("nan"),
    )


def median_abs_rpd(ours: np.ndarray, theirs: np.ndarray) -> float:
    """Median symmetric relative difference between two aligned arrays of positives."""
    mean = (ours + theirs) / 2.0
    valid = mean > 0
    if not valid.any():
        return float("nan")
    return float(np.median(np.abs((ours[valid] - theirs[valid]) / mean[valid])))


# --------------------------------------------------------------------------------------
# Miss classification
# --------------------------------------------------------------------------------------


def classify_area_miss(row: AreaRow, min_area: float) -> str:
    """Name the failure class of an area row that fell outside the agreement band.

    Tested in priority order. A row is reported under the first class that explains it.
    """
    if "absent" in row.flags:
        return "presence test rejected (shoulder on a dominant neighbour, or no peak)"
    if row.file_type == "blank":
        return "blank file (carryover / trace-level integration)"
    if row.manual < NEAR_LOQ_MULTIPLE * min_area:
        return "near-LOQ trace peak"
    return "other"


def classify_mass_miss(row: MassRow, min_mass: float) -> str:
    """Name the failure class of a mass row that fell outside the agreement band.

    The first two classes are quantification-chain failures that cannot happen at the area
    layer. The rest are the area layer's own categories, inherited through the chain.
    """
    if not row.quantified:
        return "no mass computed (no response factor, or unusable internal standard)"
    if not row.detected:
        return "presence test rejected (shoulder on a dominant neighbour, or no peak)"
    if row.file_type == "blank":
        return "blank file (carryover / trace-level integration)"
    if row.manual < NEAR_LOQ_MULTIPLE * min_mass:
        return "near-LOQ trace mass"
    return "other"


# --------------------------------------------------------------------------------------
# Joins
# --------------------------------------------------------------------------------------


def join_areas(
    peaks_by_file: dict[str, list[PeakResult]],
    file_types: dict[str, str],
    analytes: list[Analyte],
    truth_areas: pd.DataFrame,
    aliases: dict[str, str] | None = None,
) -> list[AreaRow]:
    """Join already-integrated peak areas (keyed by file id) to the analyst's manual areas.

    Pure: this function does no I/O and decides nothing about which files to integrate or
    score. That is the caller's job. ``tools/bench_areas.py`` skips integrating a file with
    no truth column at all, and :func:`generate_validation_report` integrates every reported
    file the pipeline already produced. ``peaks_by_file`` is iterated in dict order, so
    callers that care about row order should build it in that order.

    ``aliases`` maps a file id to the truth column it should be scored against (see the
    module-level callers for why this exists: a provably mis-paired ground-truth column).
    """
    aliases = aliases or {}
    by_label = {a.label: a for a in analytes}
    rows: list[AreaRow] = []
    for file_id, peaks in peaks_by_file.items():
        truth_id = aliases.get(file_id, file_id)
        if truth_id not in truth_areas.columns:
            continue
        manual_col = truth_areas[truth_id]
        shown = file_id if truth_id == file_id else f"{file_id}->{truth_id}"
        file_type = file_types.get(file_id, "")
        for peak in peaks:
            if peak.label not in manual_col.index:
                continue
            manual = manual_col.loc[peak.label]
            if manual is None or not np.isfinite(float(manual)):
                continue
            analyte = by_label.get(peak.label)
            if analyte is None:
                continue
            rows.append(
                AreaRow(
                    sample=shown,
                    file_type=file_type,
                    label=peak.label,
                    homolog=analyte.homolog or "labeled",
                    role=analyte.role,
                    manual=float(manual),
                    auto=float(peak.area),
                    flags=",".join(peak.flags) or "-",
                    method=peak.method,
                    sn=float(peak.sn),
                    rt_apex=float(peak.rt_apex),
                )
            )
    return rows


def join_masses(
    result: RunResult, truth_masses: pd.DataFrame, aliases: dict[str, str] | None = None
) -> tuple[list[MassRow], list[str], list[str]]:
    """Join ``RunResult.masses`` to the analyst's uncorrected masses, file by file.

    Only ``sample`` and ``blank`` files are joined (``result.reported_ids``): the analyst's
    workbook holds no worksheet for the calibration or reference injections, so there is
    nothing on the other side of the join. Returns ``(rows, truth_only, file_only)``.
    """
    aliases = aliases or {}
    by_label = {a.label: a for a in result.analytes}
    populated = {c for c in truth_masses.columns if truth_masses[c].notna().any()}
    rows: list[MassRow] = []
    file_only: list[str] = []
    scored: set[str] = set()

    for file_id in sorted(result.reported_ids):
        truth_id = aliases.get(file_id, file_id)
        if truth_id not in populated:
            file_only.append(file_id)
            continue
        record = result.files[file_id]
        if not record.usable:
            continue
        scored.add(truth_id)
        manual_col = truth_masses[truth_id]
        auto_col = result.masses[file_id]
        peaks = {p.label: p for p in record.peaks}
        shown = file_id if truth_id == file_id else f"{file_id}->{truth_id}"
        for label, manual in manual_col.items():
            analyte = by_label.get(label)
            if analyte is None or manual is None or not np.isfinite(float(manual)):
                continue
            peak = peaks.get(label)
            if peak is None:
                continue
            auto = auto_col.get(label, np.nan)
            quantified = bool(np.isfinite(auto))
            rows.append(
                MassRow(
                    sample=shown,
                    file_type=record.file_type,
                    label=label,
                    homolog=analyte.homolog or "labeled",
                    role=analyte.role,
                    manual=float(manual),
                    auto=float(auto) if quantified else 0.0,
                    quantified=quantified,
                    flags=",".join(peak.flags) or "-",
                    method=peak.method,
                    rt_apex=float(peak.rt_apex),
                    detected=peak.detected,
                )
            )
    truth_only = sorted(populated - scored)
    return rows, truth_only, file_only


def join_recoveries(
    result: RunResult, truth_recoveries: pd.DataFrame, aliases: dict[str, str] | None = None
) -> tuple[list[RecoveryRow], list[str]]:
    """Join ``RunResult.recoveries`` to the workbook Summary sheet's Recoveries block."""
    aliases = aliases or {}
    rows: list[RecoveryRow] = []
    unmatched: list[str] = []
    for file_id in sorted(result.reported_ids):
        truth_id = aliases.get(file_id, file_id)
        if truth_id not in truth_recoveries.columns:
            continue
        shown = file_id if truth_id == file_id else f"{file_id}->{truth_id}"
        for homolog in HOMOLOG_ORDER:
            if homolog not in truth_recoveries.index or homolog not in result.recoveries.index:
                continue
            manual = truth_recoveries.at[homolog, truth_id]
            auto = result.recoveries.at[homolog, file_id]
            if manual is None or not np.isfinite(float(manual)):
                continue
            if not np.isfinite(float(auto)):
                unmatched.append(f"{shown}/{homolog}")
                continue
            rows.append(
                RecoveryRow(sample=shown, homolog=homolog, manual=float(manual), auto=float(auto))
            )
    return rows, unmatched


def match_areas(
    cal_areas: pd.DataFrame,
    ours: dict[str, pd.Series],
    natives: set[str],
    known: set[str],
) -> list[InjectionMatch]:
    """The scoring half of :func:`match_injections`, over plain area columns.

    Split out so a sweep over integration parameters can identify the same workbook ->
    injection pairing without building a whole ``RunResult``.
    """
    matches: list[InjectionMatch] = []
    for book in cal_areas.columns:
        column = cal_areas[book]
        labels = [
            label
            for label in column.index
            if label in natives
            and label in known
            and np.isfinite(column.get(label, np.nan))
            and float(column[label]) > 0
        ]
        if len(labels) < 10:  # an unfilled or unrecognised sheet: nothing to match on
            continue
        theirs = column.loc[labels].to_numpy(dtype=float)
        scores = {
            cid: median_abs_rpd(series.loc[labels].to_numpy(dtype=float), theirs)
            for cid, series in ours.items()
        }
        scores = {cid: s for cid, s in scores.items() if np.isfinite(s)}
        if not scores:
            continue
        best = min(scores, key=lambda cid: scores[cid])
        matches.append(InjectionMatch(workbook=book, cal_id=best, scores=scores))
    return matches


def match_injections(result: RunResult, cal_areas: pd.DataFrame) -> list[InjectionMatch]:
    """Identify which calibration injection each workbook's ``cal`` sheet holds.

    Scored on *native* calibration areas alone. They are what the response factors are made
    of, and they dominate the discriminating signal. A workbook whose sheet was never filled
    in is skipped rather than matched to whichever injection happens to score least badly
    against a column of zeros.
    """
    natives = {a.label for a in result.analytes if a.role == "native"}
    usable_cal = [cid for cid in result.cal_ids if result.files[cid].usable]
    return match_areas(
        cal_areas,
        {cid: result.areas[cid] for cid in usable_cal},
        natives,
        set(result.areas.index),
    )


def analyst_rrf_spread(cal_rrfs: pd.DataFrame, natives: set[str]) -> float:
    """Median |RPD| between the analyst's own workbooks for the same analyte's RRF.

    Every workbook of a run quantified against a different calibration injection, so this is
    the analyst's own injection-to-injection reproducibility. It is the floor below which a
    disagreement stops being a statement about this pipeline's integration.
    """
    books = [b for b in cal_rrfs.columns if cal_rrfs[b].notna().any()]
    rows = [label for label in cal_rrfs.index if label in natives]
    if len(books) < 2 or not rows:
        return float("nan")
    frame = cal_rrfs.loc[rows, books]
    diffs: list[np.ndarray] = []
    for i, one in enumerate(books):
        for other in books[i + 1:]:
            a = frame[one].to_numpy(dtype=float)
            b = frame[other].to_numpy(dtype=float)
            ok = np.isfinite(a) & np.isfinite(b) & ((a + b) > 0)
            if ok.any():
                diffs.append(np.abs((a[ok] - b[ok]) / ((a[ok] + b[ok]) / 2.0)))
    return float(np.median(np.concatenate(diffs))) if diffs else float("nan")


def join_rrfs(
    result: RunResult,
    cal_rrfs: pd.DataFrame,
    cal_areas: pd.DataFrame,
    matches: list[InjectionMatch],
) -> list[RrfRow]:
    """Join the run's shipped response factors to the analyst's, analyte by analyte.

    The analyst side is the **median across the run's workbooks**. The pipeline side is an
    average over every calibration injection, and each workbook's ``cal`` sheet is one
    injection.
    """
    natives = {a.label: a for a in result.analytes if a.role == "native"}
    books = [b for b in cal_rrfs.columns if cal_rrfs[b].notna().any()]
    if not books:
        return []
    analyst = cal_rrfs[books].median(axis=1)
    anchor = matches[0].cal_id if matches else (result.cal_ids[0] if result.cal_ids else "")
    peaks = {p.label: p for p in result.files[anchor].peaks} if anchor in result.files else {}
    anchor_book = matches[0].workbook if matches else (books[0] if books else "")

    rows: list[RrfRow] = []
    for label, analyte in natives.items():
        ours = result.calibration.rrf.get(label)
        theirs = analyst.get(label, np.nan)
        if ours is None or not np.isfinite(float(theirs)) or float(theirs) <= 0:
            continue
        peak = peaks.get(label)
        rows.append(
            RrfRow(
                label=label,
                homolog=analyte.homolog or "labeled",
                ours=float(ours),
                analyst=float(theirs),
                our_cal_area=float(result.areas.at[label, anchor])
                if anchor in result.areas.columns and label in result.areas.index
                else float("nan"),
                analyst_cal_area=float(cal_areas.at[label, anchor_book])
                if anchor_book in cal_areas.columns and label in cal_areas.index
                else float("nan"),
                flags=",".join(peak.flags) if peak is not None and peak.flags else "-",
                method=peak.method if peak is not None else "-",
            )
        )
    return rows


def per_injection_rrf(
    result: RunResult, cal_rrfs: pd.DataFrame, matches: list[InjectionMatch]
) -> list[tuple[str, str, int, float, int]]:
    """One injection's RRF against the analyst's from the workbook that used it.

    Like for like, with no averaging on either side. ``quant209.core.quantify`` is imported
    at call time so that callers which never sweep do not pay for the import.
    """
    from quant209.core.quantify import build_calibration

    natives = {a.label for a in result.analytes if a.role == "native"}
    out: list[tuple[str, str, int, float, int]] = []
    for match in matches:
        if match.workbook not in cal_rrfs.columns or not result.files[match.cal_id].usable:
            continue
        model = build_calibration(result.files[match.cal_id].peaks, result.analytes)
        column = cal_rrfs[match.workbook]
        labels = [
            label
            for label in column.index
            if label in natives
            and label in model.rrf
            and np.isfinite(column.get(label, np.nan))
            and float(column[label]) > 0
        ]
        if not labels:
            continue
        ours = np.array([model.rrf[label] for label in labels])
        theirs = column.loc[labels].to_numpy(dtype=float)
        mean = (ours + theirs) / 2.0
        rpds = np.abs((ours - theirs) / np.where(mean != 0, mean, 1.0))
        out.append(
            (match.workbook, match.cal_id, len(labels), float(np.median(rpds)),
             int((rpds > BAND).sum()))
        )
    return out


def join_refs(
    result: RunResult, ref_masses: pd.DataFrame, workbook_samples: dict[str, list[str]]
) -> list[RefRow]:
    """Join the reference injections' surrogate masses to the analyst's ``ref`` sheets."""
    rows: list[RefRow] = []
    for book, samples in sorted(workbook_samples.items()):
        if book not in ref_masses.columns:
            continue
        refs = {result.reference_for[s] for s in samples if s in result.reference_for}
        if len(refs) != 1:
            continue
        ref_id = refs.pop()
        for homolog, surrogate in HOMOLOG_SURROGATE.items():
            manual = ref_masses.at[homolog, book] if homolog in ref_masses.index else np.nan
            if not np.isfinite(float(manual)):
                continue
            if surrogate not in result.masses.index or ref_id not in result.masses.columns:
                continue
            auto = result.masses.at[surrogate, ref_id]
            if not np.isfinite(float(auto)):
                continue
            rows.append(
                RefRow(ref_id=ref_id, workbook=book, homolog=homolog,
                       manual=float(manual), auto=float(auto))
            )
    return rows


# --------------------------------------------------------------------------------------
# Per-run validation report
# --------------------------------------------------------------------------------------


@dataclass
class ValidationReport:
    """What one call to :func:`generate_validation_report` measured and wrote."""

    run: str
    out_dir: Path
    area: Metrics
    mass: Metrics
    recovery: RecoveryMetrics
    rrf: RrfMetrics
    n_bias_rows: int
    area_pass: bool
    mass_pass: bool
    recovery_pass: bool
    passed: bool
    paths: dict[str, Path] = field(default_factory=dict)


def likely_cause(
    row: MassRow, min_mass: float, rrf_off: set[str], no_rrf_labels: set[str] = frozenset()
) -> str:
    """Best-guess reason a detectable native mass row fell outside ±25%.

    Checked in priority order. A row is reported under the first class that explains it.

    1. No mass at all, because this analyte never got a response factor this run. The
       cal-peak guard refused it in every calibration injection (``no_rrf_labels``). This is
       a calibration-chain failure.
    2. No mass for any other reason, meaning this file's internal standard was unusable.
       That is a per-sample quantification failure rather than a calibration one, so it
       falls to "other" rather than overloading "cal-anomaly" with a second failure mode.
    3. The presence test rejected the peak, a shoulder split off a dominant neighbour.
    4. A near-LOQ trace mass, where the analyst's own boundaries are worth tens of percent.
    5. This analyte's *own* response factor disagreed with the analyst's by more than the
       band (``rrf_off``, from :func:`join_rrfs`). An RRF error multiplies every sample of
       that analyte, so it explains the miss without anything sample-specific being wrong.
    6. Otherwise, unexplained by any of the above.
    """
    if not row.quantified:
        return "cal-anomaly" if row.label in no_rrf_labels else "other"
    if not row.detected:
        return "presence-rejected"
    if row.manual < NEAR_LOQ_MULTIPLE * min_mass:
        return "near-LOQ"
    if row.label in rrf_off:
        return "cal-anomaly"
    return "other"


_CSV_FIELDS: dict[str, list[str]] = {
    "area": ["sample", "file_type", "label", "homolog", "role", "manual", "auto", "rpd",
             "flags", "method", "sn", "rt_apex"],
    "mass": ["sample", "file_type", "label", "homolog", "role", "manual", "auto", "rpd",
             "quantified", "flags", "method", "rt_apex", "detected"],
    "recovery": ["sample", "homolog", "manual", "auto", "error"],
    "rrf": ["label", "homolog", "ours", "analyst", "rpd", "our_cal_area", "analyst_cal_area",
            "area_rpd", "flags", "method"],
    "bias": ["sample", "label", "homolog", "manual", "auto", "rpd", "flags", "method",
             "clustered", "likely_cause"],
}


def _write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _style_axes(ax) -> None:
    """Drop the top and right spines and keep ticks at 10pt or larger."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=10)


def _fig_rpd_hist(rpds: list[float], path: Path) -> None:
    """Signed RPD distribution of the native, detectable **area** rows (counts·s)."""
    fig, ax = plt.subplots(figsize=(6, 4))
    values = np.asarray([v for v in rpds if np.isfinite(v)], dtype=float) * 100.0
    if values.size:
        ax.hist(values, bins=min(20, max(5, values.size)), color="#440154", edgecolor="white")
    ax.axvline(0.0, color="#21918c", linewidth=1.5)
    ax.set_xlabel("RPD (%), automated against manual peak area (counts·s)", fontsize=11)
    ax.set_ylabel("count", fontsize=11)
    ax.set_title("Native, detectable area RPD", fontsize=12)
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(path, format="pdf")
    plt.close(fig)


def _fig_regression(rows: list[MassRow], path: Path) -> None:
    """Log-log automated vs analyst **mass** (ng), with a fitted slope and R²."""
    fig, ax = plt.subplots(figsize=(5, 5))
    manual = np.array([r.manual for r in rows if r.manual > 0 and r.auto > 0], dtype=float)
    auto = np.array([r.auto for r in rows if r.manual > 0 and r.auto > 0], dtype=float)
    if manual.size:
        ax.scatter(manual, auto, s=22, color="#31688e", alpha=0.75, edgecolor="none")
        lo, hi = min(manual.min(), auto.min()), max(manual.max(), auto.max())
        ax.plot([lo, hi], [lo, hi], color="#35b779", linewidth=1.2, linestyle="--")
        if manual.size >= 2 and np.unique(manual).size >= 2:
            logm, loga = np.log10(manual), np.log10(auto)
            slope, intercept = np.polyfit(logm, loga, 1)
            pred = slope * logm + intercept
            ss_res = float(np.sum((loga - pred) ** 2))
            ss_tot = float(np.sum((loga - loga.mean()) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot else float("nan")
            ax.text(
                0.05, 0.95, f"slope = {slope:.3f}\nR² = {r2:.3f}",
                transform=ax.transAxes, fontsize=10, va="top",
            )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("manual mass (ng)", fontsize=11)
    ax.set_ylabel("automated mass (ng)", fontsize=11)
    ax.set_title("Automated vs analyst mass", fontsize=12)
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(path, format="pdf")
    plt.close(fig)


def _fig_homolog_agreement(rows: list[MassRow], path: Path) -> None:
    """Per-homolog median |RPD| bars, native detectable **mass** rows (ng)."""
    grouped: dict[str, list[float]] = {}
    for r in rows:
        grouped.setdefault(r.homolog, []).append(abs(r.rpd))
    order = [h for h in HOMOLOG_ORDER if h in grouped]
    fig, ax = plt.subplots(figsize=(7, 4))
    if order:
        values = [100.0 * float(np.median(grouped[h])) for h in order]
        colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(order)))
        ax.bar(order, values, color=colors)
    ax.axhline(100 * BAND, color="#440154", linewidth=1.0, linestyle=":")
    ax.set_xlabel("homolog", fontsize=11)
    ax.set_ylabel("median |RPD| (%), mass (ng)", fontsize=11)
    ax.set_title("Per-homolog mass agreement", fontsize=12)
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(path, format="pdf")
    plt.close(fig)


def _pct(x: float) -> str:
    return "n/a" if not np.isfinite(x) else f"{100 * x:.1f}%"


def _num(x: float, places: int = 4) -> str:
    return "n/a" if not np.isfinite(x) else f"{x:.{places}f}"


def _table(header: list[str], rows: list[list[str]]) -> list[str]:
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return out


def _render_summary(
    run_name: str,
    area_m: Metrics,
    mass_m: Metrics,
    recovery_m: RecoveryMetrics,
    rrf_m: RrfMetrics,
    n_bias: int,
    area_pass: bool,
    mass_pass: bool,
    recovery_pass: bool,
    passed: bool,
    min_area: float,
    min_mass: float,
) -> str:
    lines = [
        f"# Validation report: {run_name}",
        "",
        "Generated by `quant209.report.validation.generate_validation_report`. Every number "
        "below is recomputed from the CSVs written alongside this file "
        "(`area_agreement.csv`, `mass_agreement.csv`, `recovery_agreement.csv`, "
        "`rrf_agreement.csv`); `bias_diagnosis.csv` and the three figures characterize the "
        "detectable rows that fell outside +/-25%.",
        "",
        "## Areas (native analytes)",
        "",
    ]
    lines += _table(
        ["metric", "value"],
        [
            ["joined rows", f"{area_m.n_joined}"],
            ["detectable rows (manual area >= threshold)", f"{area_m.n_detectable}"],
            ["median |RPD|", _pct(area_m.median_abs_rpd)],
            ["within +/-25%", _pct(area_m.within_band)],
            ["area-weighted total ratio (sum auto / sum manual)", _num(area_m.total_ratio)],
            ["median signed RPD (bias)", _pct(area_m.median_rpd)],
        ],
    )
    lines += ["", f"Detectability threshold: manual area >= {min_area:g} counts·s.", ""]

    lines += ["## Masses (native analytes, uncorrected)", ""]
    lines += _table(
        ["metric", "value"],
        [
            ["joined rows", f"{mass_m.n_joined}"],
            ["detectable rows (manual mass >= threshold)", f"{mass_m.n_detectable}"],
            ["median |RPD|", _pct(mass_m.median_abs_rpd)],
            ["within +/-25%", _pct(mass_m.within_band)],
            ["mass-weighted total ratio (sum auto / sum manual)", _num(mass_m.total_ratio)],
            ["median signed RPD (bias)", _pct(mass_m.median_rpd)],
            ["detectable rows with no computed mass", f"{mass_m.n_unquantified}"],
        ],
    )
    lines += ["", f"Detectability threshold: manual mass >= {min_mass:g} ng.", ""]

    lines += ["## Recoveries", ""]
    lines += _table(
        ["metric", "value"],
        [
            ["joined (sample, homolog) pairs", f"{recovery_m.n}"],
            ["mean absolute error (recovery fraction)", _num(recovery_m.mae)],
            ["median signed error (bias)", _num(recovery_m.median_error)],
            ["max absolute error", _num(recovery_m.max_abs_error)],
            ["within +/-0.05", _pct(recovery_m.within_005)],
        ],
    )

    lines += ["", "## Response factors (RRFs, native analytes)", ""]
    lines += _table(
        ["metric", "value"],
        [
            ["native analytes compared", f"{rrf_m.n}"],
            ["median |RPD|", _pct(rrf_m.median_abs_rpd)],
            ["within +/-25%", _pct(rrf_m.within_band)],
            ["outside +/-25%", f"{rrf_m.n_off} ({rrf_m.n_high} read high, {rrf_m.n_low} read low)"],
            ["median |RPD| of the calibration areas behind them", _pct(rrf_m.median_abs_area_rpd)],
        ],
    )
    lines += [
        "",
        f"Detectable native mass rows outside +/-25% (see `bias_diagnosis.csv`): **{n_bias}**.",
        "",
        "## Verdict",
        "",
    ]
    lines += _table(
        ["level", "gate", "result"],
        [
            [
                "area", f"median \\|RPD\\| <= {_pct(GATE_MEDIAN_RPD)} and >= "
                f"{_pct(GATE_WITHIN_FRAC)} within +/-25%",
                "PASS" if area_pass else "FAIL",
            ],
            [
                "mass", f"median \\|RPD\\| <= {_pct(GATE_MEDIAN_RPD)} and >= "
                f"{_pct(GATE_WITHIN_FRAC)} within +/-25%",
                "PASS" if mass_pass else "FAIL",
            ],
            ["recovery", f"MAE <= {GATE_RECOVERY_MAE:g}", "PASS" if recovery_pass else "FAIL"],
        ],
    )
    lines += ["", f"**Overall: {'PASS' if passed else 'FAIL'}**", ""]
    return "\n".join(lines).rstrip() + "\n"


def generate_validation_report(
    result: RunResult,
    truth: RunTruth,
    out_dir: str | Path,
    run_name: str,
    *,
    min_area: float = DEFAULT_MIN_AREA,
    min_mass: float = DEFAULT_MIN_MASS,
    aliases: dict[str, str] | None = None,
) -> ValidationReport:
    """Write one run's validation report to ``out_dir`` and return its metrics.

    ``result`` is the production pipeline's output (``core.pipeline.run_analysis``) and
    ``truth`` is the analyst's ground truth for the same run (``report.truth.read_run_truth``).
    Writes ``summary.md``, four ``*_agreement.csv`` files, ``bias_diagnosis.csv`` and three
    PDF figures to ``out_dir`` (created if it does not exist), and returns the metrics that
    went into them plus the per-level and overall verdict.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    aliases = aliases or {}

    # ---- Areas: every reported file's peaks, joined against manual areas ------------------
    peaks_by_file = {
        fid: result.files[fid].peaks for fid in result.reported_ids if result.files[fid].usable
    }
    file_types = {fid: result.files[fid].file_type for fid in peaks_by_file}
    area_rows = join_areas(peaks_by_file, file_types, result.analytes, truth.areas, aliases)
    area_native = [r for r in area_rows if r.role == "native"]
    area_m = metrics(area_native, min_area)

    # ---- Masses -----------------------------------------------------------------------------
    mass_rows, _mass_truth_only, _mass_file_only = join_masses(result, truth.masses, aliases)
    mass_native = [r for r in mass_rows if r.role == "native"]
    mass_m = metrics(mass_native, min_mass)

    # ---- Recoveries -------------------------------------------------------------------------
    recovery_rows, _recovery_unmatched = join_recoveries(result, truth.recoveries, aliases)
    recovery_m = recovery_metrics(recovery_rows)

    # ---- RRFs -------------------------------------------------------------------------------
    matches = match_injections(result, truth.cal_areas)
    rrf_rows = join_rrfs(result, truth.cal_rrfs, truth.cal_areas, matches)
    rrf_m = rrf_metrics(rrf_rows)

    # ---- Bias diagnosis: detectable native mass rows outside the band ----------------------
    rrf_off = {r.label for r in rrf_rows if np.isfinite(r.rpd) and abs(r.rpd) > BAND}
    no_rrf_labels = {
        a.label for a in result.analytes
        if a.role == "native" and a.label not in result.calibration.rrf
    }
    detectable_native_mass = [r for r in mass_native if r.manual >= min_mass]
    bias_rows = [r for r in detectable_native_mass if abs(r.rpd) > BAND]

    paths: dict[str, Path] = {}

    paths["area_csv"] = out_dir / "area_agreement.csv"
    _write_csv(
        paths["area_csv"], _CSV_FIELDS["area"],
        [
            {"sample": r.sample, "file_type": r.file_type, "label": r.label,
             "homolog": r.homolog, "role": r.role, "manual": r.manual, "auto": r.auto,
             "rpd": r.rpd, "flags": r.flags, "method": r.method, "sn": r.sn,
             "rt_apex": r.rt_apex}
            for r in area_rows
        ],
    )
    paths["mass_csv"] = out_dir / "mass_agreement.csv"
    _write_csv(
        paths["mass_csv"], _CSV_FIELDS["mass"],
        [
            {"sample": r.sample, "file_type": r.file_type, "label": r.label,
             "homolog": r.homolog, "role": r.role, "manual": r.manual, "auto": r.auto,
             "rpd": r.rpd, "quantified": r.quantified, "flags": r.flags, "method": r.method,
             "rt_apex": r.rt_apex, "detected": r.detected}
            for r in mass_rows
        ],
    )
    paths["recovery_csv"] = out_dir / "recovery_agreement.csv"
    _write_csv(
        paths["recovery_csv"], _CSV_FIELDS["recovery"],
        [
            {"sample": r.sample, "homolog": r.homolog, "manual": r.manual, "auto": r.auto,
             "error": r.error}
            for r in recovery_rows
        ],
    )
    paths["rrf_csv"] = out_dir / "rrf_agreement.csv"
    _write_csv(
        paths["rrf_csv"], _CSV_FIELDS["rrf"],
        [
            {"label": r.label, "homolog": r.homolog, "ours": r.ours, "analyst": r.analyst,
             "rpd": r.rpd, "our_cal_area": r.our_cal_area,
             "analyst_cal_area": r.analyst_cal_area, "area_rpd": r.area_rpd, "flags": r.flags,
             "method": r.method}
            for r in rrf_rows
        ],
    )
    paths["bias_csv"] = out_dir / "bias_diagnosis.csv"
    _write_csv(
        paths["bias_csv"], _CSV_FIELDS["bias"],
        [
            {"sample": r.sample, "label": r.label, "homolog": r.homolog, "manual": r.manual,
             "auto": r.auto, "rpd": r.rpd, "flags": r.flags, "method": r.method,
             "clustered": "clustered" in r.flags.split(","),
             "likely_cause": likely_cause(r, min_mass, rrf_off, no_rrf_labels)}
            for r in bias_rows
        ],
    )

    paths["rpd_hist"] = out_dir / "rpd_hist.pdf"
    _fig_rpd_hist([r.rpd for r in area_native if r.manual >= min_area], paths["rpd_hist"])
    paths["regression"] = out_dir / "regression.pdf"
    _fig_regression(detectable_native_mass, paths["regression"])
    paths["homolog"] = out_dir / "homolog_agreement.pdf"
    _fig_homolog_agreement(detectable_native_mass, paths["homolog"])

    area_pass = (
        np.isfinite(area_m.median_abs_rpd)
        and area_m.median_abs_rpd <= GATE_MEDIAN_RPD
        and area_m.within_band >= GATE_WITHIN_FRAC
    )
    mass_pass = (
        np.isfinite(mass_m.median_abs_rpd)
        and mass_m.median_abs_rpd <= GATE_MEDIAN_RPD
        and mass_m.within_band >= GATE_WITHIN_FRAC
    )
    recovery_pass = np.isfinite(recovery_m.mae) and recovery_m.mae <= GATE_RECOVERY_MAE
    passed = area_pass and mass_pass and recovery_pass

    paths["summary"] = out_dir / "summary.md"
    paths["summary"].write_text(
        _render_summary(
            run_name, area_m, mass_m, recovery_m, rrf_m, len(bias_rows),
            area_pass, mass_pass, recovery_pass, passed, min_area, min_mass,
        ),
        encoding="utf-8",
    )

    return ValidationReport(
        run=run_name,
        out_dir=out_dir,
        area=area_m,
        mass=mass_m,
        recovery=recovery_m,
        rrf=rrf_m,
        n_bias_rows=len(bias_rows),
        area_pass=area_pass,
        mass_pass=mass_pass,
        recovery_pass=recovery_pass,
        passed=passed,
        paths=paths,
    )
