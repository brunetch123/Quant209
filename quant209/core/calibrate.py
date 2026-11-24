"""Derive this run's per-analyte retention-time anchors from its calibration injection.

Retention times drift between column builds, maintenance events and instruments, so the
canonical table's ``cal_rt`` values are a starting point rather than an address. The
calibration standard contains every congener. It is therefore the one injection where each
analyte's peak is guaranteed present. Each analyte is located near its canonical ``cal_rt``,
refined to the local apex, and that observed RT becomes the anchor for every other file in
the run. Every sample is then integrated against retention times measured on the same sequence.

Relative-response factors are not computed here. They need the calibration file's integrated
peaks, so they belong to quantification, which obtains them by calling
``integrate_run(cal_chroms, analytes, anchors)`` itself. That is the same call any sample
gets.

Two problems that the anchor step alone cannot solve are handled here as well.

**Retention-time alignment** (:func:`align_rt_map`). ``build_anchors`` searches ±0.10 min
around each canonical ``cal_rt``. That is the right window for column drift within one GC
method, and far too narrow when a run was acquired under a *different* oven program. On one
such run congener 1 elutes at 8.34 min against a canonical 19.56, and every analyte anchors
onto whatever noise happened to sit in its window. Anchoring is therefore preceded by what an
analyst does by eye when the method changes. The calibration injection's own peak pattern is
laid against the canonical retention-time sequence, and the straight line mapping one onto the
other is fitted channel by channel. The line is the correspondence-finder. Once it says which
observed peak is which congener, that peak's retention time is the search centre.

**A run-fit guard** (:func:`rt_fit_quality`). Both failures above are silent by construction:
integration always returns *a* number, so a run scored against the wrong retention times
produces a plausible-looking report rather than an error. The guard checks the calibration
injection against the table it was anchored with. On a run whose RT table fits, nearly
every analyte in the calibration standard is a resolved peak, so a low found-fraction there
is direct evidence that the table does not describe this run.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field, replace

import numpy as np

from .config import Analyte, analytes_by_transition
from .mzml import Chromatogram, get_transition
from .peaks import (
    IntegrationParams,
    global_noise,
    local_maxima,
    peak_prominences,
    refine_apexes,
    smooth,
)

# --------------------------------------------------------------------------------------
# Alignment tuning constants
#
# These are alignment-only and kept out of ``IntegrationParams``: that dataclass is
# calibrated against analyst manual integration and its defaults must not move, whereas
# these govern a pre-step that either fires (different GC method) or does not (same method,
# where the alignment is required to be exactly the identity).
# --------------------------------------------------------------------------------------

#: Fewest canonical analytes on a transition for its own linear map to be fitted. Below this
#: there is not enough of a retention-time *pattern* to identify against: a single labelled
#: standard on a quiet channel could be matched to any peak at all. Such channels inherit
#: instead (see :func:`align_rt_map`).
MIN_FIT_ANALYTES = 6
#
# Six, not three. The difference is large. A channel fits its own line by matching its
# analytes to detected peaks. With only three analytes there is nothing to contradict a wrong
# match, so the channel can place all three on the wrong peaks and report a perfect fit. On a
# run from a different instrument, where the whole chromatogram is re-timed and canonical
# 19.6 min elutes near 9.6, the three-analyte mono channel fitted a=0.95 and put congener 1
# at 12.4 min against the analyst's 9.6, while the twenty-analyte tri channel's line
# extrapolated into the same region gives 9.59. Measured over six such runs, requiring six
# analytes takes anchors more than 0.2 min from the analyst's own retention times from
# 88/170/61/63/65/117 (of 185) down to 18/19/0/14/15/15. Anything from 5 to 10 gives the same
# result. At 12 the well-populated channels start being excluded too and it degrades again.

#: Largest displacement (minutes, over the transition's own canonical RT span) at which a
#: fitted map is still reported as the identity and *not applied*.
#:
#: Fitting a line to peak apexes always returns something slightly off identity. On the run
#: the canonical table was extracted from, the fitted maps displace by at most 0.008 min,
#: which is noise in the apex positions rather than a method difference. Half the search
#: window is the threshold below which applying the map could only perturb an
#: already-correct search. Declaring identity there keeps a same-method run's anchors
#: bit-for-bit identical to running without alignment.
IDENTITY_TOL = 0.05

#: Slope prior on ``a``. A GC method change re-times the chromatogram. It does not compress or
#: stretch it by more than a factor of two. The extremes measured on this dataset are 0.791
#: (the mono channel of an alternate-method run) and 1.494 (its hepta channel), with the octa
#: channel's *unconstrained* fit reaching 1.514 before consensus pulls it back to 1.361.
#: Bounding the slope is what keeps the RANSAC search from "explaining" a channel with a
#: degenerate near-flat map that piles every congener onto one peak.
SLOPE_RANGE = (0.5, 2.0)

#: A candidate pair of canonical analytes must span at least this fraction of the transition's
#: canonical RT range before it is used to propose a slope. Two peaks close together give a
#: slope estimate dominated by their apex noise.
MIN_PAIR_SPAN_FRAC = 0.5

#: Detected peaks kept per canonical analyte on a transition, plus :data:`PEAK_POOL_PAD`. The
#: calibration standard's own congeners are not the only peaks on a channel (column bleed,
#: interferences, and neighbouring homologs all show up), so the pool has to be larger than the
#: analyte count. However, an unbounded pool is both slow and an invitation to spurious
#: matches.
PEAK_POOL_FACTOR = 3
PEAK_POOL_PAD = 10

#: Two transitions that share retention-time ground must agree, where they share it, to within
#: this many minutes. Retention shift is a property of the oven program, not of the MRM
#: channel, so neighbouring channels cannot disagree about it. This is a *gross mis-lock*
#: detector, not a fit-quality measure: on the reference runs correct neighbouring fits agree
#: to within 0.7 min even where the oven ramp knees and the straight-line model starts to
#: strain, while a mis-locked fit is out by nine minutes.
CONSENSUS_TOL = 1.0

#: Two transitions are compared only if their canonical RT ranges overlap, or fall short of
#: overlapping by less than this (consecutive homolog channels sometimes merely abut). Beyond
#: it, the comparison would be pure extrapolation of both lines and says nothing.
MAX_PROBE_GAP = 1.0

#: Radius (minutes) within which a channel with no fit of its own is re-anchored onto its own
#: dominant peak. Wide enough to absorb the borrowed prediction's error, narrow enough that
#: the peak found is the labelled standard rather than some other feature.
INHERIT_SEARCH = 1.0

#: Slack (minutes) allowed on the slope test that prunes the pooled retention-time table (see
#: :func:`_pool`). Two identifications a hair apart in canonical RT say nothing about the local
#: slope. Their apparent slope is dominated by where the two apexes were digitised, so the
#: test has to tolerate one search window's worth of play before it means anything.
CURVE_SLOPE_SLACK = 0.10

#: Cap on the identify -> pool -> re-identify loop, which settles in two or three passes.
IDENTIFY_PASSES = 6

#: Largest Q1 gap (Da) across which a transition may inherit a neighbour's map. Every labelled
#: channel sits 12-13 Da above its own native homolog channel. The widest native-to-native gap
#: in the canonical table is 34 Da (nona -> deca).
MAX_INHERIT_DA = 60.0

#: Run-fit guard thresholds, applied to the *calibration* injection. The calibration standard
#: contains every analyte, so on a run whose retention-time table fits, essentially all of them
#: are resolved peaks there (measured: 182 of 185 on a canonical-method run). These floors sit
#: far below that and far above what a mismatched table produces, which is near-total loss of
#: the labelled channels.
MIN_NATIVE_FOUND_FRAC = 0.70
MIN_LABELED_FOUND_FRAC = 0.50

#: Message the guard reports when a run fails it. Written for the analyst who will read it in
#: the GUI, not for the developer: the actionable fact is the GC method, not the code path.
FIT_FAIL_REASON = (
    "retention-time table does not fit this run's calibration injection. "
    "Check the GC method and the RT alignment"
)


# --------------------------------------------------------------------------------------
# Alignment model
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class RtMap:
    """One transition's retention-time map: ``observed_rt ~= a * canonical_rt + b``.

    A straight line is the whole model. It is what an analyst sketches when a run comes off a
    different oven program: "everything is about twelve minutes early, and the late half is
    squeezed". It is auditable at a glance: two numbers per channel, both printable in a
    report, both checkable by hand against any two peaks in the calibration injection.

    ``identity`` is not "a == 1 and b == 0". It is the verdict that this channel needs no
    alignment. ``a``/``b`` keep the *fitted* values either way so a reviewer can see what the
    fit actually said, while :meth:`__call__` returns the canonical RT untouched. That split is
    what lets a same-method run pass through this module without its anchors moving by a
    sample: the fit on such a run lands at ``a = 0.9999, b = +0.003``, which is apex noise, and
    applying it would be a change with no evidence behind it.
    """

    a: float
    b: float
    identity: bool  # displacement below IDENTITY_TOL: fitted, but not applied
    #: How this channel came by its line, one of :data:`MAP_SOURCES`:
    #:
    #: ``fitted``      : this channel's own peaks, possibly refitted under a neighbour's
    #:     consistency constraint (which is still a fit of this channel).
    #: ``constrained`` : no line satisfying that constraint existed, so the conflicting
    #:     neighbour's line is carried verbatim. ``a``/``b`` are the neighbour's, not this
    #:     channel's, and ``n_matched`` is 0.
    #: ``inherited``   : too few analytes here to identify a pattern. The nearest fitted
    #:     channel's line by Q1.
    #: ``global``      : no fitted channel near enough to inherit from. One line through every
    #:     matched pair in the run.
    #: ``none``        : nothing was fitted anywhere. The identity.
    source: str
    n_analytes: int  # canonical analytes on this transition
    n_matched: int  # of those, how many the fit placed on a detected peak

    def __call__(self, canonical_rt: float) -> float:
        """The retention time this run is expected to elute ``canonical_rt`` at."""
        if self.identity:
            return float(canonical_rt)
        return float(self.a * canonical_rt + self.b)

    def displacement(self, canonical_rt: float) -> float:
        """Minutes the fitted line moves ``canonical_rt`` by (reported even when identity)."""
        return float(self.a * canonical_rt + self.b - canonical_rt)


#: Every value :attr:`RtMap.source` may take. See the attribute for what each means.
MAP_SOURCES = frozenset({"fitted", "constrained", "inherited", "global", "none"})

#: Sources whose line came from *another* channel, so this one contributed no evidence of its
#: own. Only these are allowed the wide own-peak snap in stage 4 of :func:`align_rt_map`: a
#: channel that fitted its own peaks has already identified them at the search tolerance, and
#: re-snapping it at :data:`INHERIT_SEARCH` would let a crowded channel grab a neighbour's peak.
BORROWED_SOURCES = frozenset({"inherited", "global"})

#: The map used where nothing could be fitted or inherited: leave the canonical table alone.
IDENTITY_MAP = RtMap(a=1.0, b=0.0, identity=True, source="none", n_analytes=0, n_matched=0)


@dataclass(frozen=True)
class RtAlignment:
    """Where this run's calibration injection says each canonical retention time has moved to.

    ``maps`` is the reported model: one :class:`RtMap` per canonical transition key, as
    grouped by ``core.config.analytes_by_transition``.

    ``matched`` is the evidence: analyte label -> the retention time of the calibration peak
    the fit identified as that analyte. Where it has an entry, that observed apex is used as
    the analyte's expected RT in preference to the line, because the line's job was only to
    work out *which peak is which*. Fitting a straight line to a chromatogram that a different
    oven program has re-timed leaves a systematic misfit of a few tenths of a minute, enough
    to lose a congener in a crowded homolog window, while the peak the line identified is
    right where the chromatogram says it is.

    ``interpolated`` covers the analytes the calibration injection did not pin down
    individually. Retention time belongs to the oven program rather than to the mass channel,
    so the matched peaks from *every* channel together form one run-wide table of "canonical
    RT -> this run's RT". An unmatched analyte is placed by interpolating between the two
    identified peaks that bracket it. This is what carries the late-eluting congeners, where a
    single straight line per channel visibly strains: the oven ramp flattens around 58-60 min
    canonical, and the octa channel straddles that knee.

    Anything in neither dict falls back to its transition's line.

    Transitions whose map came out as the identity contribute *no* ``matched`` or
    ``interpolated`` entries: a run that needs no alignment must leave the anchor search
    exactly where it was.
    """

    maps: dict[tuple[float, float], RtMap] = field(default_factory=dict)
    matched: dict[str, float] = field(default_factory=dict)
    interpolated: dict[str, float] = field(default_factory=dict)

    @property
    def is_identity(self) -> bool:
        """True when no transition was displaced: this run matches the canonical GC method."""
        return all(m.identity for m in self.maps.values())

    def map_for(self, analyte: Analyte) -> RtMap:
        return self.maps.get((round(analyte.q1, 1), round(analyte.q3, 1)), IDENTITY_MAP)

    def expected_rt(self, analyte: Analyte) -> float:
        """This analyte's expected retention time in this run (minutes)."""
        rt = self.matched.get(analyte.label)
        if rt is None:
            rt = self.interpolated.get(analyte.label)
        if rt is not None:
            return float(rt)
        return self.map_for(analyte)(analyte.cal_rt)

    def source_of(self, analyte: Analyte) -> str:
        """How this analyte's expected RT was arrived at: matched | interpolated | line."""
        if analyte.label in self.matched:
            return "matched"
        if analyte.label in self.interpolated:
            return "interpolated"
        return "line"


# --------------------------------------------------------------------------------------
# Peak detection on a calibration trace
# --------------------------------------------------------------------------------------


def _detect_apexes(
    times: np.ndarray, intensities: np.ndarray, params: IntegrationParams, max_peaks: int
) -> np.ndarray:
    """Retention times of the most prominent peaks on one calibration trace, in RT order.

    Detection reuses the integrator's presence rule, topographic prominence on the smoothed
    trace against ``max(min_prominence, presence_sn * noise)``, so a peak the alignment is
    allowed to anchor to is a peak the integrator would also call present. The pool is then
    capped at ``max_peaks``, keeping the most prominent.
    """
    detect = smooth(intensities, params.smooth_window, params.smooth_poly)
    candidates = local_maxima(detect)
    if candidates.size == 0:
        return np.empty(0, dtype=float)
    proms = peak_prominences(detect, list(candidates))
    noise = global_noise(intensities.astype(float))
    keep = proms >= max(params.min_prominence, params.presence_sn * noise)
    candidates, proms = candidates[keep], proms[keep]
    if candidates.size > max_peaks:
        strongest = np.argsort(proms)[::-1][:max_peaks]
        candidates = candidates[strongest]
    return np.sort(times[candidates].astype(float))


# --------------------------------------------------------------------------------------
# Fitting one transition's line
# --------------------------------------------------------------------------------------


#: How many fewer matched peaks a candidate line may explain than the best and still be
#: preferred for agreeing better with the calibration standard's composition (see
#: :func:`_fit_line`). In a crowded channel a line shifted by one peak across part of the
#: channel puts each congener beside its neighbour's peak and matches about as many, and on
#: one run's second calibration two more than the right line (31 against 28 after
#: refinement), which a slack of one let stand.
SIZE_COUNT_SLACK = 2

#: How much better (correlation of log certified mass with log peak height over the matched
#: congeners) a candidate must agree with the standard's composition to displace the
#: best-matching line. Measured on a re-run whose hexa channel two of three calibrations
#: identified one peak off: the right line scored +0.36 and +0.37, the shifted one -0.11 and
#: -0.02, while the shifted one matched one peak more.
SIZE_CORR_MARGIN = 0.3

#: Fewest matched congeners for the composition to be read at all. On the small channels (the
#: di and octa channels carry 11) a correlation over a handful of masses that hardly differ is
#: noise: there it once preferred a line that re-timed the octa channel by 8.5 min and lost the
#: internal standard on it.
SIZE_MIN_MATCHED = 12


def _size_agreement(
    a: float, b: float, canonical: np.ndarray, observed: np.ndarray, tol: float,
    weights: np.ndarray, heights: np.ndarray,
) -> float:
    """Correlation of log certified mass with log matched-peak height under one line.

    A calibration standard holds known amounts: a co-elution group of three congeners makes a
    peak about three times its neighbours'. A correct identification puts the big groups on
    the big peaks. One shifted by a peak scrambles that. ``nan`` where too few congeners match
    or their masses hardly differ.
    """
    picks = _assign(a * canonical + b, observed, tol)
    sel = (picks >= 0) & np.isfinite(weights) & (weights > 0)
    if int(sel.sum()) < SIZE_MIN_MATCHED:
        return float("nan")
    lw = np.log(weights[sel])
    lh = np.log(np.maximum(heights[picks[sel]], 1.0))
    if float(np.std(lw)) < 0.1 or float(np.std(lh)) <= 0:
        return float("nan")
    return float(np.corrcoef(lw, lh)[0, 1])


def _nearest_distance(mapped: np.ndarray, observed: np.ndarray) -> np.ndarray:
    """Distance from each mapped RT to the closest observed apex (``mapped`` may be 2-D).

    An empty ``observed`` means no peak is anywhere near, which is ``inf`` rather than an
    IndexError. Callers already screen for it today. The guard is here regardless, so
    correctness does not depend on that always being true.
    """
    if observed.size == 0:
        return np.full(mapped.shape, np.inf)
    pos = np.searchsorted(observed, mapped)
    left = np.clip(pos - 1, 0, observed.size - 1)
    right = np.clip(pos, 0, observed.size - 1)
    return np.minimum(np.abs(mapped - observed[left]), np.abs(mapped - observed[right]))


def _assign(mapped: np.ndarray, observed: np.ndarray, tol: float) -> np.ndarray:
    """Order-preserving nearest-peak assignment: index into ``observed``, or -1 for no match.

    Both sequences are in retention-time order and an increasing map cannot reorder them, so
    the assignment is walked left to right and never allowed to go backwards. Without that
    constraint two neighbouring congeners in a crowded window can both claim the same peak,
    which reads as a confident match and is a swapped identification.
    """
    out = np.full(mapped.size, -1, dtype=int)
    cursor = 0
    for i, centre in enumerate(mapped):
        while cursor < observed.size and observed[cursor] < centre - tol:
            cursor += 1
        best, best_distance = -1, tol
        k = cursor
        while k < observed.size and observed[k] <= centre + tol:
            distance = abs(observed[k] - centre)
            if distance < best_distance:
                best, best_distance = k, distance
            k += 1
        if best >= 0:
            out[i] = best
            cursor = best + 1
    return out


def _score_candidates(
    a: np.ndarray, b: np.ndarray, canonical: np.ndarray, observed: np.ndarray, tol: float
) -> tuple[np.ndarray, np.ndarray]:
    """(matched-peak count, summed matched distance) for a batch of candidate lines."""
    mapped = a[:, None] * canonical[None, :] + b[:, None]
    distance = _nearest_distance(mapped, observed)
    hit = distance <= tol
    return hit.sum(axis=1), np.where(hit, distance, 0.0).sum(axis=1)


def _fit_line(
    canonical: np.ndarray,
    observed: np.ndarray,
    tol: float,
    probe: tuple[float, float] | None = None,
    weights: np.ndarray | None = None,
    heights: np.ndarray | None = None,
) -> tuple[float, float, int] | None:
    """RANSAC-lite fit of ``observed ~= a * canonical + b``, returning ``(a, b, n_matched)``.

    Every pair of canonical analytes that spans enough of the channel, crossed with every pair
    of detected apexes, proposes one candidate line. Each candidate is scored by how many of
    the channel's canonical retention times it lands within ``tol`` of a detected peak. That
    is the same question the analyst answers by eye ("does this line put every congener on a
    peak?"). ``tol`` is always ``IntegrationParams.rt_window``, threaded from the caller rather
    than fixed here: a match means "the anchor search would have found this peak", so widening
    that search has to widen matching with it or the alignment starts claiming correspondences
    the anchor step will not honour.
    Ties go to the candidate with the smaller total matched distance. The winner is then
    refined by least squares over its own matched pairs, twice, so the reported line is the
    best fit through the peaks it identified rather than the two peaks that proposed it.

    ``probe`` is an optional ``(retention_time, value)`` consistency constraint: candidates
    that do not put ``retention_time`` within :data:`CONSENSUS_TOL` of ``value`` are discarded
    before scoring. It is how a neighbouring channel's already-accepted map is imposed on this
    one (see :func:`align_rt_map`), and it is the difference between fitting a crowded channel
    correctly and locking onto a plausible-looking shifted solution.

    ``weights`` (each analyte's certified mass in the standard) and ``heights`` (each detected
    apex's height) break a near tie. Among candidates within :data:`SIZE_COUNT_SLACK` matches
    of the best, one whose matched peaks agree with the standard's composition by at least
    :data:`SIZE_CORR_MARGIN` more than the best's do is taken instead (:func:`_size_agreement`):
    the analyst's own check that the big co-elution group sits on the big peak. The tie is a
    tie between near neighbours: the candidate may move no position by more than
    :data:`CONSENSUS_TOL` from the best line's, and the best line must itself be readable
    (:data:`SIZE_MIN_MATCHED`). Without them the choice is by matches alone, as before.

    Returns ``None`` when the channel has too little to fit.
    """
    n, m = canonical.size, observed.size
    if n < MIN_FIT_ANALYTES or m < 2:
        return None
    span = float(canonical[-1] - canonical[0])
    if span <= 0:
        return None

    ci, cj = np.triu_indices(n, 1)
    wide = (canonical[cj] - canonical[ci]) >= MIN_PAIR_SPAN_FRAC * span
    ci, cj = ci[wide], cj[wide]
    pi, pj = np.triu_indices(m, 1)
    observed_gap = observed[pj] - observed[pi]
    # The mapped canonical range must land inside the acquired window (plus a minute of
    # slack for the segment edges): a line that puts half the congeners off the end of the
    # trace has not explained the channel, whatever it scores.
    window_lo, window_hi = observed[0] - 1.0, observed[-1] + 1.0
    slope_lo, slope_hi = SLOPE_RANGE

    best_count, best_distance, best_line = -1, np.inf, None
    candidates: list[tuple[int, float, float]] = []
    for i, j in zip(ci, cj):
        gap = float(canonical[j] - canonical[i])
        a = observed_gap / gap
        keep = (a >= slope_lo) & (a <= slope_hi)
        if not keep.any():
            continue
        a = a[keep]
        b = observed[pi[keep]] - a * canonical[i]
        inside = ((a * canonical[0] + b) >= window_lo) & ((a * canonical[-1] + b) <= window_hi)
        if probe is not None:
            probe_rt, probe_value = probe
            inside &= np.abs(a * probe_rt + b - probe_value) <= CONSENSUS_TOL
        if not inside.any():
            continue
        a, b = a[inside], b[inside]
        count, distance = _score_candidates(a, b, canonical, observed, tol)
        k = int(np.lexsort((distance, -count))[0])
        if count[k] > best_count or (count[k] == best_count and distance[k] < best_distance):
            best_count = int(count[k])
            best_distance = float(distance[k])
            best_line = (float(a[k]), float(b[k]))
        if weights is not None and heights is not None:
            near = count >= best_count - SIZE_COUNT_SLACK
            candidates.extend(zip(count[near].tolist(), a[near].tolist(), b[near].tolist()))

    if best_line is None:
        return None
    a, b = best_line
    own = (_size_agreement(a, b, canonical, observed, tol, weights, heights)
           if weights is not None and heights is not None and candidates else float("nan"))
    if np.isfinite(own):
        # The best-matching line, and its near ties, judged by the standard's composition.
        a0, b0 = a, b
        own_matched = int((_assign(a0 * canonical + b0, observed, tol) >= 0).sum())
        seen: set[tuple[float, float]] = set()
        best_size = own
        for count_c, a_c, b_c in candidates:
            key = (round(a_c, 3), round(b_c, 2))
            if count_c < best_count - SIZE_COUNT_SLACK or key in seen:
                continue
            seen.add(key)
            if float(np.max(np.abs((a_c - a0) * canonical + (b_c - b0)))) > CONSENSUS_TOL:
                continue  # not a neighbour of the best line: a re-timing, not a slipped peak
            if int((_assign(a_c * canonical + b_c, observed, tol) >= 0).sum()) < own_matched - SIZE_COUNT_SLACK:
                continue
            size = _size_agreement(a_c, b_c, canonical, observed, tol, weights, heights)
            if np.isfinite(size) and size >= own + SIZE_CORR_MARGIN and size > best_size:
                best_size, a, b = size, a_c, b_c
    for _ in range(2):
        picks = _assign(a * canonical + b, observed, tol)
        selected = picks >= 0
        if int(selected.sum()) < 2:
            break
        x, y = canonical[selected], observed[picks[selected]]
        design = np.vstack([x, np.ones_like(x)]).T
        (refit_a, refit_b), *_ = np.linalg.lstsq(design, y, rcond=None)
        if probe is not None and abs(refit_a * probe[0] + refit_b - probe[1]) > CONSENSUS_TOL:
            break  # the refinement drifted out of agreement with the neighbour: keep the winner
        a, b = float(refit_a), float(refit_b)
    n_matched = int((_assign(a * canonical + b, observed, tol) >= 0).sum())
    return a, b, n_matched


# --------------------------------------------------------------------------------------
# Alignment over the whole calibration injection
# --------------------------------------------------------------------------------------


def _probe_point(one: tuple[float, float], other: tuple[float, float]) -> float:
    """Retention time at which two canonical RT ranges are compared.

    The midpoint of their overlap, or (for two channels that merely abut, as consecutive
    homologs sometimes do) the midpoint of the gap between them. Both fall out of the same
    expression: ``max(lows)`` and ``min(highs)`` cross over when the ranges are disjoint.
    """
    return 0.5 * (max(one[0], other[0]) + min(one[1], other[1]))


def _make_map(
    a: float, b: float, canonical: np.ndarray, source: str, n_analytes: int, n_matched: int
) -> RtMap:
    """Wrap a fitted line, deciding whether it is a real displacement or apex noise."""
    displacement = np.abs((a - 1.0) * canonical + b).max() if canonical.size else abs(b)
    return RtMap(
        a=a,
        b=b,
        identity=bool(displacement < IDENTITY_TOL),
        source=source,
        n_analytes=n_analytes,
        n_matched=n_matched,
    )


def align_rt_map(
    cal_chroms: dict[tuple[float, float], Chromatogram],
    analytes: list[Analyte],
    params: IntegrationParams | None = None,
) -> RtAlignment:
    """Fit this run's retention-time map from its own calibration injection.

    The calibration standard contains every analyte, so every canonical peak must be somewhere
    on the corresponding channel of this file. That makes the alignment a pattern-matching
    problem rather than a search: lay the canonical retention-time sequence against the peaks
    the calibration injection actually shows and find the straight line that puts one on the
    other. It is the same reasoning an analyst applies by eye after a GC method change, and it
    reduces to two auditable numbers per channel.

    The fit proceeds in four stages.

    1. **Per-channel fit.** Every transition carrying at least :data:`MIN_FIT_ANALYTES`
       canonical analytes gets its own RANSAC-lite line (:func:`_fit_line`).

    2. **Cross-channel consensus.** Retention shift is a property of the oven program, so
       channels that share retention-time ground cannot disagree about it. Channels are
       accepted in decreasing order of the fraction of themselves their fit explained, so a
       confident channel constrains a doubtful one and never the reverse. Any channel that
       disagrees with an already-accepted neighbour by more than :data:`CONSENSUS_TOL` at
       their shared retention time is refitted under that constraint.

       Measured on an alternate-method run whose 32-analyte tetra channel is the case this
       stage exists for: its unconstrained line explains 22 of 32 analytes (69%) and sits nine
       minutes off at the probe, and the di channel (10 of 11, 91%) is the most confident
       neighbour that contradicts it. The refit moves the channel's intercept from -10.58 to
       -20.29. Most of that nine minutes never reaches the anchors either way, because stage 4
       re-identifies against the pooled table rather than the line. What the constraint is
       worth end to end is 20 anchors (tri through octa) moving by up to 0.63 min, and their
       worst error against the analyst's own retention times falling from 0.64 to 0.29 min.
       This is smaller than the nine-minute figure suggests on its own.

    3. **Fallback chain**, for channels with too few analytes to identify a pattern: every
       labelled-standard channel (one analyte) and the deca channel (one congener).

       a. inherit the line of the nearest *fitted* transition by Q1 (for a labelled channel
          that is always its own native homolog channel, 12-13 Da below).
       b. failing that (no fitted transition within :data:`MAX_INHERIT_DA`), a single global
          line least-squares-fitted through every matched pair from every fitted channel.
       c. failing that (nothing fitted anywhere), the identity: the canonical table, unmoved.

    4. **Identify, pool, place.** The lines say which observed peak is which congener. From
       there the chromatogram is the measurement. Matched peaks from every channel are pooled
       into one run-wide ``canonical RT -> this run's RT`` table, every channel is re-identified
       against that table, and analytes still unmatched are placed by interpolating it. A
       channel with no fit of its own then snaps its single analyte onto its own nearest peak
       within :data:`INHERIT_SEARCH`, since a borrowed line can identify a labelled standard
       but cannot time it. See :class:`RtAlignment` for why the straight line is not the last
       word: the oven ramp flattens around 58-60 min canonical, the octa channel straddles
       that knee, and no line describes it. However, its neighbours' identified peaks bracket
       every one of its congeners.

    A map whose largest displacement across its own channel is under :data:`IDENTITY_TOL` is
    reported but **not applied**, and contributes no matched or interpolated retention times:
    see :class:`RtMap`. On a run acquired under the canonical method every channel lands there,
    the alignment is the identity, and the anchors that come out of :func:`build_anchors` are
    unchanged to the last digit.
    """
    params = params or IntegrationParams()
    grouped = analytes_by_transition(analytes)

    # ---- stage 1: per-channel fit -----------------------------------------------------
    ordered: dict[tuple[float, float], list[Analyte]] = {}
    canonicals: dict[tuple[float, float], np.ndarray] = {}
    observations: dict[tuple[float, float], np.ndarray] = {}
    masses: dict[tuple[float, float], np.ndarray] = {}
    heights: dict[tuple[float, float], np.ndarray] = {}
    fits: dict[tuple[float, float], tuple[float, float, int]] = {}

    for transition, transition_analytes in grouped.items():
        anchored = sorted(
            (a for a in transition_analytes if not np.isnan(a.cal_rt)), key=lambda a: a.cal_rt
        )
        if not anchored:
            continue
        chrom = get_transition(cal_chroms, *transition)
        if chrom is None or chrom.times.size == 0:
            continue
        ordered[transition] = anchored
        canonical = np.array([a.cal_rt for a in anchored], dtype=float)
        canonicals[transition] = canonical
        pool = PEAK_POOL_FACTOR * len(anchored) + PEAK_POOL_PAD
        observed = _detect_apexes(chrom.times, chrom.intensities, params, pool)
        observations[transition] = observed
        # What the standard holds and how tall each detected peak is: the near-tie breaker.
        masses[transition] = np.array([a.cal_mass for a in anchored], dtype=float)
        heights[transition] = np.interp(observed, chrom.times, chrom.intensities.astype(float))
        fit = _fit_line(canonical, observed, params.rt_window,
                        weights=masses[transition], heights=heights[transition])
        if fit is not None:
            fits[transition] = fit

    maps: dict[tuple[float, float], RtMap] = {}
    if not fits:
        return RtAlignment(maps={t: IDENTITY_MAP for t in grouped})

    # ---- stage 2: cross-channel consensus ---------------------------------------------
    ranges = {t: (float(c[0]), float(c[-1])) for t, c in canonicals.items()}
    # Channels are accepted in decreasing order of how much of themselves their fit explained,
    # so a confident channel constrains a doubtful one and never the other way round. On
    # one alternate-method run the 32-analyte tetra channel is the doubtful one. Its
    # unconstrained line explains 22 of its 32 analytes (69%) and is wrong by nine minutes at
    # the probe, and the di channel, 10 of 11 (91%), is the most confident neighbour that
    # contradicts it.
    # Stage 4 recovers most of that nine minutes on its own. What the constraint is worth in
    # anchors is measured in ``align_rt_map``'s docstring.
    confidence = {t: (fits[t][2] / len(ordered[t]), fits[t][2]) for t in fits}
    accepted: dict[tuple[float, float], tuple[float, float, int, str]] = {}
    for transition in sorted(fits, key=lambda t: confidence[t], reverse=True):
        a, b, n_matched = fits[transition]
        conflict: tuple[tuple[float, int], float, float, tuple[float, float]] | None = None
        for other, (other_a, other_b, _n, _s) in accepted.items():
            if not _comparable(ranges[other], ranges[transition]):
                continue
            probe_rt = _probe_point(ranges[other], ranges[transition])
            target = other_a * probe_rt + other_b
            if abs(a * probe_rt + b - target) > CONSENSUS_TOL:
                if conflict is None or confidence[other] > conflict[0]:
                    conflict = (confidence[other], probe_rt, target, (other_a, other_b))
        source = "fitted"
        if conflict is not None:
            _, probe_rt, target, donor = conflict
            refit = _fit_line(
                canonicals[transition],
                observations[transition],
                params.rt_window,
                probe=(probe_rt, target),
                weights=masses[transition],
                heights=heights[transition],
            )
            if refit is not None:
                a, b, n_matched = refit
            else:
                # No line at all satisfies the constraint: carry the neighbour's, which at
                # least keeps this channel in the same retention-time frame as the rest of the
                # run. Recorded as "constrained" rather than "fitted": the numbers are the
                # neighbour's, and a reader of the report is entitled to know that.
                a, b, n_matched, source = donor[0], donor[1], 0, "constrained"
        accepted[transition] = (a, b, n_matched, source)

    for transition, (a, b, n_matched, source) in accepted.items():
        maps[transition] = _make_map(
            a, b, canonicals[transition], source, len(ordered[transition]), n_matched
        )

    fitted_keys = list(accepted)

    # ---- stage 3: fallback chain for channels that could not be fitted ----------------
    # Every labelled-standard channel (one analyte) and the deca channel land here.
    global_line = _global_line(accepted, canonicals, observations, params.rt_window)
    for transition, anchored in ordered.items():
        if transition in maps:
            continue
        canonical = canonicals[transition]
        donor = min(accepted, key=lambda t: abs(t[0] - transition[0]), default=None)
        if donor is not None and abs(donor[0] - transition[0]) <= MAX_INHERIT_DA:
            a, b = accepted[donor][0], accepted[donor][1]
            source = "inherited"
        elif global_line is not None:
            a, b = global_line
            source = "global"
        else:
            maps[transition] = IDENTITY_MAP
            continue
        maps[transition] = _make_map(a, b, canonical, source, len(anchored), 0)

    for transition in grouped:
        maps.setdefault(transition, IDENTITY_MAP)
    maps = _settle_identity(maps, fitted_keys)

    # ---- stage 4: identify peaks, pool them, and place the rest ------------------------
    # A channel reported as the identity is left strictly alone here: it needs no alignment,
    # and moving its search by a sample would change anchors that are already right.
    active = {t: m for t, m in maps.items() if not m.identity and t in ordered}
    matched: dict[str, float] = {}
    interpolated: dict[str, float] = {}
    if active:
        # The first pass identifies peaks off each channel's own line. Those matches, pooled
        # across every channel and pruned for smoothness, give a run-wide "canonical RT ->
        # this run's RT" table far better than any single line, because retention time belongs
        # to the oven program and not to the mass channel. Re-identifying against that table
        # picks up the analytes the lines missed and discards the ones they mis-claimed, which
        # improves the table in turn, so the pass is repeated until nothing moves.
        predicted = {t: np.array([maps[t](rt) for rt in canonicals[t]]) for t in active}
        previous: dict[str, float] = {}
        for _ in range(IDENTIFY_PASSES):
            matched, curve = _pool(
                _identify(active, ordered, observations, predicted, params.rt_window), analytes
            )
            if curve is None or matched == previous:
                break
            previous = matched
            predicted = {t: _predict(canonicals[t], curve, maps[t]) for t in active}
        # A channel with no fit of its own (every labelled-standard channel, and deca) is one
        # analyte on a quiet trace. Nothing on it can be identified by pattern, so the borrowed
        # prediction does the identifying and the channel's own chromatogram supplies the
        # retention time: the nearest peak within INHERIT_SEARCH, if there is one.
        for transition, rt_map in active.items():
            if rt_map.source not in BORROWED_SOURCES:
                continue
            observed = observations[transition]
            for analyte, prediction in zip(ordered[transition], predicted[transition]):
                if analyte.label in matched or observed.size == 0:
                    continue
                nearest = int(np.argmin(np.abs(observed - prediction)))
                if abs(observed[nearest] - prediction) <= INHERIT_SEARCH:
                    matched[analyte.label] = float(observed[nearest])
        for transition in active:
            for analyte, prediction in zip(ordered[transition], predicted[transition]):
                if analyte.label not in matched:
                    interpolated[analyte.label] = float(prediction)

    return RtAlignment(maps=maps, matched=matched, interpolated=interpolated)


def _settle_identity(
    maps: dict[tuple[float, float], RtMap], fitted: list[tuple[float, float]]
) -> dict[tuple[float, float], RtMap]:
    """Take the identity verdict at *run* level, from the channels that were actually fitted.

    Only fitted channels carry evidence about whether this run needs aligning. A borrowed line
    is evidence about its donor. And a borrowed line is judged over the *borrower's* retention
    range, which can sit well past the donor's: a slope a thousandth off identity displaces a
    channel eluting at 22 min by nothing and one eluting at 60 min by half a search window. So
    a same-method run could otherwise have every fitted channel declare identity while a
    labelled channel that borrowed one of those lines went active, and the bit-identical-anchors
    guarantee would hold for 173 analytes and quietly fail for the other 12.

    The fitted ``a``/``b``/``source`` are preserved either way: a reviewer asking "what did the
    fit actually say on this run?" is asking about the numbers, not the decision taken on them.
    """
    if not fitted or not all(maps[t].identity for t in fitted):
        return maps
    return {t: replace(m, identity=True) for t, m in maps.items()}


def _comparable(one: tuple[float, float], other: tuple[float, float]) -> bool:
    """Whether two canonical RT ranges share enough ground to be checked against each other."""
    return min(one[1], other[1]) - max(one[0], other[0]) >= -MAX_PROBE_GAP


def _predict(
    canonical: np.ndarray, curve: tuple[np.ndarray, np.ndarray], rt_map: RtMap
) -> np.ndarray:
    """Place canonical RTs on the run-wide table, falling back to the channel's line outside it."""
    placed = np.interp(canonical, curve[0], curve[1], left=np.nan, right=np.nan)
    line = rt_map.a * canonical + rt_map.b
    return np.where(np.isnan(placed), line, placed)


def _identify(
    active: dict[tuple[float, float], RtMap],
    ordered: dict[tuple[float, float], list[Analyte]],
    observations: dict[tuple[float, float], np.ndarray],
    predicted: dict[tuple[float, float], np.ndarray],
    tol: float,
) -> dict[str, float]:
    """Analyte label -> the calibration peak this prediction identifies as that analyte."""
    matched: dict[str, float] = {}
    for transition in active:
        observed = observations[transition]
        if observed.size == 0:
            continue
        picks = _assign(np.asarray(predicted[transition], dtype=float), observed, tol)
        for analyte, pick in zip(ordered[transition], picks):
            if pick >= 0:
                matched[analyte.label] = float(observed[pick])
    return matched


def _pool(
    matched: dict[str, float], analytes: list[Analyte]
) -> tuple[dict[str, float], tuple[np.ndarray, np.ndarray] | None]:
    """The run-wide ``canonical RT -> observed RT`` table, from every identified peak.

    Pooling across channels is not a convenience: it is the physical claim that the oven
    program, not the mass channel, sets retention time. It is also what makes the late,
    sparsely-populated channels tractable: the octa channel straddles the point where the
    ramp flattens and no straight line describes it, but its neighbours' identified peaks
    bracket every one of its congeners.

    Not every identification survives into the table. A retention-time map is monotone and its
    local slope stays inside :data:`SLOPE_RANGE`: peaks cannot swap order, and the run cannot
    stand still for two minutes and then leap. The rule keeps the **largest set of identifications
    that together form a physically possible map** and discard the rest: a longest-chain
    selection over the pooled points, where a step from one identification to the next is
    admissible when both coordinates increase and the implied slope is in range (with
    :data:`CURVE_SLOPE_SLACK` of play, since two peaks a hair apart say nothing about slope).

    Selecting by consistency rather than by residual is what makes the sparse late end work.
    The last four minutes of the run carry three congeners and a surrogate. A smoothness test
    there has to extrapolate over the gap and throws away the very points that define it,
    whereas the chain test only asks whether they fit a possible map, which they do. The peak
    the octa channel's overshooting line claims a minute past congener 205 does not.

    A rejected point is dropped from ``matched`` as well, not merely from the table: the
    verdict is that the peak is not that analyte, and an analyte with no credible
    identification belongs in the interpolated population. Returns ``(kept matches, table)``.
    """
    triples = sorted(
        (a.cal_rt, matched[a.label], a.label) for a in analytes if a.label in matched
    )
    if len(triples) < 2:
        return {}, None
    xs = np.array([t[0] for t in triples], dtype=float)
    ys = np.array([t[1] for t in triples], dtype=float)
    chain = _longest_consistent(xs, ys)
    if chain.size < 2:
        return {}, None
    kept = {triples[i][2]: float(ys[i]) for i in chain}
    return kept, (xs[chain], ys[chain])


def _longest_consistent(xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    """Indices of the longest chain of points forming an admissible retention-time map.

    A standard longest-chain dynamic program over points already sorted by canonical RT: a
    step ``i -> j`` is admissible when both coordinates strictly increase and the implied slope
    lies in :data:`SLOPE_RANGE`, allowing :data:`CURVE_SLOPE_SLACK` minutes of play so that two
    nearly co-eluting analytes are judged on order alone. O(n^2) on the ~180 identified peaks
    of a run costs nothing in practice, and it gives the *maximal* consistent explanation
    rather than a greedy one: the difference between dropping one wrong identification and
    dropping the three correct ones that disagree with it.
    """
    n = xs.size
    length = np.ones(n, dtype=int)
    previous = np.full(n, -1, dtype=int)
    slope_lo, slope_hi = SLOPE_RANGE
    for j in range(1, n):
        # ``dx > 0`` is strict, so two analytes sharing a canonical ``cal_rt`` exactly can
        # never both be kept: the chain takes whichever gives the longer explanation. No
        # such pair exists in the shipped table (all 185 ``cal_rt`` values are distinct), but
        # a future table that reports two channels' analytes at one retention time would
        # silently lose one identification per collision rather than mis-order the map.
        dx = xs[j] - xs[:j]
        dy = ys[j] - ys[:j]
        admissible = (
            (dx > 0)
            & (dy > 0)
            & (dy >= slope_lo * dx - CURVE_SLOPE_SLACK)
            & (dy <= slope_hi * dx + CURVE_SLOPE_SLACK)
        )
        if not admissible.any():
            continue
        reachable = np.where(admissible, length[:j], 0)
        best = int(np.argmax(reachable))
        if reachable[best] > 0:
            length[j] = reachable[best] + 1
            previous[j] = best
    chain: list[int] = []
    index = int(np.argmax(length))
    while index >= 0:
        chain.append(index)
        index = int(previous[index])
    return np.array(chain[::-1], dtype=int)


def _global_line(
    accepted: dict[tuple[float, float], tuple[float, float, int, str]],
    canonicals: dict[tuple[float, float], np.ndarray],
    observations: dict[tuple[float, float], np.ndarray],
    tol: float,
) -> tuple[float, float] | None:
    """One line least-squares-fitted through every matched pair on every fitted channel.

    The last resort before giving up and leaving the canonical table alone. A single line
    cannot describe the whole run well: the shift is markedly retention-time dependent, and
    on an alternate-method run a global fit misses by up to 4 min against the analyst's table.
    But it is the difference between "roughly where the run put things" and "the canonical
    method's retention times", and it is only reached by a channel with no fitted neighbour
    within :data:`MAX_INHERIT_DA`.
    """
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    for transition, (a, b, _n, _s) in accepted.items():
        canonical = canonicals[transition]
        observed = observations[transition]
        picks = _assign(a * canonical + b, observed, tol)
        selected = picks >= 0
        if selected.any():
            xs.append(canonical[selected])
            ys.append(observed[picks[selected]])
    if not xs:
        return None
    x, y = np.concatenate(xs), np.concatenate(ys)
    if x.size < 2:
        return None
    design = np.vstack([x, np.ones_like(x)]).T
    (a, b), *_ = np.linalg.lstsq(design, y, rcond=None)
    return float(a), float(b)


# --------------------------------------------------------------------------------------
# Anchors
# --------------------------------------------------------------------------------------




# --------------------------------------------------------------------------------------
# The run's retention-time warp, read off the labelled standards
# --------------------------------------------------------------------------------------

#: How much taller a labelled channel's biggest peak must be than the next one before it is
#: accepted as that standard's peak. The labelled standards are spiked at high concentration
#: and each sits alone on its own transition, so the real peak towers over whatever else the
#: channel picks up: measured across this archive it is typically 20-500x the runner-up.
#: Requiring a clear margin is what makes these anchors trustworthy without any RT prior at
#: all, which is the entire point of using them.
STANDARD_DOMINANCE = 5.0

#: Fewest labelled standards that must be identified this way before their warp is used. They
#: span mono to deca, so a handful already pins the whole retention range. Below this there
#: is not enough of a curve to trust and the per-channel machinery is left to do the work.
MIN_STANDARD_ANCHORS = 6


def standards_warp(
    cal_chroms: dict[tuple[float, float], Chromatogram],
    analytes: list[Analyte],
    params: IntegrationParams,
) -> list[tuple[float, float]]:
    """``[(canonical RT, observed RT)]`` for every labelled standard found unambiguously.

    This is the one measurement in a calibration injection that needs no retention-time prior.
    A native transition carries twenty or thirty congeners in a row of similar peaks, so
    matching them to expectations is a puzzle with many near-solutions: slide the
    expectations along by a few peaks and everything still lines up. That is exactly how a
    whole run ends up anchored minutes off with nothing to notice it. A labelled channel
    carries *one* analyte and, because that standard is spiked at high concentration, one
    towering peak. There is nothing to confuse it with.

    The twelve standards span mono to deca, so the pairs they give trace the run's whole
    retention-time warp, including the curvature that defeats a straight line. On a 2025 run
    from a different instrument the local slope runs 1.05 early to 1.48 late, and a line
    fitted anywhere lands ~4 min out at the other end.

    Returns the pairs sorted by canonical RT, keeping only a strictly increasing subset: a
    warp that doubles back is not a chromatogram, it is a mis-identification.
    """
    pairs: list[tuple[float, float]] = []
    for transition, group in analytes_by_transition(analytes).items():
        if any(a.role == "native" for a in group) or len(group) != 1:
            continue
        analyte = group[0]
        if not np.isfinite(analyte.cal_rt):
            continue
        chrom = get_transition(cal_chroms, *transition)
        if chrom is None or chrom.times.size == 0:
            continue
        detected = detect_cal_peaks(chrom, params)
        if not detected:
            continue
        heights = [float(chrom.intensities[i]) for i in detected]
        order = sorted(range(len(detected)), key=lambda i: heights[i], reverse=True)
        best = order[0]
        runner_up = heights[order[1]] if len(order) > 1 else 0.0
        if runner_up > 0 and heights[best] < STANDARD_DOMINANCE * runner_up:
            continue  # two comparable peaks: this channel cannot say which is the standard
        pairs.append((float(analyte.cal_rt), float(chrom.times[detected[best]])))

    pairs.sort()
    increasing: list[tuple[float, float]] = []
    for canonical, observed in pairs:
        if not increasing or observed > increasing[-1][1]:
            increasing.append((canonical, observed))
    return increasing


def warp_rt(pairs: list[tuple[float, float]], canonical: float) -> float:
    """Map one canonical retention time through the standards' warp.

    Piecewise linear between the anchors, and continued at the end segments' own slope
    outside them: the standards bracket almost the whole range, so extrapolation is short.
    """
    if not pairs:
        return canonical
    if len(pairs) == 1:
        return canonical + (pairs[0][1] - pairs[0][0])
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    if canonical <= xs[0]:
        slope = (ys[1] - ys[0]) / (xs[1] - xs[0]) if xs[1] > xs[0] else 1.0
        return ys[0] + slope * (canonical - xs[0])
    if canonical >= xs[-1]:
        slope = (ys[-1] - ys[-2]) / (xs[-1] - xs[-2]) if xs[-1] > xs[-2] else 1.0
        return ys[-1] + slope * (canonical - xs[-1])
    return float(np.interp(canonical, xs, ys))



def seed_from_standards(
    cal_chroms: dict[tuple[float, float], Chromatogram],
    analytes: list[Analyte],
    params: IntegrationParams | None = None,
) -> list[Analyte]:
    """The analyte table with its retention times moved onto this run's own clock.

    The labelled standards say where this run's chromatogram sits without needing any prior
    (:func:`standards_warp`). Re-expressing the canonical table through their warp gives the
    per-channel alignment a starting point that is already in the right neighbourhood, which
    is what stops a densely populated channel settling into a shifted match: with twenty
    similar peaks in a row, a channel starting tens of minutes out has many local optima and
    no way to tell them apart. On a 2025 run it chose one that was ~4 min late across
    mono..tetra while hexa..octa were exact.

    Seeding rather than replacing: the standards give a dozen anchors, while the native
    channels between them give well over a hundred, so the alignment still does the fine work.
    It is simply no longer free to do it in the wrong place. A run whose standards sit
    where the table already says they do is returned unchanged.
    """
    params = params or IntegrationParams()
    pairs = standards_warp(cal_chroms, analytes, params)
    if len(pairs) < MIN_STANDARD_ANCHORS:
        return analytes
    if max(abs(observed - canonical) for canonical, observed in pairs) <= IDENTITY_TOL:
        return analytes
    return [
        replace(analyte, cal_rt=warp_rt(pairs, float(analyte.cal_rt)))
        if np.isfinite(analyte.cal_rt)
        else analyte
        for analyte in analytes
    ]


# --------------------------------------------------------------------------------------
# Elution-order anchor assignment
# --------------------------------------------------------------------------------------

#: How far an analyte's predicted retention time may sit from a detected calibration peak
#: and still be matched to it (minutes). This is far larger than
#: ``IntegrationParams.rt_window``: the window is how far the *integrator* will hunt once an
#: anchor is fixed, whereas this is how wrong the *prediction* is allowed to be before the
#: order constraint stops being able to rescue it. Measured need on this archive: the
#: hepta/octa region of one 2024 run predicts up to 1.3 min off, and a tolerance of
#: 1.2 min is what recovers the analyst's own assignment there (PCB204 at 54.53 against the
#: analyst's 54.49, where the unassigned anchor sat at 53.85 and integrated noise).
ORDER_MATCH_TOL = 1.2

#: Score charged for leaving an expected analyte unmatched, in units of a perfect match.
#: Above ~0.25 the search prefers explaining every analyte, which is the wanted behaviour
#: on a calibration injection: the standard contains all 209 congeners, so an unmatched
#: analyte means a missed peak far more often than a genuinely absent one.
ORDER_GAP_EXPECTED = 0.35

#: Score charged for skipping a detected peak. Small: extra peaks are ordinary (isotope
#: bleed from a neighbouring channel, column contamination), and refusing to skip them
#: would force real analytes onto them.
ORDER_GAP_DETECTED = 0.02

#: Below this many detected peaks per expected analyte the assignment is not attempted:
#: too little of the window is visible for order alone to identify anything, and the
#: aligned prediction is the better estimate.
ORDER_MIN_COVERAGE = 0.5

#: Fewest expected analytes on a transition before order-based assignment is used at all.
#:
#: Order is the whole basis of the method: it identifies a peak by *which* peak it is in a
#: sequence. A channel carrying one or two analytes has no sequence, so the assignment
#: degenerates to "snap to the nearest peak within ``ORDER_MATCH_TOL``", with a tolerance
#: far wider than a peak, chosen for windows where order does the constraining.
#: The labelled-standard channels are exactly this case: one 13C surrogate on its own
#: transition, and nothing to order it against. Left unguarded it cost a whole run's hepta
#: recoveries: 13C 180 in one reference injection snapped 0.93 min onto a neighbouring peak,
#: and every recovery measured against that reference came out ~80x instead of ~1.
#: Below this count the aligned prediction stands and the apex search refines it directly.
ORDER_MIN_ANALYTES = 3

_ORDER_NEG = -1e18


def assign_in_elution_order(
    predicted: list[float],
    detected: list[float],
    tol: float = ORDER_MATCH_TOL,
    gap_expected: float = ORDER_GAP_EXPECTED,
    gap_detected: float = ORDER_GAP_DETECTED,
) -> list[int | None]:
    """Match expected analytes to detected calibration peaks, preserving elution order.

    Within one MRM transition the congeners' elution order is fixed by their chemistry: a
    change of column or oven programme moves the peaks and can change their *spacing*, but
    it cannot make PCB 204 elute before PCB 201. So the expected analytes and the detected
    peaks are two ordered sequences to be matched monotonically. That is exactly what an
    analyst does by eye when a run's retention times have shifted.

    This is strictly stronger than moving each analyte to its own nearest peak. Nearest-peak
    matching has no memory, so in a crowded window a prediction that is off by more than the
    peak spacing hands every analyte its neighbour's peak, and nothing detects the collision.
    Here two analytes can never claim the same peak, and the order they come out in is a
    constraint rather than an outcome.

    Needleman-Wunsch over the two sequences: a match scores ``1 - |Δ|/tol`` (and is
    disallowed beyond ``tol``), an unmatched analyte costs ``gap_expected``, a skipped peak
    costs ``gap_detected``. Because a match is worth at most 1 and a gap costs a fraction of
    that, the optimum explains as many analytes as the order allows and uses closeness only
    to choose between equally complete explanations.

    Returns, for each predicted analyte, the index of its detected peak, or None.
    """
    n, m = len(predicted), len(detected)
    if n == 0 or m == 0:
        return [None] * n
    score = np.full((n + 1, m + 1), _ORDER_NEG, dtype=float)
    # 0 = match, 1 = leave this analyte unmatched, 2 = skip this detected peak
    back = np.zeros((n + 1, m + 1), dtype=np.int8)
    score[0, 0] = 0.0
    for i in range(n + 1):
        for j in range(m + 1):
            if i == 0 and j == 0:
                continue
            best, move = _ORDER_NEG, 1
            if i > 0 and score[i - 1, j] > _ORDER_NEG:
                cand = score[i - 1, j] - gap_expected
                if cand > best:
                    best, move = cand, 1
            if j > 0 and score[i, j - 1] > _ORDER_NEG:
                cand = score[i, j - 1] - gap_detected
                if cand > best:
                    best, move = cand, 2
            if i > 0 and j > 0 and score[i - 1, j - 1] > _ORDER_NEG:
                delta = abs(predicted[i - 1] - detected[j - 1])
                if delta <= tol:
                    cand = score[i - 1, j - 1] + (1.0 - delta / tol)
                    if cand > best:
                        best, move = cand, 0
            score[i, j] = best
            back[i, j] = move
    out: list[int | None] = [None] * n
    i, j = n, m
    while i > 0 or j > 0:
        move = back[i, j]
        if move == 0 and i > 0 and j > 0:
            out[i - 1] = j - 1
            i, j = i - 1, j - 1
        elif move == 1 and i > 0:
            i -= 1
        else:
            j -= 1
    return out


def detect_cal_peaks(chrom: Chromatogram, params: IntegrationParams) -> list[int]:
    """Indices of the calibration trace's resolved peaks, by the integrator's own rule.

    The same prominence-against-noise test the integrator uses to decide presence, so an
    anchor can never be assigned to something the integrator would then refuse to call a
    peak.
    """
    detect = smooth(chrom.intensities, params.smooth_window, params.smooth_poly)
    maxima = list(local_maxima(detect))
    if not maxima:
        return []
    proms = peak_prominences(detect, maxima)
    floor = max(params.min_prominence, params.presence_sn * global_noise(chrom.intensities))
    return [i for i, prom in zip(maxima, proms) if prom >= floor]


#: A calibration peak whose apex lies this close (minutes) to the apex of a peak on the channel
#: of the homolog one or two chlorines higher can be that congener's in-source copy (a ghost
#: peak, see :mod:`core.ghosts`): the copy's apex falls within one acquisition cycle of its
#: parent's (0.30 to 0.43 s on the archive's methods, 0.005 to 0.007 min), and this allows two.
COPY_RT_TOL = 0.012
#: ... and it is taken for the copy when its prominence is no more than this many times the copy
#: the parent is expected to leave (the largest ratio of that homolog and step in
#: ``ghost_fractions.csv``, times PCB 204's default fraction). A congener's own calibration peak
#: is several times larger than any copy, since the standard holds every congener at comparable
#: amounts. Measured on a run under the shorter oven programme: PCB 170's copy on the hexa channel
#: had a prominence of 928 against PCB 169's own 44,331, and PCB 208's on the hepta channel 11,500
#: against PCB 189's 84,200.
COPY_SIZE = 2.0
#: An anchor this far (minutes) from the drift interpolated from the other transitions' anchors
#: sends its transition back for a second assignment (see :func:`_drift_second_look`). On runs
#: under the shorter oven programme the misplaced anchors sat 0.3 to 1.7 min from that drift;
#: the agreeing ones sit within a few hundredths of a minute.
DRIFT_OUTLIER = 0.15


def copy_peaks(
    cal_chroms: dict[tuple[float, float], Chromatogram],
    analytes: list[Analyte],
    params: IntegrationParams,
) -> dict[tuple[float, float], set[int]]:
    """Detected calibration peaks that are copies of a more chlorinated congener's peak.

    A PCB losing HCl or Cl2 in the ion source leaves a copy of its peak on the channel one or two
    homologs down, at its own retention time. In a calibration standard every congener is present,
    so these copies sit among the analytes' own peaks, and an analyte whose predicted retention
    time is off can be matched to a copy instead of its own peak: on runs acquired under a shorter
    oven programme, PCB 169 was matched to the copy of PCB 170 and PCB 189 to that of PCB 208.

    Returns, per native transition, the indices (as :func:`detect_cal_peaks` reports them) of the
    peaks that coincide with a peak on a parent homolog's channel and are no larger than a copy.
    """
    from .config import load_ghost_fractions
    from .ghosts import DEFAULT_IS_FRACTION, HOMOLOGS, STEPS

    by_transition = analytes_by_transition(analytes)
    homolog_of: dict[tuple[float, float], str] = {}
    for transition, members in by_transition.items():
        natives = [a.homolog for a in members if a.role in ("native", "IS_204") and a.homolog]
        if natives:
            homolog_of[transition] = max(set(natives), key=natives.count)
    transition_of = {h: t for t, h in homolog_of.items()}
    ratios = load_ghost_fractions()
    largest: dict[tuple[str, int], float] = {}
    for a in analytes:
        if a.role == "native" and a.homolog:
            for step, r in ratios.get(a.label, {}).items():
                largest[(a.homolog, step)] = max(largest.get((a.homolog, step), 0.0), r)
        elif a.role == "IS_204":
            # The ratios are relative to PCB 204's own copies, so 204 counts as 1.
            for step in (1, 2):
                largest[(a.homolog, step)] = max(largest.get((a.homolog, step), 0.0), 1.0)

    peaks: dict[tuple[float, float], tuple[list[int], np.ndarray, np.ndarray]] = {}
    for transition in homolog_of:
        chrom = get_transition(cal_chroms, *transition)
        if chrom is None or chrom.times.size == 0:
            continue
        idx = detect_cal_peaks(chrom, params)
        if not idx:
            continue
        detect = smooth(chrom.intensities, params.smooth_window, params.smooth_poly)
        peaks[transition] = (idx, chrom.times[idx], np.asarray(peak_prominences(detect, idx)))

    out: dict[tuple[float, float], set[int]] = {}
    for transition, homolog in homolog_of.items():
        if transition not in peaks:
            continue
        idx, times, proms = peaks[transition]
        rank = HOMOLOGS.index(homolog)
        for step in (1, 2):
            if rank + step >= len(HOMOLOGS):
                continue
            parent = HOMOLOGS[rank + step]
            parent_tr = transition_of.get(parent)
            if parent_tr not in peaks or step not in STEPS.get(parent, ()):
                continue
            ratio = largest.get((parent, step))
            if not ratio:
                continue
            limit = COPY_SIZE * ratio * DEFAULT_IS_FRACTION[step]
            _pidx, ptimes, pproms = peaks[parent_tr]
            for i, t, prom in zip(idx, times, proms):
                near = np.abs(ptimes - t) <= COPY_RT_TOL
                if near.any() and prom <= limit * float(pproms[near].max()):
                    out.setdefault(transition, set()).add(int(i))
    return out


def build_anchors(
    cal_chroms: dict[tuple[float, float], Chromatogram],
    analytes: list[Analyte],
    params: IntegrationParams | None = None,
    alignment: RtAlignment | None = None,
) -> dict[str, float]:
    """Map analyte label -> observed apex RT (minutes) measured from the cal injection.

    Each transition's analytes are laid out in canonical ``cal_rt`` order and refined through
    :func:`core.peaks.refine_apexes`, the same apex definition the integrator applies to every
    other file, shared so an anchor cannot be measured under one rule and consumed under
    another. The apex is recorded unconditionally: no presence test gates an anchor. A weak
    or missing cal peak still yields the best available RT estimate, and the *integration*
    step is where presence is judged and flagged. Gating anchors here would instead drop the
    analyte from every file in the run.

    Where the cal file has no chromatogram for a transition at all (a truncated conversion,
    or a method that never acquired it), that transition's analytes fall back to their
    canonical ``cal_rt`` so they are still integrated, at the table's retention time, with
    whatever confidence flags the integrator then assigns.

    ``alignment`` (from :func:`align_rt_map`) moves each analyte's *expected* retention time
    before that search, for a run acquired under a different GC method. It changes nothing
    else: the apex is still measured by ``refine_apexes`` on the calibration trace, so the
    alignment only decides where to look, never what is reported. Passing ``None`` (the
    default) leaves each analyte's expected retention time unchanged, and an alignment that
    came out as the identity produces byte-identical anchors to passing ``None``, by
    construction.

    Two guards keep the assignment on the analytes' own peaks. The copies that more
    chlorinated congeners leave on a channel are never offered to it (:func:`copy_peaks`), and
    a transition whose anchors disagree with the drift the other transitions show is assigned a
    second time from that drift (:func:`_drift_second_look`). On runs under a shorter oven
    programme, a transition's straight line placed its last analytes up to 1.7 min early, and
    PCB 169 and PCB 189 were anchored on the copies of PCB 170 and PCB 208 that sat there.

    Returns a plain ``label -> RT`` mapping. The calibration file's ``PeakResult`` records
    that RRF derivation needs come from integrating the cal file with these anchors,
    not from this function.
    """
    params = params or IntegrationParams()
    expected_rt = (
        (lambda a: float(a.cal_rt)) if alignment is None else alignment.expected_rt
    )
    # Copies of more chlorinated congeners are never offered to the elution-order assignment
    # (see copy_peaks): a copy is not an analyte's calibration peak.
    copies = copy_peaks(cal_chroms, analytes, params)
    anchors: dict[str, float] = {}
    state: dict = {}
    for transition, transition_analytes in analytes_by_transition(analytes).items():
        anchored = [a for a in transition_analytes if not np.isnan(a.cal_rt)]
        # Sorted on the *aligned* retention times, since those are what the search and the
        # territory midpoints are laid out along. An increasing map cannot reorder a
        # transition, so on an unaligned (or identity) run this is the canonical order.
        anchored.sort(key=expected_rt)
        if not anchored:
            continue
        chrom = get_transition(cal_chroms, *transition)
        if chrom is None or chrom.times.size == 0:
            # No cal trace for this transition: keep the canonical RTs rather than dropping
            # these analytes out of the run entirely.
            for analyte in anchored:
                anchors[analyte.label] = expected_rt(analyte)
            continue
        expected = [expected_rt(a) for a in anchored]

        # Stage 1: identify this transition's peaks and match them to the expected
        # congeners in elution order. This is what rescues a window whose *spacing*
        # changed, not just its offset. See ``assign_in_elution_order``. Analytes the
        # assignment cannot place keep their aligned prediction and are refined below,
        # exactly as before, so this stage can add information but never remove it.
        excluded = copies.get(transition, set())
        detected_idx = [i for i in detect_cal_peaks(chrom, params) if i not in excluded]
        assigned: list[int | None] = [None] * len(anchored)
        if (
            len(anchored) >= ORDER_MIN_ANALYTES
            and len(detected_idx) >= ORDER_MIN_COVERAGE * len(anchored)
        ):
            detected_rt = [float(chrom.times[i]) for i in detected_idx]
            assigned = assign_in_elution_order(expected, detected_rt)

        # Stage 2: the apex search, from wherever stage 1 left each analyte expecting to
        # find it. refine_apexes is the integrator's own apex definition, shared so that an
        # anchor cannot be measured under one rule and consumed under another.
        searched = [
            float(chrom.times[detected_idx[j]]) if j is not None else expected[i]
            for i, j in enumerate(assigned)
        ]
        _detect, _bounds, apexes = refine_apexes(chrom.times, chrom.intensities, searched, params)
        for analyte, apex in zip(anchored, apexes):
            anchors[analyte.label] = float(chrom.times[apex])
        state[transition] = (anchored, chrom, detected_idx, assigned)

    # Stage 3: a second assignment for any transition whose anchors disagree with the drift the
    # other transitions show (see _drift_second_look).
    if state:
        anchors = _drift_second_look(anchors, state, params)
    return anchors


def _drift_curve(points: list[tuple[float, float]]):
    """Piecewise-linear shift (observed minus canonical) against canonical retention time."""
    if len(points) < 2:
        return None
    cs = np.array([p[0] for p in points], dtype=float)
    ss = np.array([p[1] for p in points], dtype=float)
    order = np.argsort(cs)
    cs, ss = cs[order], ss[order]
    uc, inv = np.unique(cs, return_inverse=True)
    us = np.bincount(inv, weights=ss) / np.bincount(inv)
    return lambda c: float(np.interp(c, uc, us))


def _drift_second_look(
    anchors: dict[str, float],
    state: dict,
    params: IntegrationParams,
) -> dict[str, float]:
    """Second look at a transition whose anchors disagree with the drift of the other transitions.

    A straight line per transition places the last analytes of a long transition too early on a
    run whose oven programme compresses the end of the chromatogram differently, and the
    elution-order assignment then hands each of them an earlier peak. The run's other transitions
    still show where the chromatography went at that retention time. Each anchor is compared with
    the drift interpolated from the anchors of every other transition. In a transition with
    disagreeing anchors, those analytes (and any left unassigned) are predicted from the drift of
    all agreeing anchors, the others keep their own anchors as their prediction, and the
    transition is assigned again. A transition whose anchors all agree is left exactly as it was,
    and a run whose transitions all agree comes out unchanged.

    ``state`` holds, per transition, ``(analytes, chromatogram, detected peak indices, first
    assignment)`` as :func:`build_anchors` left them.
    """
    points: dict = {}
    for tr, (anchored, _chrom, _idx, assigned) in state.items():
        lst = [(float(a.cal_rt), anchors[a.label] - float(a.cal_rt), a.label)
               for a, j in zip(anchored, assigned) if j is not None and a.role == "native"]
        if lst:
            points[tr] = lst
    outliers: dict = {}
    for tr, lst in points.items():
        f = _drift_curve([(c, s) for t2, l2 in points.items() if t2 != tr for (c, s, _l) in l2])
        if f is None:
            continue
        for c, s, lab in lst:
            if abs(s - f(c)) > DRIFT_OUTLIER:
                outliers.setdefault(tr, set()).add(lab)
    if not outliers:
        return anchors
    curve = _drift_curve([(c, s) for tr, lst in points.items() for (c, s, lab) in lst
                          if lab not in outliers.get(tr, set())])
    if curve is None:
        return anchors
    out = dict(anchors)
    for tr, labels in outliers.items():
        anchored, chrom, detected_idx, assigned = state[tr]
        pred = [float(a.cal_rt) + curve(float(a.cal_rt)) if (a.label in labels or j is None) else anchors[a.label]
                for a, j in zip(anchored, assigned)]
        detected_rt = [float(chrom.times[i]) for i in detected_idx]
        again = assign_in_elution_order(pred, detected_rt)
        searched = [float(chrom.times[detected_idx[j]]) if j is not None else pred[i] for i, j in enumerate(again)]
        _detect, _bounds, apexes = refine_apexes(chrom.times, chrom.intensities, searched, params)
        for a, apex in zip(anchored, apexes):
            out[a.label] = float(chrom.times[apex])
    return out



# --------------------------------------------------------------------------------------
# Cross-calibration reconciliation
# --------------------------------------------------------------------------------------

#: How far one calibration's anchor for an analyte may sit from what the run's other
#: calibrations say, once the smooth drift between injections has been removed, before it is
#: treated as a mis-identification rather than as chromatography (minutes).
#:
#: Retention times drift smoothly across a sequence, and that drift is very nearly the same
#: for every analyte in a neighbourhood, so once it is subtracted, one calibration's anchor
#: for a given congener should agree with the others to well within a peak width. A
#: *disagreement* of half a minute is not drift. It is one injection's peak identification
#: having slipped by a position. Measured on one re-run: its third
#: calibration placed congener 66 at 30.33 min where the other three agree on 30.76-30.84
#: (a spread of 0.08), and every file anchored to it then reported the huge 61+70+74+76 peak
#: as congener 66, a 35 million counts·s error in a single row.
RECONCILE_TOL = 0.15

#: Fewest calibration injections needed to reconcile at all: with two, a disagreement has no
#: majority to be judged against and both readings stand.
RECONCILE_MIN_CALS = 3


def reconcile_anchors(
    anchors_by_cal: dict[str, dict[str, float]],
    order: list[str] | None = None,
    tol: float = RECONCILE_TOL,
) -> tuple[dict[str, dict[str, float]], list[str]]:
    """Repair anchors that one calibration injection places where the others do not.

    Each injection is anchored independently, so each can independently mis-identify a peak.
    The elution-order assignment constrains the order *within* an injection, but nothing
    stops a whole window sliding by one position if that injection's prediction was already
    a peak out. Comparing the injections with each other is what catches that, because a
    slip shows up as one reading disagreeing with its peers by far more than the drift
    between them.

    For each analyte the per-injection anchors are compared after removing that injection's
    *global* drift (the median shift of all its anchors from the reference injection's), so
    genuine drift is not mistaken for disagreement. A reading whose residual exceeds ``tol``
    while the majority agree is replaced by the majority's estimate, carried to that
    injection by its own drift.

    Returns the repaired mapping and a list of human-readable repairs, which the pipeline
    records: a silently corrected anchor would be exactly the kind of invisible decision this
    project refuses to make.
    """
    names = order or list(anchors_by_cal)
    names = [n for n in names if n in anchors_by_cal]
    if len(names) < RECONCILE_MIN_CALS:
        return anchors_by_cal, []

    reference = names[0]
    ref_anchors = anchors_by_cal[reference]

    # One number per injection: how far it has drifted from the reference overall.
    drift: dict[str, float] = {}
    for name in names:
        deltas = [
            anchors_by_cal[name][label] - ref_anchors[label]
            for label in anchors_by_cal[name]
            if label in ref_anchors
        ]
        drift[name] = float(np.median(deltas)) if deltas else 0.0

    repaired = {name: dict(anchors_by_cal[name]) for name in names}
    notes: list[str] = []
    labels = set().union(*(set(anchors_by_cal[n]) for n in names))
    for label in sorted(labels):
        # Each injection's opinion of this analyte, expressed on the reference's clock.
        opinions = {
            name: anchors_by_cal[name][label] - drift[name]
            for name in names
            if label in anchors_by_cal[name]
        }
        if len(opinions) < RECONCILE_MIN_CALS:
            continue
        consensus = float(np.median(list(opinions.values())))
        outliers = {n: v for n, v in opinions.items() if abs(v - consensus) > tol}
        # A majority must agree before a minority is overruled. If half the injections
        # disagree there is no consensus to repair towards.
        if not outliers or len(outliers) * 2 >= len(opinions):
            continue
        for name, value in outliers.items():
            was = repaired[name][label]
            repaired[name][label] = consensus + drift[name]
            notes.append(
                f"{name}: anchor for {label} moved {was:.3f} -> {repaired[name][label]:.3f} min "
                f"({value - consensus:+.3f} from the other calibrations' consensus)"
            )
    return repaired, notes


# --------------------------------------------------------------------------------------
# Run-fit guard
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class RunFitReport:
    """Whether this run's retention-time table actually describes its calibration injection.

    ``fit_ok`` is the verdict. ``reason`` is empty when it passes and carries an
    analyst-readable explanation when it does not.
    """

    fit_ok: bool
    native_found_frac: float
    labeled_found_frac: float
    n_native: int
    n_native_found: int
    n_labeled: int
    n_labeled_found: int
    reason: str = ""

    @property
    def summary(self) -> str:
        """One line for a report or a log."""
        return (
            f"native {self.n_native_found}/{self.n_native} "
            f"({100 * self.native_found_frac:.0f}%), labeled "
            f"{self.n_labeled_found}/{self.n_labeled} ({100 * self.labeled_found_frac:.0f}%)"
            + (f": {self.reason}" if self.reason else "")
        )


def rt_fit_quality(
    cal_chroms: dict[tuple[float, float], Chromatogram],
    analytes: list[Analyte],
    anchors: dict[str, float],
    params: IntegrationParams | None = None,
) -> RunFitReport:
    """Ask the calibration injection whether the retention-time table it was anchored with fits.

    Pointing the engine at the wrong retention times is silent: integration always returns a
    number, the report fills in, and nothing distinguishes "these areas are wrong" from "these
    areas are small". The calibration injection is the one file where that can be checked
    without any ground truth, because it is the one file where every analyte is *known* to be
    present: it is the standard they were all spiked into. It is integrated with the
    anchors in question and count how many analytes come back as resolved peaks.

    On a run whose table fits, essentially all of them do (measured: 182 of 185 on a
    canonical-method run). On the same engine pointed at a table from a different GC method,
    the labelled-standard channels collapse almost completely: they are single peaks on quiet
    channels, so when the search window is in the wrong place there is nothing there at all.
    Hence two floors rather than one: :data:`MIN_NATIVE_FOUND_FRAC` on the natives, which
    degrade gradually, and the stricter-in-effect :data:`MIN_LABELED_FOUND_FRAC` on the
    labelled standards, which are the sharp indicator.

    Either floor breached sets ``fit_ok=False`` with :data:`FIT_FAIL_REASON`. The guard is
    advisory to the caller: the pipeline surfaces it rather than refusing to run, since an
    analyst looking at a flagged run is exactly who can tell a method change from a bad
    injection.
    """
    from .integrate import integrate_run  # local: integrate imports nothing from calibrate

    params = params or IntegrationParams()
    with warnings.catch_warnings():
        # Analytes with no anchor or no trace are already counted as not-found below. The
        # integrator's warning about them would be a second report of the same fact.
        warnings.simplefilter("ignore", RuntimeWarning)
        results = integrate_run(cal_chroms, analytes, anchors, params)
    found = {r.label for r in results if r.detected}

    native = [a for a in analytes if a.role == "native"]
    labeled = [a for a in analytes if a.role != "native"]
    n_native_found = sum(1 for a in native if a.label in found)
    n_labeled_found = sum(1 for a in labeled if a.label in found)
    native_frac = n_native_found / len(native) if native else 1.0
    labeled_frac = n_labeled_found / len(labeled) if labeled else 1.0

    fit_ok = native_frac >= MIN_NATIVE_FOUND_FRAC and labeled_frac >= MIN_LABELED_FOUND_FRAC
    return RunFitReport(
        fit_ok=fit_ok,
        native_found_frac=native_frac,
        labeled_found_frac=labeled_frac,
        n_native=len(native),
        n_native_found=n_native_found,
        n_labeled=len(labeled),
        n_labeled_found=n_labeled_found,
        reason="" if fit_ok else FIT_FAIL_REASON,
    )
