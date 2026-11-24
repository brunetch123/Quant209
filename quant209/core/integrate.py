"""Orchestrate RT-anchored integration across a whole file, transition by transition.

Each MRM transition carries an ordered list of expected analytes (single congeners,
inseparable co-elution groups, or labelled standards). For every analyte, the apex is
refined inside its own RT territory, the peak is bounded by the valley shared with a
resolved neighbour (perpendicular drop) or by the local baseline minimum when it stands
alone, and the region is integrated above a straight baseline joining the two boundary
points.

The per-run RT anchors come from the calibration injection (see ``core.calibrate``), so
retention times track each column build instead of being hardcoded.

The decision tree below is calibrated against analyst manual integration, which it agrees
with at a median |RPD| of ~2.5%. The presence test (prominence above *both* an absolute
floor and an S/N-scaled floor) and the valley-to-valley boundary choice were each set by
comparison with analyst areas, and changing either moves results away from that
calibration.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from .config import Analyte, analytes_by_transition, load_aroclor_abundance
from .deconvolve import PairFit, RiderFit, deconvolve_rider, fit_pair
from .mzml import Chromatogram, get_transition
from .peaks import (
    IntegrationParams,
    PeakResult,
    _local_min,
    global_noise,
    hill_climb,
    integrate_above_baseline,
    integrate_region,
    find_shoulder,
    integrate_skim,
    peak_prominences,
    refine_apexes,
    valley_index,
)

#: Every boundary-selection method a ``PeakResult`` may report.
#:
#: ``valley_drop``: at least one neighbouring analyte was a resolved peak, so this peak's
#:     boundary on that side is the shared valley (a perpendicular drop between the two).
#: ``tangent_skim``: this analyte is a rider on a larger neighbour's flank, found by its
#:     curvature rather than by a local maximum. Its boundaries are the inflection points
#:     that bracket the bulge, and its baseline is the chord between them, which on a flank is
#:     the parent peak's local tangent.
#: ``rt_window``: no resolved neighbour bounded this peak. An isolated peak is bounded by
#:     the local minima within a half-width of its apex, a peak beside an *absent* neighbour
#:     extends through it to its own baseline edge, and a shoulder/absent analyte is bounded
#:     valley-to-valley inside its own midpoint territory.
#: ``deconvolved``: this analyte is a rider whose area is the fitted area of its own
#:     component in a shape fit of the envelope it shares with its neighbours
#:     (``core.deconvolve``). The area is therefore NOT the trace above the persisted baseline
#:     between the persisted boundaries. Those describe where the rider was looked for.
KNOWN_METHODS = frozenset({"valley_drop", "rt_window", "tangent_skim", "deconvolved"})

#: Every flag a ``PeakResult`` may carry, listed in the order a peak carries them. The
#: first three are confidence verdicts. ``flags`` is the detection verdict and
#: ``area``/``height``/``sn`` are conditional on it (see ``PeakResult.detected``). The last
#: five record how the number was produced.
#:
#: ``absent``: the presence test rejected this analyte. Its apex prominence did not clear
#:     both the absolute floor and the noise-scaled floor. The reported area is then a bounded
#:     estimate of whatever sits in the analyte's RT territory, NOT a peak measurement, and
#:     its ``sn``/``height`` are not meaningful (the apex may lie outside the integrated
#:     region). This is the flag downstream code must gate on.
#: ``low_sn``: the peak's baseline-corrected apex height is below ``min_sn`` times the
#:     trace noise. Real but near the detection limit, so the area is imprecise. Only
#:     informative on a detected peak. An absent analyte can report a huge ``sn`` and so
#:     will usually NOT carry this flag.
#: ``narrow``: fewer than 3 samples between the integration boundaries. A trapezoidal area
#:     over one or two intervals is dominated by where the boundaries landed rather than by
#:     the peak shape, so the number should not be trusted as a quantity.
#: ``shoulder``: the prominence test rejected this analyte but its neighbour's flank
#:     bulges at its retention time, so it was integrated as a rider (``tangent_skim``).
#:     Not a confidence verdict. A skimmed area depends on where the inflections fell, which
#:     is a softer judgement than a valley, so a reviewer should check it.
#: ``clustered``: this area was measured above a baseline spanning a *multi-peak cluster*
#:     rather than the chord between this peak's own two boundaries (see
#:     :func:`_cluster_spans`). Not a confidence verdict: it says which baseline
#:     ``baseline_start``/``baseline_end`` describe, which is what anyone recomputing the area
#:     by hand needs to know. Boundaries are unaffected, so ``method`` still reports how they
#:     were chosen.
#: ``forced``: neither the presence test nor the shoulder test found this analyte in this
#:     file, and it was integrated anyway because a rule names it
#:     (``IntegrationParams.forced_riders`` / ``forced_present``). Every other peak is
#:     evidence read off this trace. This one is a rule applied to a label, so it is the
#:     first thing a reviewer should look at, and the flag is what lets them find it.
#:     ``method`` says which rule: ``deconvolved`` for a fitted rider, ``valley_drop`` or
#:     ``rt_window`` for an analyte bounded as if it were present.
#: ``displaced``: this analyte's apex lies outside its own retention-time territory. Its
#:     peak drifted into a minor neighbour's territory and was given back to it on the
#:     strength of the published Aroclor profiles (see :func:`_reclaim_displaced`). The
#:     peak is real and measured from the trace. Which congener it belongs to is the
#:     judgement a reviewer should check.
#: ``apportioned``: this analyte and a neighbour of comparable abundance had merged into one
#:     envelope, and its area is its share by a rigid two-peak fit (see
#:     :func:`_reclaim_displaced`). Both members of the pair carry the flag.
#: ``unresolved``: this run's calibration injection could not separate this analyte from its
#:     neighbour (see :func:`unresolved_pairs`), and in this file the two share one envelope.
#:     Its area is its component of a two-peak fit at the spacing the calibration measured
#:     between them (``core.deconvolve.fit_pair``), and its record spans the envelope. Both
#:     members of the pair carry the flag.
KNOWN_FLAGS = frozenset(
    {"absent", "low_sn", "narrow", "clustered", "shoulder", "forced", "displaced", "apportioned",
     "unresolved"}
)

#: A congener must reach this weight fraction in at least one common Aroclor before a peak
#: displaced across its boundary can be given to it. Only a major component makes a peak
#: large enough to broaden and drift that far.
MAJOR_FRACTION = 0.01


@dataclass(frozen=True)
class _Forced:
    """One analyte a named rule integrated: where, by which rule, and the fit if it was one."""

    start: int
    apex: int
    end: int
    method: str
    fit: RiderFit | None  # None when the analyte was bounded as if present


def _halfwidth_points(times: np.ndarray, max_halfwidth: float) -> int:
    """Convert a max half-width in minutes to a sample count for this trace.

    One count for the whole trace, from its median sampling interval: what a rider's envelope
    (``_fit_rider``) and the displaced-peak drift cap use. A peak's own boundaries measure the
    minutes on the trace's times instead (``_reach_limits``), because an MRM time segment can
    sample a transition faster or slower than the rest of its trace.
    """
    if times.size < 2:
        return 5
    dt = float(np.median(np.diff(times)))
    return max(2, int(round(max_halfwidth / dt))) if dt > 0 else 5


def _reach_limits(times: np.ndarray, apex: int, max_halfwidth: float) -> tuple[int, int]:
    """The furthest indices an unbounded side may reach: ``max_halfwidth`` minutes either side of ``apex``.

    Measured on this trace's own times, never as one sample count for the whole trace. A
    transition's MRM time segments can sample it at different rates: the hexa channel is sampled
    every 0.549 s up to 47-51 min and every 0.305 s after, so a count taken from its median
    interval (27 samples) reached only 0.137 min where 156+157 and PCB 169 elute, and those peaks
    were cut on their own tails. Measured on the 142-run dataset before the change: 169's
    calibration area a median 0.67 of the analysts' on the newer oven program and 0.78 on the
    older one, 156+157 in newer-program samples 0.80. Where a segment is sampled more slowly than
    the median the count reached too far instead (0.32 min on the hepta channel's early segment).
    At least one sample either side, so a boundary search is never empty.
    """
    t0 = float(times[apex])
    lo = int(np.searchsorted(times, t0 - max_halfwidth, side="left"))
    hi = int(np.searchsorted(times, t0 + max_halfwidth, side="right")) - 1
    return max(0, min(lo, apex - 1)), min(times.size - 1, max(hi, apex + 1))


def _resolve_boundaries(
    detect: np.ndarray,
    bounds: list[tuple[int, int]],
    apexes: list[int],
    present: list[bool],
    i: int,
    times: np.ndarray,
    max_halfwidth: float,
) -> tuple[int, int, str]:
    """Choose analyte ``i``'s integration boundaries and the method string describing them.

    Boundaries and method are decided together, from one evaluation of each neighbour
    predicate, so a reported ``method`` can never describe boundaries other than the ones
    returned beside it.

    Three cases:

    * The analyte failed the presence test. It is a shoulder on a neighbour, or nothing at
      all. It gets its own midpoint territory, split valley-to-valley about the territory
      centre. (A cluster baseline anchored to the floor over-assigns a big peak's underlying
      area to a small shoulder. The analyst skims shoulders instead.)
    * A neighbour is itself a resolved peak. The boundary on that side is the valley the two
      share, a perpendicular drop that gives each peak the area on its own side.
    * The neighbour is absent or there is none. This peak extends *through* that neighbour's
      territory to its own baseline edge, capped at ``max_halfwidth`` minutes from the apex
      (``_reach_limits``, on the trace's own ``times``).
    """
    apex = apexes[i]

    if not present[i]:
        lo_b, hi_b = bounds[i]
        mid = (lo_b + hi_b) // 2
        return _local_min(detect, lo_b, mid, lo_b), _local_min(detect, mid, hi_b, hi_b), "rt_window"

    left_resolved = i > 0 and present[i - 1]
    right_resolved = i < len(apexes) - 1 and present[i + 1]
    reach_lo, reach_hi = _reach_limits(times, apex, max_halfwidth)

    left_limit = valley_index(detect, apexes[i - 1], apex) if left_resolved else reach_lo
    i_start = _local_min(detect, left_limit, apex, left_limit)

    right_limit = valley_index(detect, apex, apexes[i + 1]) if right_resolved else reach_hi
    i_end = _local_min(detect, apex, right_limit, right_limit)

    # One resolved neighbour is enough for this peak to have been split at a shared valley.
    has_resolved_neighbour = left_resolved or right_resolved
    return i_start, i_end, "valley_drop" if has_resolved_neighbour else "rt_window"


def _cluster_spans(
    detect: np.ndarray,
    resolved: list[tuple[int, int, str]],
    apexes: list[int],
    present: list[bool],
    overlap_frac: float,
) -> list[tuple[int, int]]:
    """Per-analyte ``(base_left, base_right)``: the span its integration baseline rests on.

    Perpendicular-drop boundaries tile a crowded window edge to edge, and drawing each slice's
    baseline as the chord between *its own* two boundaries puts that baseline on top of the
    elevated valleys either side of it. Where several congeners sit on one unresolved bed,
    such as a calibration injection where the standard contains all 209 at comparable
    concentration, the summed area of the slices then misses most of the bed. On one
    calibration injection the five-congener penta window at 39.0-39.9 min summed to 271,316
    counts·s against the analyst's 768,411, a **65%** loss across the window. Its worst single
    member, '87+125', integrated to 29,839 against 298,886, a **90%** loss. With the cluster
    baseline those become 754,427 (-1.8%) and 351,012 (+17.4%): the rule reduces the error on
    the window and slightly overshoots on the single worst congener.

    So consecutive analytes are grouped into a cluster, and every member of a cluster is
    integrated above **one** baseline spanning the whole cluster
    (:func:`core.peaks.integrate_above_baseline`) while keeping its own valley-drop
    boundaries. Two neighbours join a cluster when both of them:

    * are **present**. An analyte the presence test rejected breaks the run: drawing a
      floor-anchored baseline under a *shoulder* hands it the dominant neighbour's underlying
      area, whereas the analyst skims shoulders. That case is untouched. Only runs of
      genuinely resolved peaks are merged.
    * do not **return to base** between their apexes. The valley between them stands more
      than ``overlap_frac`` of the smaller apex's height above the lower of the pair's outer
      boundaries. Two peaks whose valley reaches the chromatographic floor are separately
      resolved, need no shared baseline, and are left alone whatever ``overlap_frac`` is. The
      parameter is only the tolerance on "did it reach the floor".

    The threshold is therefore a soft one, but what is stable across it differs from what is
    not. The **merge set is not**: on one run's ``Cal01`` it shrinks monotonically from 37
    merged pairs at 0.02 to 31 / 22 / 20 at 0.05 / 0.10 / 0.15 (each a strict subset of the
    last), and 42 analytes' areas differ between the two endpoints, by as much as -22%. What is
    stable is the **outcome**: the count of response factors more than 25% from the analyst's
    is flat across 0.02-0.15 on five of the six (injection, workbook) pairs of the two
    reference runs, and moves by one on the sixth. The plateau is in the metric, not in the
    mechanism, and 0.05 sits on the conservative, fewer-merges, half of a range whose
    endpoints both pass. Regenerate with
    ``python tools/bench_quant.py --sweep overlap_frac=10,0.02,0.05,0.10,0.15``. The table is
    in ``docs/METHODS.md`` section 2.4.

    Clusters are short in practice: on ``Cal01`` at the shipped 0.05 there are 26 of them, 21
    of two members and 5 of three, and none longer anywhere in either reference run's
    calibration injections. That bounds the span of a single baseline chord, which is what
    keeps the straight-line assumption under a cluster safe. A chord drawn across a long,
    sloping bed would under-cut one end. Nothing here enforces a maximum length, because
    nothing in the data has needed one. The two guards actually relied on are the ``np.clip``
    in :func:`core.peaks.integrate_above_baseline`, which prevents a mis-drawn chord from
    contributing negative area, and the two regression tests that pin clustering to never
    lower an area and never cross an absent analyte. A future dataset with longer clusters
    should re-check that before trusting this.
    """
    n = len(apexes)
    group = list(range(n))
    for i in range(n - 1):
        if not (present[i] and present[i + 1]):
            continue
        valley = valley_index(detect, apexes[i], apexes[i + 1])
        # The pair's own outer boundaries are the best local estimate of where the
        # chromatographic floor is under it. They are where the drop found the trace lowest.
        base = min(float(detect[resolved[i][0]]), float(detect[resolved[i + 1][1]]))
        lower_apex = min(float(detect[apexes[i]]), float(detect[apexes[i + 1]]))
        if lower_apex <= base:
            continue
        if (float(detect[valley]) - base) / (lower_apex - base) > overlap_frac:
            group[i + 1] = group[i]
    # group[i] is always at most i, so one left-to-right pass fully resolves each chain. Merging is
    # transitive by construction here: three mutually unresolved peaks become one cluster with
    # one baseline chord spanning all three. That is the intended reading of "the trace never
    # returned to base", and it is also the only place cluster length can grow without bound.
    # See the docstring for why the observed maximum of three members keeps the straight-chord
    # assumption safe on this data.
    for i in range(1, n):
        group[i] = group[group[i]]

    spans: dict[int, tuple[int, int]] = {}
    for i in range(n):
        lo, hi = resolved[i][0], resolved[i][1]
        current = spans.get(group[i])
        if current is not None:
            lo, hi = min(current[0], lo), max(current[1], hi)
        spans[group[i]] = (lo, hi)
    return [spans[group[i]] for i in range(n)]


def _reclaim_displaced(
    times: np.ndarray,
    detect: np.ndarray,
    integ: np.ndarray,
    ordered: list[Analyte],
    expected: list[float],
    bounds: list[tuple[int, int]],
    apexes: list[int],
    present: list[bool],
    params: IntegrationParams,
) -> tuple[dict[int, int], dict[int, tuple[int, PairFit]]]:
    """Give a dominant peak back to the analyte whose apex search it drifted out of.

    In a concentrated sample a dominant peak broadens, and its top can drift across the
    boundary it shares with a neighbour. The analyte's apex search, confined to its own
    territory, then stops on the peak's rising flank. That point has no prominence, so the
    presence test rejects the largest peak in the chromatogram, and the neighbour's search
    finds the top and claims the whole peak. The trace alone cannot say whose peak it is.
    An analyst decides by knowing which congener is the big one in these mixtures, and so
    does this rule (``core.config.load_aroclor_abundance``).

    Analyte ``i`` takes the top neighbour ``j`` claimed only when all of these hold:

    * ``i`` failed the presence test, and its apex sits on its own territory edge with the
      trace still climbing across that edge into ``j``'s territory;
    * climbing on from there, within ``max_reach`` of ``i``'s anchor, ends exactly on ``j``'s
      apex, and ``j`` passed the presence test there;
    * ``i`` is a major congener (``MAJOR_FRACTION``) and at least ``dominance_ratio`` times
      more abundant than ``j`` in the published Aroclor profiles. Where the table marks
      ``j``'s abundance unknown, the top must instead lie more than half the apex search
      window from ``j``'s own anchor.

    ``apexes`` and ``present`` are updated in place: ``i`` becomes present at the top, and
    ``j`` absent at its own anchor. Every boundary in the transition is then resolved with
    one consistent verdict. Returns ``{i: j}`` for every peak reclaimed.

    Where the drift is just as clear but the two are both major and neither is
    ``dominance_ratio`` times the other, the profiles cannot say whose peak it is, because it
    is both: the two congeners merged. The envelope is then split by a rigid two-peak fit
    (``core.deconvolve.fit_pair``), and nothing in ``apexes`` or ``present`` changes. Those
    come back as the second mapping, ``{i: (j, fit)}``.
    """
    reclaimed: dict[int, int] = {}
    shared: dict[int, tuple[int, PairFit]] = {}
    if params.dominance_ratio <= 0:
        return reclaimed, shared
    abundance = load_aroclor_abundance()
    last = detect.size - 1
    reach = _halfwidth_points(times, params.max_reach)
    for i, analyte in enumerate(ordered):
        mine = abundance.get(analyte.label, 0.0)
        # Written so that an unknown (nan) abundance fails it: a congener whose abundance
        # the profiles do not give can lose a peak to a dominant neighbour, never take one.
        if present[i] or not mine >= MAJOR_FRACTION:
            continue
        lo_b, hi_b = bounds[i]
        for j, edge, step in ((i + 1, hi_b, 1), (i - 1, lo_b, -1)):
            if not (0 <= j < len(ordered)) or not present[j]:
                continue
            theirs_j = abundance.get(ordered[j].label, 0.0)
            if (
                not np.isnan(theirs_j)
                and theirs_j >= MAJOR_FRACTION
                and mine < params.dominance_ratio * theirs_j
                and theirs_j < params.dominance_ratio * mine
            ):
                # Comparable congeners: the fit, with its own quality gates, is the evidence.
                if j not in {k for k, _ in shared.values()} and j not in shared and j not in reclaimed:
                    pair = fit_pair(times, integ, expected[i], expected[j], params.max_reach)
                    if pair is not None:
                        shared[i] = (j, pair)
                        break
                continue
            if apexes[i] != edge:
                continue
            if j in reclaimed.values() or j in reclaimed:
                continue
            if not (0 <= edge + step <= last) or detect[edge + step] <= detect[edge]:
                continue
            theirs = abundance.get(ordered[j].label, 0.0)
            unknown = np.isnan(theirs)
            if not unknown and mine < params.dominance_ratio * theirs:
                continue
            centre = int(np.searchsorted(times, expected[i]))
            top = hill_climb(detect, edge, max(0, centre - reach), min(last, centre + reach))
            if top != apexes[j]:
                continue
            # A congener whose abundance the profiles leave unknown (the table marks it nan) is
            # compared by position instead: the peak must sit away from that congener's own
            # retention time, as a drifted peak does and the congener's own peak does not.
            if unknown and abs(float(times[top]) - expected[j]) <= params.rt_window / 2:
                continue
            j_lo, j_hi = bounds[j]
            apexes[i], present[i], present[j] = top, True, False
            apexes[j] = min(max(int(np.searchsorted(times, expected[j])), j_lo), j_hi)
            reclaimed[i] = j
            break
    return reclaimed, shared


def _fit_rider(
    times: np.ndarray,
    integ: np.ndarray,
    detect: np.ndarray,
    apexes: list[int],
    present: list[bool],
    i: int,
    rider_rt: float,
    halfwidth: int,
) -> RiderFit | None:
    """Fit analyte ``i`` as a rider on the resolved peaks around it, or None.

    The envelope is every resolved peak within three half-widths of the rider. Its parent,
    the peak whose measured width seeds the fit, is the tallest of them. With no resolved
    peak in reach there is nothing for a rider to be riding on, and nothing is fitted.
    """
    near = [
        j for j in range(len(apexes))
        if present[j] and abs(apexes[j] - apexes[i]) <= 3 * halfwidth
    ]
    if not near:
        return None
    parent = max(near, key=lambda j: detect[apexes[j]])
    return deconvolve_rider(
        times, integ, detect,
        rider_rt=rider_rt,
        neighbour_rts=[float(times[apexes[j]]) for j in near],
        parent_apex=apexes[parent],
    )


def _force(
    times: np.ndarray,
    integ: np.ndarray,
    detect: np.ndarray,
    ordered: list[Analyte],
    expected: list[float],
    bounds: list[tuple[int, int]],
    apexes: list[int],
    present: list[bool],
    rescued: set[int],
    halfwidth: int,
    params: IntegrationParams,
) -> dict[int, _Forced]:
    """Apply the named rules to the analytes neither detection test found.

    Returns the analytes a rule integrated, by position in ``ordered``. An analyte the rule
    could not integrate honestly is left out and stays ``absent``: a fit that does not
    reproduce the trace, a boundary that lands on the apex, a region with nothing above its
    baseline. A forced peak is reported as detected, so it has to meet the same record
    contract as any other detected peak (start < apex < end, and a positive area).

    Like the shoulder rescue, this may add a peak but never alters one. Every other analyte
    is still resolved against the unmodified presence verdict and the unmodified apex list.
    """
    forced: dict[int, _Forced] = {}
    last = times.size - 1
    for i, analyte in enumerate(ordered):
        if present[i] or i in rescued:
            continue
        if analyte.label in params.forced_riders:
            # The rider is sought where the run's own calibration injection put it. The
            # refined apex is no guide here: with no maximum of its own to stop at, the search
            # ends on the territory edge, half way to the parent, and a fit started there
            # credits the rider with the parent's flank. Measured on the archive, starting at
            # the anchor brings the fitted total for both named riders from about 1.3 times
            # the analyst's to within 7% of it, at the same median agreement.
            fit = _fit_rider(times, integ, detect, apexes, present, i, expected[i], halfwidth)
            if fit is None:
                continue
            apex = int(np.argmin(np.abs(times - fit.centre)))
            lo_b, hi_b = bounds[i]
            start, end = max(0, min(lo_b, apex - 1)), min(last, max(hi_b, apex + 1))
            if start < apex < end:
                forced[i] = _Forced(start, apex, end, "deconvolved", fit)
        elif analyte.label in params.forced_present:
            as_if = list(present)
            as_if[i] = True
            start, end, method = _resolve_boundaries(detect, bounds, apexes, as_if, i, times,
                                                     params.max_halfwidth)
            # Where the search ended on a slope, its territory cut the climb short: the top
            # of the peak lies further on. The boundaries above are the rule's, and the apex
            # recorded is the top actually reached inside them, so the height and S/N
            # describe the peak and not a point on its flank.
            apex = hill_climb(detect, apexes[i], start, end)
            if start < apex < end and integrate_region(times, integ, start, end)[0] > 0.0:
                forced[i] = _Forced(start, apex, end, method, None)
    return forced


def integrate_transition(
    times: np.ndarray,
    intensities: np.ndarray,
    analytes: list[Analyte],
    anchors: dict[str, float],
    params: IntegrationParams,
    src_q1: float = -1.0,
    src_q3: float = -1.0,
    unresolved: frozenset[tuple[str, str]] = frozenset(),
) -> list[PeakResult]:
    """Integrate every analyte that shares one transition trace.

    ``src_q1``/``src_q3`` identify the acquired transition these arrays came from and are
    recorded on each ``PeakResult`` so a re-plot can find the trace again. The canonical
    transition on the analyte may not be a reader key. See ``core.mzml.get_transition``.

    ``unresolved`` holds the label pairs this run's calibration injection could not separate
    (:func:`unresolved_pairs`). Where only one of such a pair has a peak in this file, the
    envelope is fitted as two components and each takes its own (flag ``unresolved``), unless
    one peak explains it as well. Empty, the default, changes nothing.

    ``times``/``intensities`` are never written to. The same trace is integrated once per
    analyte sharing the transition, so any in-place edit would corrupt its neighbours.
    """
    if times.size == 0:
        return []

    ordered = [a for a in analytes if not np.isnan(anchors.get(a.label, np.nan))]
    ordered.sort(key=lambda a: anchors[a.label])
    if not ordered:
        return []

    # 1) Refine each apex within its midpoint-bounded RT territory. Detection runs on the
    #    smoothed trace (stable apexes/valleys). The area is measured on whichever trace
    #    ``integrate_on_smoothed`` selects, raw by default, MassHunter-like.
    expected = [anchors[a.label] for a in ordered]
    detect, bounds, apexes = refine_apexes(times, intensities, expected, params)
    integ = detect if params.integrate_on_smoothed else intensities.astype(float)
    halfwidth = _halfwidth_points(times, params.max_halfwidth)

    # 2) Presence test: topographic prominence must clear an absolute floor *and* a noise-scaled
    #    floor. Both halves matter. Raw-trace prominence alone can exceed ``min_prominence`` on
    #    a flat noisy trace, and the S/N term alone lets a tiny bump on a quiet trace masquerade
    #    as a resolved peak. "Present" here means "resolved enough to split a boundary against".
    #    It also drives the confidence flags below.
    #
    #    The asymmetry below is intentional: prominence is measured on the *smoothed* detect
    #    trace (so a single noisy sample cannot manufacture prominence), while noise is
    #    measured on the *integration* trace. Setting ``integrate_on_smoothed`` therefore does
    #    not merely change which trace is integrated. It also lowers the noise estimate and so
    #    lowers the presence threshold, admitting peaks that the default settings reject. The
    #    two parameters are coupled. Retune presence if you change it.
    #
    #    ``peak_prominences`` returns 0 for an apex at either end of the trace, so an analyte
    #    whose apex lands on a trace boundary is forced absent by construction. Not observed on
    #    the reference data, since territories are interior, but it is the behaviour by design.
    noise = global_noise(integ)
    prominences = peak_prominences(detect, apexes)
    prominence_threshold = max(params.min_prominence, params.presence_sn * noise)
    present = [bool(prominences[i] >= prominence_threshold) for i in range(len(apexes))]

    #    2a) A dominant peak displaced across its territory edge is given back to the analyte
    #    whose search it drifted out of, before anything else reads ``present`` or ``apexes``
    #    (see ``_reclaim_displaced``). The neighbour that had claimed it is then absent, so the
    #    shoulder test or a named rule may still integrate it as a rider on the reclaimed peak.
    reclaimed, shared = _reclaim_displaced(
        times, detect, integ, ordered, expected, bounds, apexes, present, params
    )
    partners = {j: i for i, (j, _) in shared.items()}

    #    2a') A pair this run's calibration could not resolve (``unresolved``). Where only one of
    #    the two has a peak in this file, that peak is an envelope the two may share. It is
    #    fitted as the two congeners at the spacing the calibration measured between them
    #    (their anchors, from ``core.deconvolve.fit_unresolved_pair`` on the calibration), with
    #    one common drift (``core.deconvolve.fit_pair``), and each takes its own component, as
    #    a comparable pair does. A congener alone fits as one component and nothing, which the
    #    fit refuses (``MIN_PAIR_SHARE``), so a single peak is never split. Like a comparable
    #    pair, the absent one is then not offered to the shoulder test: a skim models a rider
    #    on a parent, not two congeners the calibration showed eluting together.
    unresolved_members: set[int] = set()
    if unresolved:
        index = {a.label: k for k, a in enumerate(ordered)}
        taken = set(shared) | set(partners) | set(reclaimed) | set(reclaimed.values())
        for one, two in sorted(unresolved):
            if one not in index or two not in index or abs(index[one] - index[two]) != 1:
                continue
            for i, j in ((index[one], index[two]), (index[two], index[one])):
                if present[i] or not present[j] or {i, j} & (taken | unresolved_members):
                    continue
                # The envelope is the present one's peak, bounded as it would be without the rule.
                j_start, j_end, _ = _resolve_boundaries(detect, bounds, apexes, present, j, times,
                                                        params.max_halfwidth)
                pair = fit_pair(times, integ, expected[i], expected[j], params.max_reach,
                                window=(j_start, j_end))
                if pair is None:
                    continue
                shared[i] = (j, pair)
                partners[j] = i
                unresolved_members |= {i, j}
                break


    #    2b) Shoulder rescue. An analyte eluting on the flank of a much larger neighbour has
    #    no local maximum of its own, so its prominence is ~0 and the test above rejects it.
    #    Yet the analyst integrates it, because the flank visibly bulges. Downward curvature
    #    is that bulge (see ``core.peaks.find_shoulder``). A rescued analyte is *present* for
    #    every purpose except that its area is skimmed off its parent's tangent rather than
    #    dropped to a baseline of its own, and it is flagged ``shoulder`` so a reviewer can
    #    see which rule produced the number.
    #
    #    Only analytes the prominence test rejected are examined, so this can add peaks but
    #    never remove or alter one. A run with no shoulders is untouched.
    shoulders: dict[int, tuple[int, int, int, float]] = {}
    if params.shoulder_curvature_sn > 0:
        min_curvature = params.shoulder_curvature_sn * noise
        # A neighbour that lost a reclaimed peak sits on that peak's flank, where the flank's
        # own curvature scales with the peak's height and clears any noise-scaled threshold.
        # The bulge test says nothing there, so it is not asked. A named rider is still fitted.
        robbed = set(reclaimed.values())
        for i in range(len(apexes)):
            if present[i] or i in robbed or i in shared:
                continue
            lo_b, hi_b = bounds[i]
            found = find_shoulder(detect, lo_b, hi_b, min_curvature)
            if found is None:
                continue
            i_left, i_apex, i_right = found
            # The rider must sit on a flank, not in flat noise: one side of its window has to
            # belong to a peak that cleared the presence test. Without this, a noise ripple in
            # an empty window would be promoted to a peak.
            beside_real_peak = (i > 0 and present[i - 1]) or (i + 1 < len(apexes) and present[i + 1])
            if not beside_real_peak:
                continue
            # Bound the rider: its foot on the side away from its parent, and a
            # perpendicular drop on the parent side. The parent is whichever neighbour
            # carries the taller peak. The drop falls half way from the rider's own apex to
            # the midpoint it already shares with that parent, which is where the parent's
            # rise has visibly taken over. Both are measured, not assumed. See
            # ``docs/validation/`` for the traces these were read off.
            left_h = detect[apexes[i - 1]] if i > 0 else 0.0
            right_h = detect[apexes[i + 1]] if i + 1 < len(apexes) else 0.0
            if right_h >= left_h:
                foot = lo_b + int(np.argmin(detect[lo_b:i_apex + 1])) if i_apex > lo_b else lo_b
                drop = min(hi_b, i_apex + max(1, (hi_b - i_apex) // 2))
                start, end = foot, drop
            else:
                foot = i_apex + int(np.argmin(detect[i_apex:hi_b + 1])) if hi_b > i_apex else hi_b
                drop = max(lo_b, i_apex - max(1, (i_apex - lo_b) // 2))
                start, end = drop, foot
            # A rescued rider is reported as a detected peak, so it must satisfy the same
            # auditability contract as any other: start < apex < end, or it cannot be
            # re-plotted from its own record. Where the drop lands on the apex there is no
            # rider to integrate, and the analyte stays absent.
            if not (start < i_apex < end):
                continue
            shoulders[i] = (start, i_apex, end, float(detect[foot]))

    #    2c) Named rules, for the few analytes that both tests above miss in most files and
    #    the analyst integrates anyway (``IntegrationParams.forced_riders`` /
    #    ``forced_present``). Applied last and only to what is still unfound, so a file in
    #    which the trace itself shows the peak is integrated on that evidence, not by rule.
    forced = _force(
        times, integ, detect, ordered, expected, bounds, apexes, present,
        set(shoulders) | set(shared), halfwidth, params,
    )

    # 3) Per-analyte boundaries (see ``_resolve_boundaries``), then the span each analyte's
    #    baseline rests on: its own boundaries, or a whole unresolved cluster's
    #    (see ``_cluster_spans``). Boundaries are decided first and are not affected by
    #    clustering. Only the baseline drawn under them is.
    #    A rescued shoulder takes its own inflection bounds and a forced analyte its rule's,
    #    but every *other* analyte is resolved against the unmodified presence verdict and
    #    the unmodified apex list. A rescue must not move a neighbour's valley: these
    #    branches may add a peak, never silently re-cut one that already agreed with the
    #    analyst. Neither kind joins a cluster, because ``present`` is not changed for them.
    resolved = [
        (shoulders[i][0], shoulders[i][2], "tangent_skim")
        if i in shoulders
        else (forced[i].start, forced[i].end, forced[i].method)
        if i in forced
        else _resolve_boundaries(detect, bounds, apexes, present, i, times, params.max_halfwidth)
        for i in range(len(ordered))
    ]
    #    A neighbour that lost a reclaimed peak is absent, and an absent analyte is bounded
    #    inside its own territory, which the reclaimed peak now crosses. Its window stops where
    #    the reclaimed peak's starts, so no area is reported twice. A rider integrated on the
    #    reclaimed peak's flank keeps its own region, as every rider does.
    for i, j in reclaimed.items():
        if j in shoulders or j in forced:
            continue
        j_start, j_end, j_method = resolved[j]
        i_start, i_end, _ = resolved[i]
        if j > i:
            j_start = max(j_start, i_end)
            j_end = max(j_end, j_start)
        else:
            j_end = min(j_end, i_start)
            j_start = min(j_start, j_end)
        resolved[j] = (j_start, j_end, j_method)
    spans = _cluster_spans(detect, resolved, apexes, present, params.overlap_frac)


    results: list[PeakResult] = []
    for i, analyte in enumerate(ordered):
        if i in shoulders:
            apex = shoulders[i][1]
        elif i in forced:
            apex = forced[i].apex
        else:
            apex = apexes[i]
        i_start, i_end, method = resolved[i]
        base_l, base_r = spans[i]
        fit = forced[i].fit if i in forced else None
        if i in shared:
            # The analyte that was absent takes its share of the envelope it merged into. Its
            # record spans that envelope, and its apex is where the fit placed it.
            fit = shared[i][1].first
            i_start, i_end = resolved[shared[i][0]][0], resolved[shared[i][0]][1]
            apex = min(max(int(np.argmin(np.abs(times - fit.centre))), i_start + 1), i_end - 1)
            method = "deconvolved"
        elif i in partners:
            fit = shared[partners[i]][1].second
            method = "deconvolved"

        clustered = (base_l < i_start or base_r > i_end) and i not in shoulders and fit is None
        if fit is not None:
            # A forced rider. The fit supplies the area. The baseline recorded is the chord
            # across the window it was looked for in, which is what a reviewer is shown.
            _, baseline_start, baseline_end = integrate_region(times, integ, i_start, i_end)
            area = fit.area
        elif i in shoulders:
            area, baseline_start, baseline_end = integrate_skim(
                times, integ, i_start, i_end, shoulders[i][3]
            )
            # A skim approximates a rider. Fitting its shape against its parent's actually
            # separates them, where the fit reproduces the trace. ``deconvolve_rider``
            # returns None when it does not, and the skim above stands.
            if params.deconvolve_shoulders:
                refit = _fit_rider(
                    times, integ, detect, apexes, present, i, float(times[apex]), halfwidth
                )
                if refit is not None:
                    area = refit.area
                    method = "deconvolved"
        elif clustered:
            area, baseline_start, baseline_end = integrate_above_baseline(
                times, integ, i_start, i_end, base_l, base_r
            )
        else:
            area, baseline_start, baseline_end = integrate_region(times, integ, i_start, i_end)
        if fit is not None:
            # A forced rider has no apex on the trace to measure: any one sample there is
            # mostly its parent. Its height is its fitted component's, so its S/N grades the
            # rider and not the flank it sits on.
            height = fit.height
        else:
            # The cluster baseline is a straight line, so interpolating between its values at
            # this peak's two boundaries evaluates it exactly at the apex. The height, and
            # therefore ``sn`` and the ``low_sn`` flag, are measured above whichever baseline
            # the area was.
            height = float(
                integ[apex]
                - np.interp(
                    times[apex],
                    [times[i_start], times[i_end]],
                    [baseline_start, baseline_end],
                )
            )
        sn = height / noise  # global_noise never returns 0 (it floors at 1.0)
        # Appended in the order ``KNOWN_FLAGS`` documents, so two peaks carrying the same
        # flags always spell them the same way.
        flags: list[str] = []
        if not present[i] and i not in shoulders and i not in forced and i not in shared:
            flags.append("absent")
        if sn < params.min_sn:
            flags.append("low_sn")
        if i_end - i_start < 3:
            flags.append("narrow")
        if clustered:
            flags.append("clustered")
        if i in shoulders:
            flags.append("shoulder")
        if i in forced:
            flags.append("forced")
        if i in reclaimed:
            flags.append("displaced")
        if (i in shared or i in partners) and i not in unresolved_members:
            flags.append("apportioned")
        if i in unresolved_members:
            flags.append("unresolved")
        results.append(
            PeakResult(
                label=analyte.label,
                q1=analyte.q1,
                q3=analyte.q3,
                rt_apex=float(times[apex]),
                rt_start=float(times[i_start]),
                rt_end=float(times[i_end]),
                area=area,
                height=height,
                baseline_start=baseline_start,
                baseline_end=baseline_end,
                sn=sn,
                method=method,
                flags=flags,
                i_start=i_start,
                i_apex=apex,
                i_end=i_end,
                src_q1=src_q1,
                src_q3=src_q3,
            )
        )
    return results


def integrate_run(
    chroms: dict[tuple[float, float], Chromatogram],
    analytes: list[Analyte],
    anchors: dict[str, float],
    params: IntegrationParams | None = None,
    unresolved: frozenset[tuple[str, str]] = frozenset(),
) -> list[PeakResult]:
    """Integrate every analyte in one data file against a per-run RT anchor table.

    ``unresolved``: the label pairs this run's calibration could not separate
    (:func:`unresolved_pairs`), passed to every transition (see :func:`integrate_transition`).

    ``chroms`` is a reader dict as returned by :func:`core.mzml.read_srm_mzml`. ``anchors``
    maps analyte label -> anchor RT (minutes), normally from :func:`core.calibrate.build_anchors`.
    Analytes with no anchor, or whose transition is absent from the file, are skipped (with a
    warning) rather than reported at a guessed retention time.

    An inseparable co-elution group is a single row in the canonical table with a single RT,
    so it falls out of this loop as exactly one summed peak. No special case is needed, and
    none may be added. Splitting a group would report congeners the chromatography cannot
    resolve.

    **Reading the results.** ``flags`` is the detection verdict. ``area``, ``height`` and
    ``sn`` are conditional on it. A result flagged ``absent`` was rejected by the presence
    test: its area is a bounded estimate of whatever occupies the analyte's RT territory, its
    apex may lie outside the integrated region, and its ``sn`` is not meaningful. Absent
    peaks frequently report S/N well above ``min_sn`` and so are usually not also flagged
    ``low_sn``. Gate on :attr:`PeakResult.detected`, never on ``sn``. Use ``src_q1``/``src_q3``
    (not ``q1``/``q3``) to look the source trace back up in ``chroms``.

    ``area``, ``height`` and both baseline endpoints refer to the raw trace unless
    ``params.integrate_on_smoothed`` is set, in which case every one of them is measured on
    the smoothed trace instead. Areas are counts·seconds.

    Results are only reproducible alongside the ``IntegrationParams`` that produced them.
    This function does not attach them, so callers that persist results must record the
    params used.
    """
    params = params or IntegrationParams()
    results: list[PeakResult] = []
    for transition, transition_analytes in analytes_by_transition(analytes).items():
        chrom = get_transition(chroms, *transition)
        if chrom is None:
            continue
        results.extend(
            integrate_transition(
                chrom.times, chrom.intensities, transition_analytes, anchors, params,
                src_q1=chrom.q1, src_q3=chrom.q3, unresolved=unresolved,
            )
        )

    # A run's output silently shrinking against the canonical table is the kind of loss that
    # only shows up as a blank column in a report weeks later. Say so at the point it happens.
    skipped = sorted({a.label for a in analytes} - {r.label for r in results})
    if skipped:
        warnings.warn(
            f"{len(skipped)} of {len(analytes)} analytes were not integrated "
            f"(no RT anchor, or transition missing from this file): {', '.join(skipped)}",
            RuntimeWarning,
            stacklevel=2,
        )
    return results


def unresolved_pairs(
    cal_peaks: list[PeakResult],
    analytes: list[Analyte],
    anchors: dict[str, float],
    params: IntegrationParams | None = None,
) -> frozenset[tuple[str, str]]:
    """Neighbouring congeners this run's calibration injection could not separate.

    The calibration standard holds every congener, so in its injection an analyte with no
    apex of its own is not absent: it is inside a neighbour's peak. Where its anchor lies
    within that neighbour's integrated peak, and within the apex search window of the
    neighbour's apex, the two eluted as one envelope, and the chromatography of this run does
    not resolve them. Integrated as it stands, the whole envelope goes to one of the pair and
    the other gets nothing: no response factor, and so no mass in any file of the run.

    ``cal_peaks`` is the calibration file integrated at its own ``anchors`` with no pairs
    given. Returns ``(earlier, later)`` label pairs by anchor, for :func:`integrate_run`.
    A run whose calibration resolves every pair returns an empty set and integrates exactly
    as before.
    """
    params = params or IntegrationParams()
    by_label = {p.label: p for p in cal_peaks}
    pairs: set[tuple[str, str]] = set()
    for _transition, group in analytes_by_transition(analytes).items():
        ordered = sorted(
            (a for a in group if a.role == "native" and a.label in by_label and a.label in anchors),
            key=lambda a: anchors[a.label],
        )
        for first, second in zip(ordered, ordered[1:]):
            for lost, kept in ((first, second), (second, first)):
                p_lost, p_kept = by_label[lost.label], by_label[kept.label]
                if p_lost.detected or not p_kept.detected or "absent" not in p_lost.flags:
                    continue
                at = anchors[lost.label]
                if p_kept.rt_start < at < p_kept.rt_end and abs(at - p_kept.rt_apex) <= params.rt_window:
                    pairs.add((first.label, second.label))
    return frozenset(pairs)
