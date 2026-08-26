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

import numpy as np

from .config import Analyte, analytes_by_transition
from .deconvolve import deconvolve_rider
from .mzml import Chromatogram, get_transition
from .peaks import (
    IntegrationParams,
    PeakResult,
    _local_min,
    global_noise,
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
KNOWN_METHODS = frozenset({"valley_drop", "rt_window", "tangent_skim", "deconvolved"})

#: Every flag a ``PeakResult`` may carry. The first three are confidence verdicts.
#: ``flags`` is the detection verdict and ``area``/``height``/``sn`` are conditional on it
#: (see ``PeakResult.detected``). The last records how the area's baseline was drawn.
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
KNOWN_FLAGS = frozenset({"absent", "low_sn", "narrow", "clustered", "shoulder"})


def _halfwidth_points(times: np.ndarray, max_halfwidth: float) -> int:
    """Convert a max half-width in minutes to a sample count for this trace."""
    if times.size < 2:
        return 5
    dt = float(np.median(np.diff(times)))
    return max(2, int(round(max_halfwidth / dt))) if dt > 0 else 5


def _resolve_boundaries(
    detect: np.ndarray,
    bounds: list[tuple[int, int]],
    apexes: list[int],
    present: list[bool],
    i: int,
    halfwidth: int,
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
      territory to its own baseline edge, capped at one half-width from the apex.
    """
    n = detect.size
    apex = apexes[i]

    if not present[i]:
        lo_b, hi_b = bounds[i]
        mid = (lo_b + hi_b) // 2
        return _local_min(detect, lo_b, mid, lo_b), _local_min(detect, mid, hi_b, hi_b), "rt_window"

    left_resolved = i > 0 and present[i - 1]
    right_resolved = i < len(apexes) - 1 and present[i + 1]

    left_limit = valley_index(detect, apexes[i - 1], apex) if left_resolved else max(0, apex - halfwidth)
    i_start = _local_min(detect, left_limit, apex, left_limit)

    right_limit = (
        valley_index(detect, apex, apexes[i + 1]) if right_resolved else min(n - 1, apex + halfwidth)
    )
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
    in ``docs/METHODS.md`` section 4.

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


def integrate_transition(
    times: np.ndarray,
    intensities: np.ndarray,
    analytes: list[Analyte],
    anchors: dict[str, float],
    params: IntegrationParams,
    src_q1: float = -1.0,
    src_q3: float = -1.0,
) -> list[PeakResult]:
    """Integrate every analyte that shares one transition trace.

    ``src_q1``/``src_q3`` identify the acquired transition these arrays came from and are
    recorded on each ``PeakResult`` so a re-plot can find the trace again. The canonical
    transition on the analyte may not be a reader key. See ``core.mzml.get_transition``.

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
        for i in range(len(apexes)):
            if present[i]:
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

    # 3) Per-analyte boundaries (see ``_resolve_boundaries``), then the span each analyte's
    #    baseline rests on: its own boundaries, or a whole unresolved cluster's
    #    (see ``_cluster_spans``). Boundaries are decided first and are not affected by
    #    clustering. Only the baseline drawn under them is.
    #    A rescued shoulder takes its own inflection bounds, but every *other* analyte is
    #    resolved against the unmodified presence verdict and the unmodified apex list. A
    #    rescue must not move a neighbour's valley: the shoulder branch may add a peak, never
    #    silently re-cut one that already agreed with the analyst.
    resolved = [
        (shoulders[i][0], shoulders[i][2], "tangent_skim")
        if i in shoulders
        else _resolve_boundaries(detect, bounds, apexes, present, i, halfwidth)
        for i in range(len(ordered))
    ]
    spans = _cluster_spans(detect, resolved, apexes, present, params.overlap_frac)

    results: list[PeakResult] = []
    for i, analyte in enumerate(ordered):
        apex = shoulders[i][1] if i in shoulders else apexes[i]
        i_start, i_end, method = resolved[i]
        base_l, base_r = spans[i]

        clustered = (base_l < i_start or base_r > i_end) and i not in shoulders
        if i in shoulders:
            area, baseline_start, baseline_end = integrate_skim(
                times, integ, i_start, i_end, shoulders[i][3]
            )
            # A skim approximates a rider. Fitting its shape against its parent's actually
            # separates them, where the fit reproduces the trace. ``deconvolve_rider``
            # returns None when it does not, and the skim above stands.
            if params.deconvolve_shoulders:
                neighbours = [
                    float(times[apexes[j]])
                    for j in range(len(ordered))
                    if present[j] and abs(apexes[j] - apexes[i]) <= 3 * halfwidth
                ]
                parent = max(
                    (j for j in range(len(ordered)) if present[j]),
                    key=lambda j: detect[apexes[j]] if abs(apexes[j] - apexes[i]) <= 3 * halfwidth else -1.0,
                    default=None,
                )
                if neighbours and parent is not None:
                    fitted = deconvolve_rider(
                        times, integ, detect,
                        rider_rt=float(times[shoulders[i][1]]),
                        neighbour_rts=neighbours,
                        parent_apex=apexes[parent],
                    )
                    if fitted is not None:
                        area = fitted[0]
                        method = "deconvolved"
        elif clustered:
            area, baseline_start, baseline_end = integrate_above_baseline(
                times, integ, i_start, i_end, base_l, base_r
            )
        else:
            area, baseline_start, baseline_end = integrate_region(times, integ, i_start, i_end)
        # The cluster baseline is a straight line, so interpolating between its values at this
        # peak's two boundaries evaluates it exactly at the apex. The height, and therefore
        # ``sn`` and the ``low_sn`` flag, are measured above whichever baseline the area was.
        height = float(
            integ[apex]
            - np.interp(
                times[apex],
                [times[i_start], times[i_end]],
                [baseline_start, baseline_end],
            )
        )
        sn = height / noise  # global_noise never returns 0 (it floors at 1.0)
        flags: list[str] = []
        if not present[i] and i not in shoulders:
            flags.append("absent")
        if sn < params.min_sn:
            flags.append("low_sn")
        if i_end - i_start < 3:
            flags.append("narrow")
        if clustered:
            flags.append("clustered")
        if i in shoulders:
            flags.append("shoulder")
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
) -> list[PeakResult]:
    """Integrate every analyte in one data file against a per-run RT anchor table.

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
                src_q1=chrom.q1, src_q3=chrom.q3,
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
