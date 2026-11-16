"""Core peak-integration primitives: apex detection, boundary walking, baselines, area.

These are explicit and tunable rather than fitted by an opaque curve, so every integrated
area is reconstructable. Given the persisted RT window, boundaries, and baseline, an
analyst can re-plot exactly what was integrated and why.

These primitives are numerically validated against analyst ground truth. They use numpy
only, with no scipy dependency, keeping them cheap and reliable on a memory-constrained
host. Defaults in ``IntegrationParams`` are calibrated against analyst manual integration.
Changing them moves results away from that calibration.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class IntegrationParams:
    """Tunable integration parameters (calibrated against analyst ground truth)."""

    rt_window: float = 0.10  # ± minutes for the initial apex guess around the expected RT
    max_reach: float = 0.20  # ± minutes an isolated peak's territory may extend (drift cap)
    max_halfwidth: float = 0.25  # max peak half-width (min) when no neighbour bounds the peak
    smooth_window: int = 5  # Savitzky-Golay window for apex/boundary detection, off under 3
    smooth_poly: int = 2
    min_sn: float = 3.0  # below this S/N the peak is flagged low-confidence
    presence_sn: float = 4.0  # prominence S/N above which a peak is "present" (splits neighbours)
    min_prominence: float = 150.0  # absolute prominence floor for presence
    # Adjacent *present* peaks merge into one baseline cluster when the valley between them
    # stands above this fraction of the smaller apex's height over the pair's outer bases.
    # That is, when the trace never returns to the chromatographic floor between them.
    # Consumed by ``core.integrate._cluster_spans``. See its docstring for the rule and for
    # why the present-only restriction matters. Without it, a floor-anchored cluster baseline
    # also fires on *shoulders*, handing a small shoulder its dominant neighbour's underlying
    # area, whereas the analyst skims shoulders. The shoulder case is excluded here.
    #
    # Over 0.02-0.15 the *outcome* is flat. The count of response factors more than 25% from
    # the analyst's does not move on five of the six (injection, workbook) pairs of the two
    # reference runs, and moves by one on the sixth, while the merge set itself shrinks
    # monotonically across that range (37 -> 20 merged pairs on one run's Cal01). The plateau is in
    # the metric, not the mechanism, so 0.05 is a point on the conservative half of a range
    # whose endpoints both pass, not a knife edge and not a value the data pins exactly.
    # Regenerate the evidence with
    # ``python tools/bench_quant.py --sweep overlap_frac=10,0.02,0.05,0.10,0.15``.
    # Table in docs/METHODS.md section 2.4.
    overlap_frac: float = 0.05
    integrate_on_smoothed: bool = False  # integrate the raw signal by default (MassHunter-like)
    #: Shoulder test: how strong the downward curvature at an analyte's apex must be, as a
    #: multiple of the trace noise, before a peak the prominence test rejected is treated as
    #: a rider on its neighbour's flank rather than as absent. Set to 0 to disable the
    #: shoulder branch entirely and recover the pre-shoulder engine.
    #:
    #: Congeners 83, 51 and 115 are the reason this exists: each elutes 0.02-0.06 min from a
    #: much larger neighbour and each was reported absent in 91-97% of the archive's
    #: detectable rows, against analyst areas with a median of 3,500-26,500 counts·s.
    shoulder_curvature_sn: float = 6.0
    #: Resolve a rescued rider by fitting its shape against its parent's
    #: (``core.deconvolve``) instead of skimming it off a flat baseline. A rider on a
    #: *leading* edge has no valley to drop at, so a boundary rule can only approximate it.
    #: Fitting two peak shapes to the envelope is what actually separates them. The fit is
    #: used only when it reproduces the trace. Otherwise the skim stands. Set False to
    #: disable and keep the skim everywhere.
    deconvolve_shoulders: bool = True
    #: Named rules: analytes integrated by a fixed rule in any file where neither the presence
    #: test nor the shoulder test finds them. Every other rule in this engine reads the trace
    #: and decides. These two read the label, so each list is an explicit, per-method statement
    #: about the chromatography, kept short and kept here where a different column or method
    #: can empty it. A peak either rule produced carries the ``forced`` flag. Consumed by
    #: ``core.integrate.integrate_transition``. See ``docs/METHODS.md`` section 2.7 for the
    #: evidence, including where each rule does not help.
    #:
    #: ``forced_riders`` are fitted against their neighbours (``core.deconvolve``) without the
    #: curvature evidence the shoulder test asks for. Congeners 83 and 117 ride so close to a
    #: larger neighbour that most files show no bulge at all, and were reported absent on most
    #: of the archive's detectable rows. The fit's own quality gate still applies: a fit that
    #: does not reproduce the trace is refused and the analyte stays absent.
    #:
    #: ``forced_present`` are bounded as if the presence test had accepted them. The two
    #: listed are dominant congeners of the technical mixtures. In a concentrated sample the
    #: apex of such a peak leaves its search window, the point examined in its place has no
    #: prominence, and the largest peak in the chromatogram is reported absent.
    #:
    #: Set either to ``()`` to turn that rule off.
    forced_riders: tuple[str, ...] = ("83", "117")
    forced_present: tuple[str, ...] = ("90+101+113", "118")
    #: Displaced dominant peaks. In a concentrated sample a dominant peak broadens and its top
    #: can drift across the boundary it shares with a neighbour, so the apex search stops on
    #: its rising flank and the neighbour is credited with the whole peak. When an analyte's
    #: search ends on its territory edge with the trace still climbing into the neighbour's,
    #: and that climb leads to the very top the neighbour claimed, the peak is given to the
    #: analyte if it is at least this many times more abundant than the neighbour in the
    #: published Aroclor profiles (``core.config.load_aroclor_abundance``) and is itself a
    #: major congener. The two are then resolved together, so no area is counted twice, and
    #: the analyte carries the ``displaced`` flag. Set to 0 to turn the rule off. See
    #: ``docs/METHODS.md`` section 2.8.
    dominance_ratio: float = 10.0


@dataclass
class PeakResult:
    """A fully-auditable integrated peak.

    ``flags`` is the detection verdict. ``sn``, ``height`` and ``area`` are conditional on it.
    See :attr:`sn` and :meth:`detected` before using any of the three as evidence that a peak
    is real.
    """

    label: str
    q1: float  # canonical-table precursor m/z for this analyte
    q3: float  # canonical-table product m/z for this analyte
    rt_apex: float
    rt_start: float
    rt_end: float
    area: float
    height: float
    baseline_start: float
    baseline_end: float
    # Signal-to-noise of this peak's baseline-corrected apex height.
    #
    # MEANINGFUL ONLY WHEN ``detected`` IS TRUE. For an analyte the presence test rejected
    # ("absent" in ``flags``) the integrated region is the analyte's RT territory rather than a
    # region around its apex, so the apex may sit *outside* [i_start, i_end]. ``height`` then
    # extrapolates the region's straight baseline out to that apex and can come out large and
    # positive. Absent peaks routinely report S/N far above ``min_sn`` as a result (observed up
    # to ~4.9e4 on a calibration injection), so they do NOT pick up a "low_sn" flag. Never infer
    # presence from ``sn``. Gate on ``detected`` first, then use ``sn`` to grade the peak.
    sn: float
    method: str  # boundary-selection method, see ``core.integrate.KNOWN_METHODS``
    flags: list[str] = field(default_factory=list)  # see ``core.integrate.KNOWN_FLAGS``
    # index bounds into the transition trace (for re-plotting)
    i_start: int = -1
    i_apex: int = -1
    i_end: int = -1
    # The acquired transition actually integrated, as keyed by ``core.mzml.read_srm_mzml``.
    # Use these, not q1/q3, for ``chroms[...]`` lookups when re-plotting. The canonical
    # table's labelled-standard transitions are isotope-shift-derived and can sit ~0.1 Da off
    # the acquired key, so q1/q3 alone will not join back to the reader's dict.
    #
    # -1.0 means "not recorded", following this class's ``i_start``/``i_apex``/``i_end``
    # convention. A NaN sentinel would read better but would silently break this dataclass's
    # default ``__eq__`` (NaN != NaN), which callers and tests rely on.
    src_q1: float = -1.0
    src_q3: float = -1.0

    @property
    def detected(self) -> bool:
        """Whether the presence test accepted this analyte as a resolved peak.

        The single question downstream code should ask before treating ``area``/``height``/
        ``sn`` as a measurement of a real peak. False means the reported area is a bounded
        estimate of whatever sits in the analyte's RT territory, not a peak.
        """
        return "absent" not in self.flags


def smooth(intensity: np.ndarray, window: int, poly: int) -> np.ndarray:
    """Savitzky-Golay smoothing for detection (numpy-only, raw signal left untouched).

    The SG smoother is a sliding local polynomial fit. The centre-value weights are the first
    row of the pseudo-inverse of the Vandermonde design matrix. Convolution then uses edge
    padding. Implemented directly to avoid the heavy/fragile scipy.signal import on
    constrained hosts.
    """
    y = intensity.astype(float, copy=True)
    if not window or window < 3 or window >= y.size:
        return y
    w = window if window % 2 == 1 else window + 1
    if poly >= w:
        return y
    half = w // 2
    x = np.arange(-half, half + 1)
    design = x[:, None] ** np.arange(poly + 1)[None, :]
    weights = np.linalg.pinv(design)[0]  # weights producing the smoothed centre value
    padded = np.pad(y, half, mode="edge")
    return np.convolve(padded, weights[::-1], mode="valid")


def local_maxima(signal: np.ndarray) -> np.ndarray:
    """Indices of every interior local maximum of ``signal`` (a plateau reports its left edge).

    The unguided counterpart of :func:`find_apex`. Where ``find_apex`` asks "where is the peak
    near this expected RT", this asks "where are all the peaks", which is what
    ``core.calibrate.align_rt_map`` needs when the expected retention times are the very thing
    in doubt. Callers filter the (many) returned candidates by :func:`peak_prominences`.
    """
    if signal.size < 3:
        return np.empty(0, dtype=int)
    interior = signal[1:-1]
    return np.where((interior > signal[:-2]) & (interior >= signal[2:]))[0] + 1


def peak_prominences(signal: np.ndarray, peaks: list[int]) -> np.ndarray:
    """Topographic prominence of each peak index (numpy-only port of scipy's definition).

    Prominence = peak height minus the higher of the two minima reached by descending each
    way until the signal rises above the peak (or the trace ends).
    """
    n = signal.size
    proms = np.zeros(len(peaks), dtype=float)
    for k, p in enumerate(peaks):
        if p <= 0 or p >= n - 1:
            continue
        h = signal[p]
        i = p - 1
        left_min = h
        while i >= 0 and signal[i] < h:
            left_min = min(left_min, signal[i])
            i -= 1
        j = p + 1
        right_min = h
        while j < n and signal[j] < h:
            right_min = min(right_min, signal[j])
            j += 1
        proms[k] = h - max(left_min, right_min)
    return proms


def find_apex(time: np.ndarray, signal: np.ndarray, rt_expected: float, window: float,
              lo: int, hi: int) -> int:
    """Index of the maximum within the RT window, bounded by [lo, hi] (inclusive)."""
    mask = (time >= rt_expected - window) & (time <= rt_expected + window)
    idx = np.where(mask)[0]
    idx = idx[(idx >= lo) & (idx <= hi)]
    if idx.size == 0:
        # fall back to nearest sample to expected RT within [lo, hi]
        sub = np.arange(lo, hi + 1)
        return int(sub[np.argmin(np.abs(time[sub] - rt_expected))]) if sub.size else lo
    return int(idx[np.argmax(signal[idx])])


def hill_climb(signal: np.ndarray, i: int, lo: int, hi: int) -> int:
    """Climb from ``i`` to the nearest local maximum, staying within [lo, hi].

    Handles RT drift between the calibration anchor and a sample. If the windowed apex landed
    on a rising shoulder (because drift pushed the true apex just past the search window), this
    walks up to the real peak top. It stays bounded by the analyte's territory, so it cannot
    cross into a neighbouring peak.
    """
    while lo < i < hi:
        if signal[i + 1] > signal[i]:
            i += 1
        elif signal[i - 1] > signal[i]:
            i -= 1
        else:
            break
    return i


def territory_bounds(
    times: np.ndarray, expected_rts: list[float], max_reach: float
) -> list[tuple[int, int]]:
    """Each analyte's exclusive RT territory as inclusive index bounds.

    A territory runs out to the midpoint with each neighbour, capped at ± ``max_reach`` for
    an isolated peak (the drift cap). It is the region within which an analyte's apex and
    boundaries may roam. Because the bound is the *midpoint* to a neighbour, two analytes can
    never claim the same peak, while a peak whose true apex drifted past the narrow search
    window can still be reached by the hill-climb.

    ``expected_rts`` must be sorted ascending. The caller owns that ordering.
    """
    n = times.size
    bounds: list[tuple[int, int]] = []
    for i, rt_expected in enumerate(expected_rts):
        lo_rt, hi_rt = rt_expected - max_reach, rt_expected + max_reach
        if i > 0:
            lo_rt = max(lo_rt, 0.5 * (expected_rts[i - 1] + rt_expected))
        if i < len(expected_rts) - 1:
            hi_rt = min(hi_rt, 0.5 * (rt_expected + expected_rts[i + 1]))
        lo = int(np.searchsorted(times, lo_rt, side="left"))
        hi = int(np.searchsorted(times, hi_rt, side="right")) - 1
        lo = max(0, min(lo, n - 1))
        hi = max(lo, min(hi, n - 1))
        bounds.append((lo, hi))
    return bounds


def refine_apexes(
    times: np.ndarray,
    intensities: np.ndarray,
    expected_rts: list[float],
    params: IntegrationParams,
) -> tuple[np.ndarray, list[tuple[int, int]], list[int]]:
    """Locate each expected analyte's apex on one trace: smooth, territory, window, climb.

    Returns ``(detect, bounds, apexes)``: the smoothed detection trace, each analyte's
    territory as index bounds, and its refined apex index.

    This is the *single* definition of "where an analyte's apex is", shared by
    ``core.calibrate.build_anchors`` (measuring anchors off the calibration injection) and
    ``core.integrate`` (locating that analyte in every other file). The two must not drift
    apart. An anchor measured under one apex rule and consumed under another would silently
    bias every retention time in the run.

    ``expected_rts`` must be sorted ascending and index-aligned with the caller's analytes.
    ``intensities`` is not modified. The same trace is refined once per analyte sharing the
    transition.
    """
    detect = smooth(intensities, params.smooth_window, params.smooth_poly)
    bounds = territory_bounds(times, expected_rts, params.max_reach)
    apexes = [
        hill_climb(detect, find_apex(times, detect, expected_rts[i], params.rt_window, lo, hi), lo, hi)
        for i, (lo, hi) in enumerate(bounds)
    ]
    return detect, bounds, apexes


def _local_min(signal: np.ndarray, start: int, stop: int, default: int) -> int:
    """Index of the minimum of signal[start:stop+1], or ``default`` when the range is empty."""
    if stop < start:
        return default
    seg = signal[start:stop + 1]
    return int(start + np.argmin(seg))


def valley_index(signal: np.ndarray, i_left: int, i_right: int) -> int:
    """Index of the minimum (valley) between two apex indices."""
    if i_right <= i_left + 1:
        return i_left
    return int(i_left + 1 + np.argmin(signal[i_left + 1:i_right]))


def integrate_region(time: np.ndarray, intensity: np.ndarray, i_start: int, i_end: int
                     ) -> tuple[float, float, float]:
    """Trapezoidal area above a straight baseline joining the two boundary points.

    Returns (area, baseline_start_intensity, baseline_end_intensity).
    """
    if i_end <= i_start:
        return 0.0, float(intensity[i_start]), float(intensity[i_start])
    # Integrate in counts·seconds (Agilent/MassHunter convention). RT arrays are in minutes.
    t = time[i_start:i_end + 1] * 60.0
    y = intensity[i_start:i_end + 1]
    y0, y1 = float(y[0]), float(y[-1])
    # linear baseline across the region
    baseline = np.linspace(y0, y1, t.size)
    corrected = np.clip(y - baseline, 0.0, None)
    area = float(np.trapezoid(corrected, t))
    return area, y0, y1


def integrate_above_baseline(
    time: np.ndarray, intensity: np.ndarray,
    i_start: int, i_end: int, base_l: int, base_r: int,
) -> tuple[float, float, float]:
    """Area above a baseline that spans a whole merged cluster (true perpendicular drop).

    Called by ``core.integrate.integrate_transition`` for every analyte that
    ``core.integrate._cluster_spans`` placed in a multi-peak cluster, a run of consecutive
    *present* analytes the trace never descends to base between. An analyte outside such a
    cluster is integrated by :func:`integrate_region` against the chord between its own two
    boundaries instead. See ``_cluster_spans`` for why the distinction is restricted to
    present analytes, and ``IntegrationParams.overlap_frac`` for the threshold.

    The baseline is the straight line joining the cluster's outer base points
    (``base_l``..``base_r``, at the chromatographic floor), evaluated across this analyte's
    [i_start, i_end] split range. Splitting overlapping peaks this way, rather than drawing
    the baseline up to the elevated inter-peak valley, recovers a shoulder peak's full area
    (e.g. congener 5 on the side of congener 8). Returns (area, baseline@start, baseline@end).
    """
    if i_end <= i_start:
        return 0.0, float(intensity[i_start]), float(intensity[i_start])
    xl, xr = float(time[base_l]), float(time[base_r])
    yl, yr = float(intensity[base_l]), float(intensity[base_r])
    span = (xr - xl) or 1.0
    t = time[i_start:i_end + 1]
    baseline = yl + (yr - yl) * (t - xl) / span
    corrected = np.clip(intensity[i_start:i_end + 1] - baseline, 0.0, None)
    area = float(np.trapezoid(corrected, t * 60.0))  # counts·seconds (Agilent convention)
    return area, float(baseline[0]), float(baseline[-1])



def curvature(signal: np.ndarray, span: int = 2) -> np.ndarray:
    """Discrete second derivative over a ±``span`` stencil.

    A peak riding on the flank of a much larger neighbour produces no local maximum. The
    trace is still descending through it, so topographic prominence cannot see it. It does
    produce downward curvature: the second derivative dips to a local minimum at the rider's
    apex and returns through zero at the inflections either side. That is how a
    chromatography data system finds a shoulder, and it is what an analyst is reacting to by
    eye when they skim one.

    The wider stencil (``span`` samples rather than adjacent points) avoids the
    sample-to-sample noise that dominates a three-point second difference on real data.
    """
    out = np.zeros_like(signal, dtype=float)
    if signal.size < 2 * span + 1:
        return out
    out[span:-span] = signal[2 * span:] - 2.0 * signal[span:-span] + signal[: -2 * span]
    return out


def find_shoulder(
    detect: np.ndarray,
    lo: int,
    hi: int,
    min_curvature: float,
) -> tuple[int, int, int] | None:
    """Locate a rider inside ``[lo, hi]`` as ``(i_start, i_apex, i_end)``, or None.

    The apex is the strongest downward curvature in the window. The bounds are where that
    curvature relaxes to zero on either side, which is where the rider merges back into its
    parent's flank. Those inflection points are what a tangent skim draws its baseline
    between. See :func:`integrate_skim`.
    """
    if hi - lo < 4:
        return None
    d2 = curvature(detect)
    window = d2[lo:hi + 1]
    if window.size == 0:
        return None
    apex = lo + int(np.argmin(window))
    if -d2[apex] < min_curvature:
        return None
    left = apex
    while left > lo and d2[left] < 0:
        left -= 1
    right = apex
    while right < hi and d2[right] < 0:
        right += 1
    if right - left < 3 or apex in (left, right):
        return None
    return left, apex, right


def integrate_skim(
    times: np.ndarray,
    signal: np.ndarray,
    i_start: int,
    i_end: int,
    level: float,
) -> tuple[float, float, float]:
    """Area of a rider above a *flat* baseline at ``level``, in counts·seconds.

    Why flat, and not the chord ``integrate_region`` draws. A rider on a leading edge sits on
    a trace that climbs monotonically into its parent. It has no valley on the parent side
    and no local maximum of its own. A chord between its two boundaries would run from the
    rider's foot up to a point far up the parent's flank, and because that flank is convex
    the chord lies *above* the trace almost everywhere between them, so the clipped area
    comes out to zero. For congeners 83, 51 and 115, a chord baseline reports 0.0 while the
    analyst is integrating tens of thousands of counts·s.

    The flat line at the rider's own foot reflects what the analyst draws. Everything above
    the level the trace had *before* the rider began belongs to the rider, up to where its
    parent's rise takes over (the caller's ``i_end``).
    """
    if i_end <= i_start:
        return 0.0, level, level
    x = times[i_start:i_end + 1]
    y = signal[i_start:i_end + 1].astype(float)
    above = np.clip(y - level, 0.0, None)
    return float(np.trapezoid(above, x) * 60.0), level, level

def estimate_noise(intensity: np.ndarray, i_start: int, i_end: int) -> float:
    """Rough local noise estimate from the baseline shoulders around a peak.

    RETAINED-UNUSED: the integrator uses ``global_noise`` instead, so that the presence
    threshold is one number per trace rather than a per-peak quantity that shifts as
    boundaries move. Kept for diagnostics and for a future local-noise LOQ.
    """
    n = intensity.size
    lo = max(0, i_start - 10)
    hi = min(n, i_end + 11)
    shoulders = np.concatenate([intensity[lo:i_start], intensity[i_end + 1:hi]])
    if shoulders.size >= 3:
        return float(np.std(shoulders)) or 1.0
    return 1.0


def global_noise(intensity: np.ndarray) -> float:
    """Trace-wide noise estimate (MAD of the lower-intensity samples)."""
    if intensity.size < 5:
        return 1.0
    # use the bottom 60% of samples (baseline regions) for a spread estimate
    lo = intensity[intensity <= np.quantile(intensity, 0.6)]
    if lo.size < 3:
        lo = intensity
    mad = float(np.median(np.abs(lo - np.median(lo))))
    return (mad * 1.4826) or 1.0
