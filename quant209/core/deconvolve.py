"""Resolve overlapping peaks by fitting their shapes, for riders no boundary can separate.

A congener riding the *leading* edge of a much larger neighbour cannot be separated by any
boundary rule. The trace climbs monotonically through it, so wherever a perpendicular falls
it either takes some of the parent or leaves some of the rider. There is no valley to drop
at and no local maximum to bound. What separates them is *shape*: two peaks of known form
summing to the observed envelope.

The model is a sum of Gaussians on a flat offset, fitted by Levenberg-Marquardt implemented
here in numpy (scipy is deliberately absent from this project's hot path. See METHODS).
Three constraints do the real work, and each is a statement about the chromatography rather
than a convenience:

* **One shared width.** Peaks a few hundredths of a minute apart came off the same column
  under the same conditions, so they have the same width. Letting each component choose its
  own is what lets a rider widen until it explains its parent's entire rising edge. On this
  archive an unconstrained fit returns up to 7 times the analyst's area.
* **Centres pinned near their anchors.** Where each congener elutes was already decided, from
  the run's own calibration injection. The fit is asked "how much of this envelope belongs to
  each", not "where are there peaks".
* **A fit-quality gate.** A decomposition that does not reproduce the trace is not evidence
  about anything. A deeply buried rider is genuinely ill-posed, since two components can
  trade amplitude and fit equally well. :func:`deconvolve_rider` therefore returns ``None``
  rather than a number whenever the residual, the amplitude or the width look implausible,
  and the caller keeps whatever it had.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

SQRT_2PI = float(np.sqrt(2.0 * np.pi))


@dataclass(frozen=True)
class RiderFit:
    """What an accepted fit says about the rider: how much, how well, how tall, and where.

    ``height`` and ``centre`` exist for the caller that has no apex of its own to measure.
    A rider with no local maximum and no visible bulge has nowhere on the trace to read a
    height from. The trace at any one sample there is mostly its parent. The fitted
    component's own amplitude is the only height that describes the rider.
    """

    area: float  # counts·seconds, the fitted component's analytic area
    residual: float  # ||y - model|| / ||y|| over the fitted window
    height: float  # counts, the fitted component's amplitude above the fit's flat offset
    centre: float  # minutes, the fitted component's retention time

#: Largest relative residual (||y - model|| / ||y||) an accepted fit may leave behind.
#: Above this the model is not describing the trace and its apportionment means nothing.
MAX_RESIDUAL = 0.25

#: How far a component's centre may move from its anchor, as a fraction of the peak width.
CENTRE_TOL_SIGMAS = 0.4

#: Bounds on the shared width, as a fraction of the width measured on the parent peak.
SIGMA_BOUNDS = (0.7, 1.4)

#: Half-width of the fitted window, in peak widths either side of the outermost component.
WINDOW_SIGMAS = 6.0

#: A rider whose fitted area exceeds this fraction of the whole fitted envelope is rejected:
#: a rider that is most of its own parent is a mis-identification, not a shoulder.
MAX_RIDER_SHARE = 0.5

_MAX_ITERATIONS = 40


def _model(x: np.ndarray, params: np.ndarray, n: int) -> np.ndarray:
    """Sum of ``n`` Gaussians sharing one width: ``[sigma, amp_0, c_0, amp_1, c_1, ...]``."""
    sigma = max(abs(params[0]), 1e-6)
    out = np.zeros_like(x, dtype=float)
    for k in range(n):
        out += params[1 + 2 * k] * np.exp(-0.5 * ((x - params[2 + 2 * k]) / sigma) ** 2)
    return out


def peak_sigma(times: np.ndarray, detect: np.ndarray, apex: int) -> float:
    """Gaussian width of the peak at ``apex``, from its full width at half maximum."""
    if not 0 <= apex < detect.size:
        return 0.03
    half = float(detect[apex]) / 2.0
    left = apex
    while left > 0 and detect[left] > half:
        left -= 1
    right = apex
    while right < detect.size - 1 and detect[right] > half:
        right += 1
    fwhm = float(times[right] - times[left])
    return max(fwhm / 2.355, 1e-3) if fwhm > 0 else 0.03


def fit_gaussians(
    x: np.ndarray,
    y: np.ndarray,
    centres: list[float],
    sigma0: float,
    sigma_bounds: tuple[float, float],
    centre_tol: float,
) -> tuple[np.ndarray, float]:
    """Levenberg-Marquardt fit of shared-width Gaussians. Returns ``(params, residual)``."""
    n = len(centres)
    params = np.zeros(1 + 2 * n, dtype=float)
    params[0] = sigma0
    for k, centre in enumerate(centres):
        params[1 + 2 * k] = max(float(y[int(np.argmin(np.abs(x - centre)))]), 1.0)
        params[2 + 2 * k] = centre
    lo, hi = params.copy(), params.copy()
    lo[0], hi[0] = sigma_bounds
    ceiling = max(float(y.max()) * 3.0, 1.0)
    for k in range(n):
        lo[1 + 2 * k], hi[1 + 2 * k] = 0.0, ceiling
        lo[2 + 2 * k] = centres[k] - centre_tol
        hi[2 + 2 * k] = centres[k] + centre_tol

    lam = 1e-2
    scale = max(float(np.linalg.norm(y)), 1.0)
    for _ in range(_MAX_ITERATIONS):
        model = _model(x, params, n)
        residual = y - model
        cost = float(residual @ residual)
        jac = np.zeros((x.size, params.size), dtype=float)
        for j in range(params.size):
            step = max(abs(params[j]) * 1e-4, 1e-7)
            bumped = params.copy()
            bumped[j] += step
            jac[:, j] = (_model(x, bumped, n) - model) / step
        jtj = jac.T @ jac
        jtr = jac.T @ residual
        for _attempt in range(6):
            damping = lam * np.diag(np.diag(jtj) + 1e-12)
            try:
                delta = np.linalg.solve(jtj + damping, jtr)
            except np.linalg.LinAlgError:
                lam *= 10.0
                continue
            trial = np.clip(params + delta, lo, hi)
            trial_residual = y - _model(x, trial, n)
            if float(trial_residual @ trial_residual) < cost:
                params, lam = trial, max(lam * 0.3, 1e-9)
                break
            lam *= 10.0
        else:
            break
    residual = y - _model(x, params, n)
    return params, float(np.linalg.norm(residual) / scale)


def component_area(params: np.ndarray, k: int) -> float:
    """Analytic area of component ``k``, in counts·seconds (``x`` is in minutes)."""
    return float(params[1 + 2 * k] * abs(params[0]) * SQRT_2PI * 60.0)


def deconvolve_rider(
    times: np.ndarray,
    signal: np.ndarray,
    detect: np.ndarray,
    rider_rt: float,
    neighbour_rts: list[float],
    parent_apex: int,
) -> RiderFit | None:
    """Fit the rider at ``rider_rt``, or return None if the fit cannot be trusted.

    ``neighbour_rts`` are the retention times of the resolved peaks sharing the envelope:
    the parent and anything else close enough to overlap it.
    """
    sigma0 = peak_sigma(times, detect, parent_apex)
    centres = sorted([rider_rt, *neighbour_rts])
    rider_index = centres.index(rider_rt)
    lo_t = centres[0] - WINDOW_SIGMAS * sigma0
    hi_t = centres[-1] + WINDOW_SIGMAS * sigma0
    window = np.where((times >= lo_t) & (times <= hi_t))[0]
    if window.size < 2 * len(centres) + 6:
        return None

    x = times[window]
    y = signal[window].astype(float)
    offset = float(np.percentile(y, 5))
    params, residual = fit_gaussians(
        x,
        y - offset,
        centres=centres,
        sigma0=sigma0,
        sigma_bounds=(SIGMA_BOUNDS[0] * sigma0, SIGMA_BOUNDS[1] * sigma0),
        centre_tol=max(CENTRE_TOL_SIGMAS * sigma0, 0.008),
    )
    if residual > MAX_RESIDUAL:
        return None
    area = component_area(params, rider_index)
    total = sum(component_area(params, k) for k in range(len(centres)))
    if area <= 0.0 or total <= 0.0 or area / total > MAX_RIDER_SHARE:
        return None
    return RiderFit(
        area=area,
        residual=residual,
        height=float(params[1 + 2 * rider_index]),
        centre=float(params[2 + 2 * rider_index]),
    )


@dataclass(frozen=True)
class PairFit:
    """Two comparable congeners' shares of one merged envelope."""

    first: RiderFit
    second: RiderFit
    drift: float  # minutes both centres moved from their anchors, together
    residual: float


#: The least share of a merged envelope either congener must hold for a pair fit to stand. A
#: component smaller than this is the fit spending a parameter on a shoulder of noise.
MIN_PAIR_SHARE = 0.05


def fit_pair(
    times: np.ndarray,
    signal: np.ndarray,
    first_rt: float,
    second_rt: float,
    max_drift: float,
    window: tuple[int, int] | None = None,
) -> PairFit | None:
    """Split the envelope two co-eluting congeners merged into, or None.

    The model is rigid on purpose. The two congeners keep the spacing the run's own
    calibration injection measured between them, and move together by one drift: in a
    concentrated sample the whole window drifts, it does not rearrange. They share one width,
    for the reason ``deconvolve_rider`` gives. The amplitudes are then a linear least-squares
    problem, solved on a grid of drift and width. A fit that does not reproduce the envelope
    (``MAX_RESIDUAL``), whose drift ran to the edge of what was allowed, or that leaves either
    congener under ``MIN_PAIR_SHARE`` of the envelope, says nothing about the split.

    ``window`` restricts the fit to one envelope (sample indices ``start``, ``end``), measured
    above its own baseline chord. Without it the fit spans 0.6 min either side of the pair.
    """
    if window is None:
        lo_t, hi_t = min(first_rt, second_rt) - 0.6, max(first_rt, second_rt) + 0.6
        mask = (times >= lo_t) & (times <= hi_t)
        if mask.sum() < 12:
            return None
        x = times[mask]
        y = signal[mask].astype(float)
        y = y - float(np.percentile(y, 5))
    else:
        # The envelope alone, samples ``start`` to ``end``, above the chord between them: a
        # neighbouring peak inside the default window would otherwise defeat the fit.
        start, end = window
        if end - start < 12:
            return None
        x = times[start:end + 1]
        y = signal[start:end + 1].astype(float)
        y = y - np.interp(x, [x[0], x[-1]], [y[0], y[-1]])
    scale = max(float(np.linalg.norm(y)), 1.0)
    drifts = np.arange(-max_drift, max_drift + 1e-9, 0.005)
    best = None
    for sigma in np.arange(0.02, 0.1201, 0.005):
        for drift in drifts:
            g1 = np.exp(-0.5 * ((x - (first_rt + drift)) / sigma) ** 2)
            g2 = np.exp(-0.5 * ((x - (second_rt + drift)) / sigma) ** 2)
            design = np.vstack([g1, g2]).T
            amps, *_ = np.linalg.lstsq(design, y, rcond=None)
            if (amps <= 0).any():
                continue
            residual = float(np.linalg.norm(y - design @ amps)) / scale
            if best is None or residual < best[0]:
                best = (residual, float(drift), float(sigma), amps)
    if best is None:
        return None
    residual, drift, sigma, amps = best
    if residual > MAX_RESIDUAL or abs(abs(drift) - max_drift) < 0.0051:
        return None
    k = sigma * SQRT_2PI * 60.0
    areas = amps * k
    if areas.min() < MIN_PAIR_SHARE * areas.sum():
        return None
    return PairFit(
        first=RiderFit(area=float(areas[0]), residual=residual, height=float(amps[0]), centre=first_rt + drift),
        second=RiderFit(area=float(areas[1]), residual=residual, height=float(amps[1]), centre=second_rt + drift),
        drift=drift,
        residual=residual,
    )


#: Spacing a free pair fit may give two congeners (minutes). Closer than the first, two
#: components of one width are one peak and the split between them means nothing. Farther
#: than the second, the pair would have shown two apexes and not needed a fit.
FREE_PAIR_SPACING = (0.015, 0.15)

#: A free pair fit stands only if it reproduces the envelope this much better than one peak
#: does (ratio of relative residuals). One peak is what a congener alone looks like, and two
#: components of one width can always imitate it as two halves: that is not evidence of two.
FREE_PAIR_GAIN = 0.5


def _gauss_grid(x: np.ndarray, centres: np.ndarray, sigma: float) -> np.ndarray:
    return np.exp(-0.5 * ((x[None, :] - centres[:, None]) / sigma) ** 2)


def fit_unresolved_pair(
    times: np.ndarray,
    signal: np.ndarray,
    first_rt: float,
    second_rt: float,
    start: int,
    end: int,
) -> PairFit | None:
    """Split the envelope of two congeners the run's calibration could not resolve, or None.

    Unlike :func:`fit_pair`, the two centres are free. The pair's spacing cannot be taken from
    its anchors: the congener that had no apex in the calibration was anchored at a predicted
    retention time, not a measured one. The fit covers the envelope itself, samples ``start``
    to ``end`` above the chord between them, and each centre may lie anywhere inside it, the
    earlier congener first, ``FREE_PAIR_SPACING`` apart, with one shared width. The
    amplitudes are a linear least-squares problem, solved on the grid in closed form.

    The fit stands only if it reproduces the envelope (``MAX_RESIDUAL``), clearly better than
    one peak does (``FREE_PAIR_GAIN``), with neither centre at the edge of its range and each
    congener holding at least ``MIN_PAIR_SHARE``. ``first`` is the component of the congener
    at ``first_rt``, whichever elutes first.
    """
    if end - start < 12:
        return None
    x = times[start:end + 1]
    y = signal[start:end + 1].astype(float)
    y = y - np.interp(x, [x[0], x[-1]], [y[0], y[-1]])
    lo_c, hi_c = float(x[0]), float(x[-1])
    yy = float(y @ y)
    if yy <= 0:
        return None
    scale = np.sqrt(yy)
    centres = np.arange(lo_c, hi_c + 1e-9, 0.0025)
    gaps = np.arange(FREE_PAIR_SPACING[0], FREE_PAIR_SPACING[1] + 1e-9, 0.0025)
    best_one, best_two = np.inf, None
    for sigma in np.arange(0.015, 0.1201, 0.0025):
        g = _gauss_grid(x, centres, sigma)
        gg = np.einsum("ij,ij->i", g, g)
        gy = g @ y
        # One peak: amplitude gy/gg, residual sum of squares yy - gy^2/gg.
        one = np.sqrt(np.maximum(yy - gy ** 2 / gg, 0.0)) / scale
        best_one = min(best_one, float(one[gy > 0].min(initial=np.inf)))
        for gap in gaps:
            k = int(round(gap / 0.0025))
            if k >= centres.size:
                break
            g1, g2 = g[:-k], g[k:]
            a11, a22 = gg[:-k], gg[k:]
            a12 = np.einsum("ij,ij->i", g1, g2)
            b1, b2 = gy[:-k], gy[k:]
            det = a11 * a22 - a12 ** 2
            ok = det > 1e-12 * a11 * a22
            amp1 = np.where(ok, (a22 * b1 - a12 * b2) / np.where(ok, det, 1.0), -1.0)
            amp2 = np.where(ok, (a11 * b2 - a12 * b1) / np.where(ok, det, 1.0), -1.0)
            rss = yy - (amp1 * b1 + amp2 * b2)
            rss = np.where((amp1 > 0) & (amp2 > 0), rss, np.inf)
            j = int(np.argmin(rss))
            if np.isfinite(rss[j]) and (best_two is None or rss[j] < best_two[0]):
                best_two = (float(rss[j]), float(centres[j]), float(centres[j + k]), float(sigma),
                            float(amp1[j]), float(amp2[j]))
    if best_two is None:
        return None
    rss, c1, c2, sigma, h1, h2 = best_two
    residual = float(np.sqrt(max(rss, 0.0))) / scale
    if residual > MAX_RESIDUAL or residual > FREE_PAIR_GAIN * best_one:
        return None
    if c1 - lo_c < 0.0026 or hi_c - c2 < 0.0026:
        return None  # a centre ran to the edge of its range: the fit was looking elsewhere
    a1, a2 = h1 * sigma * SQRT_2PI * 60.0, h2 * sigma * SQRT_2PI * 60.0
    if min(a1, a2) < MIN_PAIR_SHARE * (a1 + a2):
        return None
    earlier = RiderFit(area=a1, residual=residual, height=h1, centre=c1)
    later = RiderFit(area=a2, residual=residual, height=h2, centre=c2)
    first, second = (earlier, later) if first_rt <= second_rt else (later, earlier)
    return PairFit(first=first, second=second, drift=0.0, residual=residual)
