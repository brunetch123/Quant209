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

import numpy as np

SQRT_2PI = float(np.sqrt(2.0 * np.pi))

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
) -> tuple[float, float] | None:
    """Fitted area of the rider at ``rider_rt``, or None if the fit cannot be trusted.

    ``neighbour_rts`` are the retention times of the resolved peaks sharing the envelope:
    the parent and anything else close enough to overlap it. Returns
    ``(area, relative_residual)``.
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
    return area, residual
