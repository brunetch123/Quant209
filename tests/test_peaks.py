"""Tests for quant209.core.peaks, integration primitives (apex, prominence, boundaries,
baseline). Traces are synthetic (numpy-built Gaussians / noise), asserted against closed-form
analytic values rather than against the implementation's own intermediate outputs, since
these primitives are the numerical foundation the RT-anchored integration engine builds on.
"""

import math

import numpy as np
import pytest

from quant209.core.peaks import (
    IntegrationParams,
    PeakResult,
    _local_min,  # private helper, tested directly here (primitives module, no public wrapper)
    estimate_noise,
    find_apex,
    global_noise,
    hill_climb,
    integrate_above_baseline,
    integrate_region,
    peak_prominences,
    smooth,
    valley_index,
)

DT = 0.005  # minutes, a representative MRM sampling interval


def _time_axis(t_min: float, t_max: float, dt: float = DT) -> np.ndarray:
    n = int(round((t_max - t_min) / dt)) + 1
    return t_min + dt * np.arange(n)


def _gaussian(t: np.ndarray, height: float, mu: float, sigma: float) -> np.ndarray:
    return height * np.exp(-0.5 * ((t - mu) / sigma) ** 2)


def _gaussian_area_counts_seconds(height: float, sigma: float) -> float:
    """Analytic area of an isolated Gaussian peak (RT in minutes) in counts·seconds."""
    return height * sigma * math.sqrt(2 * math.pi) * 60.0


# ---------------------------------------------------------------------------
# IntegrationParams / PeakResult
# ---------------------------------------------------------------------------


def test_integration_params_legacy_defaults():
    """Defaults must match the analyst-calibrated values exactly. Reproducing analyst manual
    integration depends on it, so changing any of these would silently retune it.

    ``overlap_frac`` is the one exception and is pinned separately below: it was never one of
    the analyst-calibrated values, because until the cluster baseline shipped no code path
    read it at all.
    """
    params = IntegrationParams()
    assert params.rt_window == 0.10
    assert params.max_reach == 0.20
    assert params.max_halfwidth == 0.25
    assert params.smooth_window == 5
    assert params.smooth_poly == 2
    assert params.min_sn == 3.0
    assert params.presence_sn == 4.0
    assert params.min_prominence == 150.0
    assert params.integrate_on_smoothed is False


def test_overlap_frac_default_is_the_measured_cluster_threshold():
    """``core.integrate._cluster_spans`` reads ``overlap_frac``, and the default is the
    middle of the 0.02-0.15 plateau over which the same pairs merge on every calibration
    injection of both reference runs.

    Moving it changes how much of a crowded calibration window's area is recovered, and
    therefore every RRF in a run, so the value is pinned to the measured plateau rather than
    chosen arbitrarily.
    """
    assert IntegrationParams().overlap_frac == 0.05


def test_the_named_rules_name_four_congeners_and_no_more():
    """The two label lists are the only rules in the engine that read a label, not a trace.

    Pinned so that growing either list is a deliberate, reviewed act: each entry is a claim
    about the chromatography of one method, measured against the archive
    (``docs/METHODS.md`` section 2.7), and a label added casually would be integrated by
    rule in every file where the trace shows nothing.
    """
    params = IntegrationParams()
    assert params.forced_riders == ("83", "117")
    assert params.forced_present == ("90+101+113", "118")
    # Immutable defaults: one instance emptying a list must not empty it for the next.
    assert IntegrationParams(forced_riders=()).forced_riders == ()
    assert IntegrationParams().forced_riders == ("83", "117")


def test_a_displaced_peak_needs_a_tenfold_abundance_advantage():
    """An order of magnitude in the published Aroclor profiles before the trace's ambiguity is
    settled by composition. Pinned so that loosening it is a deliberate act."""
    assert IntegrationParams().dominance_ratio == 10.0


def test_peak_result_default_eq_and_flags():
    """No array fields, so default dataclass eq should work (unlike Chromatogram, which needs
    eq=False for its ndarray fields), and the flags list defaults to its own empty list per
    instance, not a shared mutable default."""
    kwargs = dict(
        label="8", q1=291.9, q3=222.0, rt_apex=12.3, rt_start=12.1, rt_end=12.5,
        area=1000.0, height=500.0, baseline_start=10.0, baseline_end=12.0,
        sn=25.0, method="valley_drop",
    )
    a = PeakResult(**kwargs)
    b = PeakResult(**kwargs)
    assert a == b
    assert a.flags == [] and b.flags == []
    a.flags.append("low_sn")
    assert b.flags == [], "default_factory must give each instance its own list"
    assert a.i_start == -1 and a.i_apex == -1 and a.i_end == -1
    # "not recorded" sentinels, not NaN. NaN would break the default __eq__ asserted above.
    assert a.src_q1 == -1.0 and a.src_q3 == -1.0


def test_peak_result_detected_tracks_the_absent_flag():
    """``detected`` is the presence gate downstream code must use instead of ``sn``.

    An absent analyte is integrated across its RT territory rather than around its apex, so
    its ``sn`` can be arbitrarily large. Only the flag carries the detection verdict.
    """
    kwargs = dict(
        label="8", q1=291.9, q3=222.0, rt_apex=12.3, rt_start=12.1, rt_end=12.5,
        area=1000.0, height=500.0, baseline_start=10.0, baseline_end=12.0,
        sn=25.0, method="rt_window",
    )
    assert PeakResult(**kwargs).detected is True
    assert PeakResult(**kwargs, flags=["low_sn"]).detected is True
    assert PeakResult(**kwargs, flags=["absent"]).detected is False
    assert PeakResult(**{**kwargs, "sn": 4.9e4}, flags=["absent", "narrow"]).detected is False


# ---------------------------------------------------------------------------
# Isolated Gaussian: apex, boundaries, area
# ---------------------------------------------------------------------------


def test_isolated_gaussian_area_matches_analytic():
    height, mu, sigma = 1e5, 0.50, 0.02
    time = _time_axis(0.0, 1.0)
    signal = _gaussian(time, height, mu, sigma)
    lo, hi = 0, time.size - 1

    i_apex = find_apex(time, signal, rt_expected=mu, window=0.10, lo=lo, hi=hi)
    i_apex = hill_climb(signal, i_apex, lo, hi)
    assert time[i_apex] == pytest.approx(mu)

    half_reach = int(round(0.20 / DT))  # IntegrationParams.max_reach, in samples
    i_start = _local_min(signal, max(lo, i_apex - half_reach), i_apex, lo)
    i_end = _local_min(signal, i_apex, min(hi, i_apex + half_reach), hi)

    area, base_l, base_r = integrate_region(time, signal, i_start, i_end)

    expected = _gaussian_area_counts_seconds(height, sigma)
    assert area == pytest.approx(expected, rel=0.02)
    # boundaries are ~10 sigma out: baseline should be indistinguishable from zero
    assert base_l == pytest.approx(0.0, abs=1e-3)
    assert base_r == pytest.approx(0.0, abs=1e-3)


def test_isolated_gaussian_prominence_exceeds_min_prominence():
    """A real, well-resolved peak's prominence must clear the presence floor by a wide margin
    (the flip side of the flat-noise case below, which must fall short of it)."""
    height, mu, sigma = 1e5, 0.50, 0.02
    time = _time_axis(0.0, 1.0)
    signal = _gaussian(time, height, mu, sigma)
    i_apex = find_apex(time, signal, rt_expected=mu, window=0.10, lo=0, hi=time.size - 1)

    prom = peak_prominences(signal, [i_apex])[0]

    assert prom > IntegrationParams().min_prominence
    assert prom == pytest.approx(height, rel=1e-6)  # descends essentially to zero on both sides


def test_peak_prominences_simple_triangle():
    """Direct correctness check of the topographic-prominence algorithm on a hand-built trace,
    independent of any Gaussian machinery."""
    signal = np.array([0.0, 0.0, 0.0, 5.0, 10.0, 5.0, 0.0, 0.0, 0.0])
    proms = peak_prominences(signal, [4])
    assert proms[0] == pytest.approx(10.0)

    # a smaller secondary bump nested between two higher shoulders: prominence is limited by
    # the *lower* of the two directions it can descend before the trace rises back above it.
    signal2 = np.array([0.0, 10.0, 3.0, 6.0, 3.0, 10.0, 0.0])
    proms2 = peak_prominences(signal2, [3])
    assert proms2[0] == pytest.approx(3.0)  # 6 - max(3, 3)


# ---------------------------------------------------------------------------
# Two overlapping Gaussians: valley boundary + split-baseline areas
# ---------------------------------------------------------------------------


def test_two_overlapping_gaussians_valley_between_apexes_and_areas_match():
    h1, mu1 = 1.0e5, 0.40
    h2, mu2 = 0.9e5, 0.48  # 0.08 min apart, comparable (not identical) heights
    sigma = 0.02
    time = _time_axis(0.0, 1.0)
    signal = _gaussian(time, h1, mu1, sigma) + _gaussian(time, h2, mu2, sigma)

    mid_idx = int(round(((mu1 + mu2) / 2 - 0.0) / DT))
    i_apex1 = hill_climb(
        signal, find_apex(time, signal, mu1, 0.10, 0, mid_idx), 0, mid_idx
    )
    i_apex2 = hill_climb(
        signal, find_apex(time, signal, mu2, 0.10, mid_idx, time.size - 1), mid_idx, time.size - 1
    )
    assert time[i_apex1] == pytest.approx(mu1)
    assert time[i_apex2] == pytest.approx(mu2)

    i_valley = valley_index(signal, i_apex1, i_apex2)
    assert i_apex1 < i_valley < i_apex2

    half_reach = int(round(0.20 / DT))
    i_start1 = _local_min(signal, max(0, i_apex1 - half_reach), i_apex1, 0)
    i_end2 = _local_min(signal, i_apex2, min(time.size - 1, i_apex2 + half_reach), time.size - 1)

    # Split the merged cluster at the valley, baselined against the cluster's own outer
    # (near-zero) base points rather than the elevated valley. integrate_above_baseline is
    # exactly the primitive documented for this ("recovers a shoulder peak's full area").
    area1, _, _ = integrate_above_baseline(time, signal, i_start1, i_valley, i_start1, i_end2)
    area2, _, _ = integrate_above_baseline(time, signal, i_valley, i_end2, i_start1, i_end2)

    expected1 = _gaussian_area_counts_seconds(h1, sigma)
    expected2 = _gaussian_area_counts_seconds(h2, sigma)
    assert area1 == pytest.approx(expected1, rel=0.05)
    assert area2 == pytest.approx(expected2, rel=0.05)


# ---------------------------------------------------------------------------
# Elevated flat baseline
# ---------------------------------------------------------------------------


def test_elevated_baseline_area_unchanged():
    height, mu, sigma = 1e5, 0.50, 0.02
    baseline_level = 500.0
    time = _time_axis(0.0, 1.0)
    signal = baseline_level + _gaussian(time, height, mu, sigma)
    lo, hi = 0, time.size - 1

    i_apex = find_apex(time, signal, rt_expected=mu, window=0.10, lo=lo, hi=hi)
    i_apex = hill_climb(signal, i_apex, lo, hi)

    half_reach = int(round(0.20 / DT))
    i_start = _local_min(signal, max(lo, i_apex - half_reach), i_apex, lo)
    i_end = _local_min(signal, i_apex, min(hi, i_apex + half_reach), hi)

    area, base_l, base_r = integrate_region(time, signal, i_start, i_end)

    assert base_l == pytest.approx(baseline_level, abs=1e-3)
    assert base_r == pytest.approx(baseline_level, abs=1e-3)
    expected = _gaussian_area_counts_seconds(height, sigma)
    assert area == pytest.approx(expected, rel=0.02)


# ---------------------------------------------------------------------------
# Flat noise only: absent analyte should not integrate to a real peak's area, and its
# S/N should fall short of the presence thresholds.
# ---------------------------------------------------------------------------


def test_flat_noise_absent_analyte_area_and_sn_below_thresholds():
    rt_expected = 0.5
    time = _time_axis(0.0, 1.0)
    # Seed 0 is hand-picked, not representative: it verifiably clears both thresholds by a
    # comfortable margin (measured sn ~2.16). Across candidate seeds, sn lands close to
    # min_sn=3.0 often enough to matter, with roughly 3/50 falling in [3.0, 3.45], so an
    # arbitrary seed would flake on the `sn < min_sn` assertion below about 6% of the time.
    rng = np.random.default_rng(0)
    noise = rng.normal(0.0, 50.0, size=time.size)

    max_reach = IntegrationParams().max_reach
    i_start = int(round((rt_expected - max_reach) / DT))
    i_end = int(round((rt_expected + max_reach) / DT))

    area, _, _ = integrate_region(time, noise, i_start, i_end)
    real_peak_area = _gaussian_area_counts_seconds(1e5, 0.02)
    assert abs(area) < 0.05 * real_peak_area

    i_apex = find_apex(time, noise, rt_expected, IntegrationParams().rt_window, i_start, i_end)
    noise_est = estimate_noise(noise, i_start, i_end)
    sn = noise[i_apex] / noise_est

    assert sn < IntegrationParams().min_sn
    assert sn < IntegrationParams().presence_sn


# ---------------------------------------------------------------------------
# smooth()
# ---------------------------------------------------------------------------


def test_smooth_preserves_linear_ramp_in_interior():
    """A Savitzky-Golay filter of any polynomial order ≥ 1 exactly reproduces a linear input,
    but only where the sliding window sees genuine ramp values. Edge padding (mode="edge")
    replaces the would-be continuation with a repeated boundary value, so exactness is checked
    away from the first/last `window // 2` samples."""
    window, poly = IntegrationParams().smooth_window, IntegrationParams().smooth_poly
    x = np.arange(200) * DT
    ramp = 3.0 + 7.0 * x

    smoothed = smooth(ramp, window, poly)

    half = window // 2
    np.testing.assert_allclose(smoothed[half:-half], ramp[half:-half], atol=1e-9)


def test_smooth_reduces_noise_variance():
    window, poly = IntegrationParams().smooth_window, IntegrationParams().smooth_poly
    rng = np.random.default_rng(7)
    noisy = 1000.0 + rng.normal(0.0, 80.0, size=200)

    smoothed = smooth(noisy, window, poly)

    assert np.var(smoothed) < 0.6 * np.var(noisy)


def test_smooth_gaussian_apex_curvature_fidelity():
    """The linear-ramp test above is satisfied by ANY zero-bias linear filter (e.g. a same-
    width boxcar average), since a ramp has zero curvature. It can't catch a wrong-kernel
    regression. Near a curved apex the two diverge sharply: true SG(5,2) tracks the local
    quadratic curvature closely (~-0.095% measured bias here) while a boxcar average washes
    the peak down (~-5.9% measured bias), so this pins the actual SG property."""
    window, poly = IntegrationParams().smooth_window, IntegrationParams().smooth_poly
    mu, sigma, height = 0.50, 0.02, 1e5
    t = mu + DT * np.arange(-20, 21)
    signal = _gaussian(t, height, mu, sigma)
    i_apex = 20  # center sample, exactly at mu

    smoothed = smooth(signal, window, poly)
    sg_bias = abs(smoothed[i_apex] - signal[i_apex]) / signal[i_apex]
    assert sg_bias < 0.005

    # Same-width boxcar average, built inline (not imported from peaks.py) so this test
    # doesn't depend on the module having a boxcar implementation. It's only a foil to
    # prove the assertion above actually discriminates a correct SG kernel from a wrong one.
    boxcar_weights = np.ones(window) / window
    padded = np.pad(signal, window // 2, mode="edge")
    boxcar = np.convolve(padded, boxcar_weights, mode="valid")
    boxcar_bias = abs(boxcar[i_apex] - signal[i_apex]) / signal[i_apex]
    assert boxcar_bias > 0.01  # comfortably clears the 0.5% SG tolerance above


def test_smooth_leaves_signal_untouched_when_disabled():
    """window < 3 (e.g. 0, meaning "off" per IntegrationParams docstring) returns the trace
    unmodified (as a float array), and the input itself is never mutated in place. The same
    trace gets re-integrated for neighbouring analytes."""
    y = np.array([1.0, 5.0, 2.0, 9.0, 3.0])
    y_copy = y.copy()

    out = smooth(y, 0, 2)

    np.testing.assert_array_equal(out, y_copy)
    np.testing.assert_array_equal(y, y_copy)  # not mutated


# ---------------------------------------------------------------------------
# hill_climb
# ---------------------------------------------------------------------------


def test_hill_climb_reaches_apex_from_shoulder():
    height, mu, sigma = 1e5, 0.50, 0.02
    time = _time_axis(0.0, 1.0)
    signal = _gaussian(time, height, mu, sigma)
    i_apex_true = int(round(mu / DT))
    i_shoulder = i_apex_true - 5  # on the rising side, below the apex

    result = hill_climb(signal, i_shoulder, 0, signal.size - 1)

    assert result == i_apex_true


def test_hill_climb_respects_bounds():
    height, mu, sigma = 1e5, 0.50, 0.02
    time = _time_axis(0.0, 1.0)
    signal = _gaussian(time, height, mu, sigma)
    i_apex_true = int(round(mu / DT))
    i_shoulder = i_apex_true - 5
    hi_cap = i_apex_true - 3  # wall placed before the true apex

    result = hill_climb(signal, i_shoulder, 0, hi_cap)

    assert result == hi_cap  # climbed up to (not past) the bound


# ---------------------------------------------------------------------------
# global_noise robustness
# ---------------------------------------------------------------------------


def test_global_noise_robust_to_outlier_spikes():
    """global_noise uses the bottom 60% of samples (a robust MAD estimate), so a couple of
    huge spikes (simulating real peaks sitting on the baseline) must not blow up the estimate
    the way a plain std() over the whole trace would."""
    rng = np.random.default_rng(3)
    baseline_std = 5.0
    trace = rng.normal(0.0, baseline_std, size=500)
    trace[100] = 1e5  # spike 1
    trace[300] = 8e4  # spike 2

    noise = global_noise(trace)

    assert noise == pytest.approx(baseline_std, rel=0.5)
    assert noise < np.std(trace) / 10  # unlike naive std, this must not be dominated by spikes


# ---------------------------------------------------------------------------
# Degenerate-input branches
# ---------------------------------------------------------------------------


def test_find_apex_falls_back_to_nearest_sample_when_window_is_empty():
    """rt_expected far outside the trace's time range means the ± window mask matches no
    sample at all. find_apex must fall back to the nearest sample to rt_expected within
    [lo, hi] rather than raising or returning an out-of-range/undefined index."""
    time = _time_axis(0.0, 1.0)
    signal = np.zeros(time.size)

    i = find_apex(time, signal, rt_expected=5.0, window=0.01, lo=0, hi=time.size - 1)

    assert i == time.size - 1  # t=1.0 (the trace's last sample) is nearest to rt_expected=5.0


def test_valley_index_degenerate_adjacent_apexes():
    """When the two apex indices are equal or adjacent there is no sample strictly between
    them to search. valley_index must return i_left rather than indexing an empty slice."""
    signal = np.array([1.0, 2.0, 3.0, 4.0, 5.0])

    assert valley_index(signal, i_left=2, i_right=2) == 2  # i_right ≤ i_left
    assert valley_index(signal, i_left=2, i_right=3) == 2  # i_right == i_left + 1


def test_integrate_region_zero_or_negative_width_returns_zero_area():
    """i_end ≤ i_start (a zero- or negative-width region) must not crash. It returns zero
    area with both baseline endpoints pinned to the single intensity sample at i_start."""
    time = _time_axis(0.0, 1.0)
    signal = np.full(time.size, 7.0)

    area, base_l, base_r = integrate_region(time, signal, i_start=10, i_end=10)
    assert area == 0.0
    assert base_l == base_r == 7.0

    area2, base_l2, base_r2 = integrate_region(time, signal, i_start=10, i_end=5)
    assert area2 == 0.0
    assert base_l2 == base_r2 == 7.0
