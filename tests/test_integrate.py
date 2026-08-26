"""Tests for the RT-anchored integration engine (``core.calibrate`` + ``core.integrate``).

These run against a real acquisition rather than synthetic traces. The engine's whole value
is that its decision tree (resolved / shoulder / absent / edge / co-elution) reproduces an
analyst's manual integration on actual chromatograms, a property invisible on a hand-built
Gaussian. ``tests/test_peaks.py`` covers the numerical primitives analytically. This module
covers the orchestration on real data.
"""

from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest

from quant209.core.calibrate import build_anchors
from quant209.core.config import load_canonical_peaks
from quant209.core.integrate import (
    KNOWN_FLAGS,
    KNOWN_METHODS,
    _cluster_spans,
    _resolve_boundaries,
    integrate_run,
)
from quant209.core.mzml import Chromatogram, get_transition, read_srm_mzml
from quant209.core.peaks import IntegrationParams, PeakResult

from dataset import NAMES


@lru_cache(maxsize=4)
def _chroms(path_str: str):
    """Read one mzML once per session (parsing is the slow part of these tests)."""
    return read_srm_mzml(Path(path_str))


def _find(run_dir: Path, token: str) -> Path:
    matches = [p for p in sorted((run_dir / "Converted").glob("*.mzML")) if token in p.name]
    if not matches:
        pytest.skip(f"no {token} mzML in {run_dir}")
    return matches[0]


@pytest.fixture
def analytes():
    return load_canonical_peaks()


@pytest.fixture
def params():
    return IntegrationParams()


@pytest.fixture
def cal_chroms(reference_run: Path):
    return _chroms(str(_find(reference_run, "Cal01")))


@pytest.fixture
def sample_chroms(reference_run: Path):
    return _chroms(str(_find(reference_run, NAMES.sample)))


@pytest.fixture
def anchors(cal_chroms, analytes, params):
    return build_anchors(cal_chroms, analytes, params)


def _assert_auditable(result: PeakResult) -> None:
    """Every persisted peak must be re-plottable from its own record alone.

    The strict ``rt_start < rt_apex < rt_end`` / ``i_start < i_apex < i_end`` ordering is
    asserted for every *resolved* peak, meaning everything the engine reports as a real peak.
    Peaks flagged ``absent`` are the documented exception: that branch bounds the peak by its
    RT territory rather than by its apex, so the apex may sit outside the integrated region
    (see ``test_absent_peaks_are_bounded_by_their_rt_territory``).
    """
    assert result.rt_start <= result.rt_end, f"{result.label}: RT boundary ordering"
    assert 0 <= result.i_start <= result.i_end, f"{result.label}: index boundary ordering"
    assert np.isfinite(result.baseline_start), f"{result.label}: baseline_start"
    assert np.isfinite(result.baseline_end), f"{result.label}: baseline_end"
    assert np.isfinite(result.area) and result.area >= 0.0, f"{result.label}: area"
    assert np.isfinite(result.sn), f"{result.label}: S/N"
    assert np.isfinite(result.height), f"{result.label}: height"
    assert result.method in KNOWN_METHODS, f"{result.label}: unknown method {result.method!r}"
    assert set(result.flags) <= KNOWN_FLAGS, f"{result.label}: unknown flags {result.flags!r}"
    assert result.detected == ("absent" not in result.flags), f"{result.label}: detected"
    if result.detected:
        assert result.rt_start < result.rt_apex < result.rt_end, f"{result.label}: RT ordering"
        assert 0 <= result.i_start < result.i_apex < result.i_end, f"{result.label}: index order"
        assert result.area > 0.0, f"{result.label}: resolved peak with no area"


# ---------------------------------------------------------------------------
# Tolerant transition resolution (the engine's contract with the reader)
# ---------------------------------------------------------------------------


def _chrom(q1: float, q3: float) -> Chromatogram:
    return Chromatogram(q1=q1, q3=q3, times=np.array([1.0, 2.0]), intensities=np.array([0.0, 1.0]))


def test_get_transition_prefers_an_exact_key():
    """An exact rounded key must never be displaced by a near neighbour."""
    chroms = {(200.0, 165.1): _chrom(200.0, 165.1), (200.1, 165.1): _chrom(200.1, 165.1)}
    assert get_transition(chroms, 200.0, 165.1) is chroms[(200.0, 165.1)]


def test_get_transition_falls_back_to_the_nearest_key():
    """The labelled-standard case: a derived transition ~0.1 Da off the acquired key."""
    chroms = {(200.1, 165.1): _chrom(200.1, 165.1)}
    assert get_transition(chroms, 200.0, 165.1) is chroms[(200.1, 165.1)]


def test_get_transition_refuses_a_distant_key():
    """Beyond the tolerance there is no match. A different transition is never substituted."""
    chroms = {(234.0, 164.1): _chrom(234.0, 164.1)}
    assert get_transition(chroms, 200.0, 165.1) is None


# ---------------------------------------------------------------------------
# RT anchors from the run's own calibration injection
# ---------------------------------------------------------------------------


def test_build_anchors_covers_the_analyte_table(anchors, analytes):
    """The cal injection contains every congener, so nearly every analyte gets an anchor."""
    assert len(anchors) >= 180
    assert all(np.isfinite(rt) for rt in anchors.values())


def test_anchors_stay_near_their_canonical_cal_rt(anchors, analytes, params):
    """An anchor is a *refinement* of the canonical RT, never a jump to another peak."""
    for analyte in analytes:
        if analyte.label in anchors:
            drift = abs(anchors[analyte.label] - analyte.cal_rt)
            assert drift <= params.max_reach + 1e-9, f"{analyte.label} drifted {drift:.3f} min"


# ---------------------------------------------------------------------------
# Integrating the calibration injection itself
# ---------------------------------------------------------------------------


def test_integrate_cal_file_finds_nearly_every_analyte(cal_chroms, analytes, anchors, params):
    """The cal standard holds all 209 congeners: almost every row must integrate positive."""
    results = integrate_run(cal_chroms, analytes, anchors, params)
    positive = [r for r in results if r.area > 0]
    assert len(positive) >= 170, f"only {len(positive)} of {len(analytes)} analytes integrated"


def test_cal_peaks_are_fully_auditable(cal_chroms, analytes, anchors, params):
    """Per-peak persistence is a hard product requirement, not a nicety."""
    results = integrate_run(cal_chroms, analytes, anchors, params)
    assert results
    for result in results:
        _assert_auditable(result)


def test_coelution_group_is_one_summed_peak(cal_chroms, analytes, anchors, params):
    """``129+138+163`` is chromatographically inseparable: one label, one integrated peak."""
    results = integrate_run(cal_chroms, analytes, anchors, params)
    matching = [r for r in results if r.label == "129+138+163"]
    assert len(matching) == 1
    assert matching[0].area > 0


def test_labeled_standards_integrate_via_tolerant_lookup(cal_chroms, analytes, anchors, params):
    """The table's labelled-standard transitions sit ~0.1 Da off the acquired keys.

    Without a tolerant nearest-transition lookup the surrogate/IS channels resolve to no
    chromatogram at all and the whole quantification chain loses its denominator.
    """
    labeled = [a for a in analytes if a.role != "native"]
    assert len(labeled) == 12
    by_label = {r.label: r for r in integrate_run(cal_chroms, analytes, anchors, params)}
    found = [a.label for a in labeled if by_label.get(a.label) and by_label[a.label].area > 0]
    missing = sorted({a.label for a in labeled} - set(found))
    # Hard requirement, not a tolerance: the calibration standard contains all 12 labelled
    # standards by construction, and quantification divides by the two IS responses.
    assert len(found) == 12, f"labeled standards missing from the cal injection: {missing}"


def test_peak_results_carry_the_acquired_transition(cal_chroms, analytes, anchors, params):
    """``src_q1``/``src_q3`` must join back to the reader dict even for labelled channels."""
    results = integrate_run(cal_chroms, analytes, anchors, params)
    for result in results:
        key = (result.src_q1, result.src_q3)
        assert key in cal_chroms, f"{result.label}: source transition {key} not a reader key"
    # The labelled channels are exactly the ones whose canonical key would NOT have joined.
    rekeyed = [r for r in results if (r.q1, r.q3) != (r.src_q1, r.src_q3)]
    assert rekeyed, "expected the tolerant lookup to re-key the labeled-standard channels"
    labeled_labels = {a.label for a in analytes if a.role != "native"}
    assert {r.label for r in rekeyed} <= labeled_labels


def test_resolved_neighbours_split_at_their_shared_valley():
    """Two resolved peaks must be split by a perpendicular drop, reported as ``valley_drop``.

    Pinned on a synthetic pair because it is the branch that decides how co-eluting congeners
    divide their area: without this, inverting the neighbour predicate leaves every other test
    green while silently handing each peak its neighbour's tail.
    """
    n = 300
    idx = np.arange(n, dtype=float)
    detect = (
        1000.0 * np.exp(-0.5 * ((idx - 100.0) / 8.0) ** 2)
        + 1000.0 * np.exp(-0.5 * ((idx - 200.0) / 8.0) ** 2)
    )
    bounds = [(0, 149), (150, n - 1)]
    apexes = [100, 200]

    left = _resolve_boundaries(detect, bounds, apexes, [True, True], 0, halfwidth=50)
    right = _resolve_boundaries(detect, bounds, apexes, [True, True], 1, halfwidth=50)
    assert left[2] == "valley_drop" and right[2] == "valley_drop"
    # The split is the shared valley: neither peak reaches past it into the other.
    assert 100 < left[1] < 200 and 100 < right[0] < 200
    assert left[1] <= right[0]

    # With the neighbour absent, the same peak is bounded by its own baseline edge instead.
    isolated = _resolve_boundaries(detect, bounds, apexes, [True, False], 0, halfwidth=50)
    assert isolated[2] == "rt_window"


def test_valley_drop_is_produced_on_real_data(cal_chroms, analytes, anchors, params):
    """The cal injection is dense with resolved neighbours, so the split branch must fire."""
    results = integrate_run(cal_chroms, analytes, anchors, params)
    split = [r for r in results if r.method == "valley_drop"]
    assert len(split) >= 50, f"only {len(split)} peaks split at a shared valley"
    assert all(r.detected for r in split), "an absent analyte must never claim a valley split"


def test_absent_peaks_are_bounded_by_their_rt_territory(sample_chroms, analytes, anchors, params):
    """Characterizes the one place the apex may fall outside the integrated region.

    When the presence test rejects an analyte, the engine does not integrate around the
    (meaningless) apex. It reports the valley-to-valley span of the analyte's midpoint
    territory, so the number is a bounded estimate of whatever sits at that retention time
    rather than a peak measurement. The apex index is still persisted, and may lie outside
    ``[i_start, i_end]``. Downstream consumers must therefore key off the ``absent`` flag,
    not off the apex, when deciding whether a peak is real.
    """
    results = integrate_run(sample_chroms, analytes, anchors, params)
    outside = [r for r in results if not (r.i_start < r.i_apex < r.i_end)]
    assert outside, "expected some rejected analytes in a real sample"
    for result in outside:
        assert "absent" in result.flags, f"{result.label}: resolved peak with apex out of region"
        assert result.method == "rt_window"


def test_sn_is_not_a_detection_test(sample_chroms, analytes, anchors, params):
    """Characterizes the S/N inversion that makes ``detected`` the only safe presence gate.

    An absent analyte's height is measured from an apex that may sit outside the integrated
    region, against a baseline extrapolated to reach it, so its S/N can come out arbitrarily
    large. Filtering on sn ≥ min_sn would therefore admit most rejected analytes as if they
    were real peaks. This test exists so that inversion stays visible rather than being
    rediscovered downstream.
    """
    results = integrate_run(sample_chroms, analytes, anchors, params)
    absent = [r for r in results if not r.detected]
    assert absent, "expected rejected analytes in a real sample"
    sn_would_pass = [r for r in absent if r.sn >= params.min_sn]
    assert sn_would_pass, "the S/N inversion is the documented contract; it must stay pinned"
    # Their high S/N is exactly why they are not flagged low_sn either.
    assert any("low_sn" not in r.flags for r in sn_would_pass)


def test_labels_are_unique_across_the_run(cal_chroms, analytes, anchors, params):
    """One PeakResult per canonical row, no analyte integrated twice."""
    results = integrate_run(cal_chroms, analytes, anchors, params)
    labels = [r.label for r in results]
    assert len(labels) == len(set(labels))


# ---------------------------------------------------------------------------
# Integrating a sample with the cal file's anchors
# ---------------------------------------------------------------------------


def test_integrate_sample_with_cal_anchors(sample_chroms, analytes, anchors, params):
    """A real sample genuinely lacks some congeners, but most still integrate."""
    results = integrate_run(sample_chroms, analytes, anchors, params)
    positive = [r for r in results if r.area > 0]
    assert len(positive) >= 170, f"only {len(positive)} analytes integrated in the sample"
    for result in results:
        _assert_auditable(result)


def test_sample_flags_low_confidence_peaks(sample_chroms, analytes, anchors, params):
    """Congeners absent from a real sample must be flagged, not silently reported as peaks."""
    results = integrate_run(sample_chroms, analytes, anchors, params)
    flagged = [r for r in results if r.flags]
    assert flagged, "no low-confidence peaks flagged in a real sample"
    assert any("absent" in r.flags for r in flagged)
    assert any("low_sn" in r.flags for r in flagged)
    # An absent analyte is bounded by its RT territory, never split off a neighbour's valley.
    for result in results:
        if "absent" in result.flags:
            assert result.method == "rt_window"


def test_sample_peaks_track_the_cal_anchors(sample_chroms, analytes, anchors, params):
    """Sample apexes stay inside the drift territory around their anchor RT."""
    results = integrate_run(sample_chroms, analytes, anchors, params)
    for result in results:
        drift = abs(result.rt_apex - anchors[result.label])
        assert drift <= params.max_reach + 1e-9, f"{result.label} drifted {drift:.3f} min"


# ---------------------------------------------------------------------------
# Cluster baselines on crowded windows
# ---------------------------------------------------------------------------


def _two_peak_bed() -> tuple[np.ndarray, list[tuple[int, int, str]], list[int]]:
    """Two Gaussians overlapping enough that the valley between them never reaches base.

    Note the elevation has to come from the peaks' own overlap, not from a constant offset:
    the merge test is the valley's height *above the pair's outer boundaries* as a fraction
    of the apex height above them, so a raised floor lifts both terms and changes nothing.
    That ratio is exactly what makes ``overlap_frac`` a tolerance on "did it reach the floor"
    rather than an absolute intensity.
    """
    n = 300
    idx = np.arange(n, dtype=float)
    detect = (
        1000.0 * np.exp(-0.5 * ((idx - 100.0) / 18.0) ** 2)
        + 1000.0 * np.exp(-0.5 * ((idx - 160.0) / 18.0) ** 2)
    )
    resolved = [(40, 130, "valley_drop"), (130, 220, "valley_drop")]
    return detect, resolved, [100, 160]


def test_cluster_spans_merges_two_peaks_that_never_return_to_base():
    """An elevated valley means one bed under both peaks, so they share one baseline span."""
    detect, resolved, apexes = _two_peak_bed()
    spans = _cluster_spans(detect, resolved, apexes, [True, True], overlap_frac=0.05)
    assert spans[0] == spans[1], "an unresolved pair must share one baseline span"
    assert spans[0] == (40, 220), "the shared span is the cluster's outer boundaries"


def test_cluster_spans_leaves_peaks_that_reach_the_floor_alone():
    """Two peaks whose valley reaches the chromatographic floor need no shared baseline."""
    n = 300
    idx = np.arange(n, dtype=float)
    detect = (
        1000.0 * np.exp(-0.5 * ((idx - 60.0) / 8.0) ** 2)
        + 1000.0 * np.exp(-0.5 * ((idx - 240.0) / 8.0) ** 2)
    )
    resolved = [(20, 150, "valley_drop"), (150, 280, "valley_drop")]
    spans = _cluster_spans(detect, resolved, apexes=[60, 240], present=[True, True], overlap_frac=0.05)
    assert spans[0] != spans[1], "separately resolved peaks must keep their own baselines"


def test_an_absent_analyte_breaks_the_cluster_run():
    """Present-only membership restricts which peaks may share a baseline.

    Without the restriction, a floor-anchored baseline drawn under a *shoulder* would hand it
    the dominant neighbour's underlying area. Clustering therefore never crosses an analyte
    the presence test rejected, so the shoulder/absent branch behaves exactly as the rest of
    the engine expects.
    """
    detect, resolved, apexes = _two_peak_bed()
    spans = _cluster_spans(detect, resolved, apexes, [True, False], overlap_frac=0.05)
    assert spans[0] != spans[1], "an absent analyte must never join a cluster"
    assert spans[0] == (40, 130) and spans[1] == (130, 220), "each keeps its own span"


def test_overlap_frac_above_one_reproduces_the_chord_baseline_engine(
    cal_chroms, analytes, anchors
):
    """``overlap_frac`` above 1 disables clustering, so each peak keeps its own chord baseline.

    The ratio is ``(valley - base) / (lower_apex - base)`` and a valley cannot stand higher
    than the lower of the two apexes it lies between, so it cannot exceed 1. The one way it
    could is degenerate rather than chromatographic: ``valley_index`` returns ``i_left`` when
    two apexes are index-adjacent, making the "valley" an apex and the ratio exactly 1, still
    not *greater* than 1, so still no merge, and no such pair occurs on either reference run.
    """
    off = integrate_run(cal_chroms, analytes, anchors, IntegrationParams(overlap_frac=10.0))
    assert not any("clustered" in r.flags for r in off), "clustering fired with the switch off"


def test_clustering_only_ever_raises_an_area(cal_chroms, analytes, anchors):
    """A cluster baseline sits at or below each member's own chord, so no area may fall.

    This is empirical on this data rather than structural: a chord across a long sloping bed
    could in principle under-cut one end. That is why it is pinned here, and why the observed
    maximum cluster length is recorded in ``_cluster_spans``'s docstring.
    """
    shipped = integrate_run(cal_chroms, analytes, anchors, IntegrationParams())
    off = integrate_run(cal_chroms, analytes, anchors, IntegrationParams(overlap_frac=10.0))
    by_label = {r.label: r for r in off}
    clustered = [r for r in shipped if "clustered" in r.flags]
    assert clustered, "the cal injection is crowded; clustering must fire on it"
    for result in shipped:
        chord = by_label[result.label]
        assert result.area >= chord.area - 1e-6, f"{result.label}: clustering lowered the area"
        if "clustered" not in result.flags:
            assert result.area == pytest.approx(chord.area), f"{result.label}: unclustered drift"


def test_clustered_peaks_keep_their_own_boundaries_and_method(cal_chroms, analytes, anchors):
    """Clustering changes only the baseline, never where the boundaries were drawn.

    ``method`` describes boundary selection, so it must stay identical with clustering on or
    off. Only the ``clustered`` flag and the persisted baseline differ. That is what lets an
    auditor recompute a clustered area by hand from the record alone.
    """
    shipped = integrate_run(cal_chroms, analytes, anchors, IntegrationParams())
    off = {r.label: r for r in integrate_run(cal_chroms, analytes, anchors, IntegrationParams(overlap_frac=10.0))}
    for result in shipped:
        chord = off[result.label]
        assert result.method == chord.method, f"{result.label}: method changed with clustering"
        assert (result.i_start, result.i_end) == (chord.i_start, chord.i_end), result.label
        if "clustered" in result.flags:
            assert result.detected, "only present analytes may cluster"


# ---------------------------------------------------------------------------
# Shape-fitting for riders no boundary can separate
# ---------------------------------------------------------------------------


def test_deconvolution_separates_a_rider_from_a_synthetic_parent():
    """Two Gaussians of known area, the smaller buried on the larger's leading edge.

    No boundary rule can recover this: the trace rises monotonically through the rider, so
    wherever a perpendicular falls it takes some of the parent or leaves some of the rider.
    """
    from quant209.core.deconvolve import deconvolve_rider

    times = np.arange(30.0, 32.0, 0.005)
    sigma = 0.03
    rider_rt, parent_rt = 30.90, 31.00
    rider = 2.0e4 * np.exp(-0.5 * ((times - rider_rt) / sigma) ** 2)
    parent = 4.0e5 * np.exp(-0.5 * ((times - parent_rt) / sigma) ** 2)
    signal = rider + parent
    true_area = 2.0e4 * sigma * np.sqrt(2 * np.pi) * 60.0

    parent_apex = int(np.argmin(np.abs(times - parent_rt)))
    got = deconvolve_rider(times, signal, signal, rider_rt, [parent_rt], parent_apex)
    assert got is not None, "a clean two-peak envelope must be resolvable"
    area, residual = got
    assert residual < 0.05, f"fit did not reproduce the trace (residual {residual:.3f})"
    assert abs(area / true_area - 1.0) < 0.15, f"recovered {area:,.0f} of {true_area:,.0f}"


def test_deconvolution_refuses_a_fit_it_cannot_justify():
    """Pure noise has no two-peak decomposition. The caller must be told so, not given a number."""
    from quant209.core.deconvolve import deconvolve_rider

    rng = np.random.default_rng(0)
    times = np.arange(30.0, 32.0, 0.005)
    signal = rng.normal(500.0, 200.0, times.size)
    assert deconvolve_rider(times, signal, signal, 30.9, [31.0], 200) is None


def test_deconvolution_refuses_a_rider_larger_than_its_parent():
    """A 'rider' that is most of the envelope is a mis-identification, not a shoulder."""
    from quant209.core.deconvolve import deconvolve_rider

    times = np.arange(30.0, 32.0, 0.005)
    sigma = 0.03
    big = 4.0e5 * np.exp(-0.5 * ((times - 30.90) / sigma) ** 2)
    small = 2.0e4 * np.exp(-0.5 * ((times - 31.00) / sigma) ** 2)
    signal = big + small
    parent_apex = int(np.argmin(np.abs(times - 31.00)))
    # 30.90 carries most of the envelope, so claiming it as the rider must be refused
    assert deconvolve_rider(times, signal, signal, 30.90, [31.00], parent_apex) is None
