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
from quant209.core.config import Analyte, load_canonical_peaks
from quant209.core.integrate import (
    KNOWN_FLAGS,
    KNOWN_METHODS,
    _cluster_spans,
    _resolve_boundaries,
    integrate_run,
    integrate_transition,
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
    assert got.residual < 0.05, f"fit did not reproduce the trace (residual {got.residual:.3f})"
    assert abs(got.area / true_area - 1.0) < 0.15, f"recovered {got.area:,.0f} of {true_area:,.0f}"
    # The fit also says how tall the rider is and where, which is all a caller with no apex
    # of its own to measure has to go on.
    assert abs(got.height / 2.0e4 - 1.0) < 0.15, f"fitted height {got.height:,.0f} of 20,000"
    assert abs(got.centre - rider_rt) < 0.01, f"fitted centre {got.centre:.3f} for {rider_rt}"


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


# ---------------------------------------------------------------------------
# Named rules: analytes integrated by rule where neither detection test finds them
# ---------------------------------------------------------------------------

#: The order ``KNOWN_FLAGS`` documents, and the order every peak must spell its flags in.
FLAG_ORDER = (
    "absent", "low_sn", "narrow", "clustered", "shoulder", "forced", "displaced", "apportioned",
)

NO_RULES = IntegrationParams(forced_riders=(), forced_present=())


def _analyte(label: str, rt: float) -> Analyte:
    return Analyte(
        label=label, q1=325.9, q3=255.9, homolog="penta", order=1, role="native",
        is_norm="204", cal_mass=1.0, cal_rt=rt,
    )


def _trace(*peaks: tuple[float, float, float]) -> tuple[np.ndarray, np.ndarray]:
    """Two minutes of trace: ``(height, centre, sigma)`` Gaussians on a noisy 400-count floor."""
    times = np.arange(30.0, 32.0, 0.005)
    signal = 400.0 + np.random.default_rng(7).normal(0.0, 60.0, times.size)
    for height, centre, sigma in peaks:
        signal = signal + height * np.exp(-0.5 * ((times - centre) / sigma) ** 2)
    return times, signal


def _by_label(times, signal, analytes, params) -> dict[str, PeakResult]:
    anchors = {a.label: a.cal_rt for a in analytes}
    return {r.label: r for r in integrate_transition(times, signal, analytes, anchors, params)}


def test_a_forced_rider_is_fitted_where_neither_test_finds_it():
    """A rider 3% the height of its parent and 2.5 widths ahead of it: no maximum, no bulge.

    The trace climbs straight through it, so the prominence test sees nothing and the
    curvature test sees nothing either. Without the rule the analyte is absent and the area
    reported for it is a sliver of its territory. With it, the rider is fitted against its
    parent, and the record describes the rider rather than the flank it sits on.
    """
    sigma, parent_rt = 0.03, 31.00
    rider_rt = parent_rt - 2.5 * sigma
    times, signal = _trace((1.2e4, rider_rt, sigma), (4.0e5, parent_rt, sigma))
    analytes = [_analyte("rider", rider_rt), _analyte("parent", parent_rt)]
    rule = IntegrationParams(forced_riders=("rider",), forced_present=())
    true_area = 1.2e4 * sigma * np.sqrt(2 * np.pi) * 60.0

    off = _by_label(times, signal, analytes, NO_RULES)
    on = _by_label(times, signal, analytes, rule)

    assert off["rider"].flags == ["absent"], "the case the rule exists for: found by neither test"
    assert off["rider"].area < 0.01 * true_area

    rider = on["rider"]
    assert rider.flags == ["forced"] and rider.method == "deconvolved"
    assert abs(rider.area / true_area - 1.0) < 0.10, f"fitted {rider.area:,.0f} of {true_area:,.0f}"
    _assert_auditable(rider)
    # Where the rider is and how tall, from the fit. The search's own apex is the territory
    # edge, half way to the parent, and the trace there is mostly the parent.
    assert abs(rider.rt_apex - rider_rt) < 0.01
    assert abs(rider.height / 1.2e4 - 1.0) < 0.10
    assert off["rider"].height > 5 * rider.height, "the flank's height is not the rider's"
    assert "low_sn" not in rider.flags, "a rider 200 times the noise is not near the limit"
    # A rule may add a peak. It never re-cuts one.
    assert on["parent"] == off["parent"]


def test_a_forced_rider_with_nothing_to_ride_on_stays_absent():
    """No resolved peak in reach means no envelope to apportion, so nothing is fitted."""
    times, signal = _trace()
    rule = IntegrationParams(forced_riders=("rider",), forced_present=())
    only = _by_label(times, signal, [_analyte("rider", 30.9)], rule)["rider"]
    assert "absent" in only.flags and "forced" not in only.flags
    assert only.method == "rt_window"


def test_a_named_rule_leaves_a_peak_the_trace_already_shows():
    """The rules apply to what is still unfound. A peak with its own maximum is evidence."""
    times, signal = _trace((5.0e4, 30.9, 0.03))
    analytes = [_analyte("x", 30.9)]
    both = IntegrationParams(forced_riders=("x",), forced_present=("x",))
    named = _by_label(times, signal, analytes, both)
    assert named == _by_label(times, signal, analytes, NO_RULES)
    assert named["x"].flags == []


def test_a_displaced_peak_is_recovered_by_the_forced_present_rule():
    """A concentrated sample: the peak is broad and its top lies past its search territory.

    The search stops on the rising flank, a point with no prominence, and the presence test
    rejects the largest peak on the trace. Bounded as if present, it is integrated from the
    valley it shares with its resolved neighbour to its own far foot.
    """
    sigma, top = 0.06, 31.22
    times, signal = _trace((2.0e4, 30.50, 0.03), (6.0e6, top, sigma))
    analytes = [_analyte("small", 30.50), _analyte("big", 31.00)]
    rule = IntegrationParams(forced_riders=(), forced_present=("big",))
    true_area = 6.0e6 * sigma * np.sqrt(2 * np.pi) * 60.0

    off = _by_label(times, signal, analytes, NO_RULES)
    on = _by_label(times, signal, analytes, rule)

    assert off["big"].flags == ["absent"]
    assert off["big"].area < 0.01 * true_area, "the failure: 54 million counts reported as nothing"

    big = on["big"]
    assert big.flags == ["forced"] and big.method == "valley_drop"
    assert abs(big.area / true_area - 1.0) < 0.02, f"integrated {big.area:,.0f} of {true_area:,.0f}"
    _assert_auditable(big)
    assert abs(big.rt_apex - top) < 0.01, "the apex recorded is the top, not the flank"
    assert on["small"] == off["small"]

    # ``method`` says how the boundaries were chosen. With no resolved neighbour there is no
    # shared valley, and the record must not claim one.
    alone = _by_label(times, signal, analytes[1:], rule)["big"]
    assert alone.flags == ["forced"] and alone.method == "rt_window"


NO_RECLAIM = IntegrationParams(dominance_ratio=0.0)


def _displaced_pair(first: str, second: str, top: float) -> tuple:
    """Two analytes 0.06 min apart, and one large broad peak whose top lies at ``top``.

    The first analyte's territory ends half way to the second, at 31.03 min. A top beyond
    that is where a concentrated sample's dominant peak drifts: the first analyte's search
    stops on the rising flank at its edge, and the second claims the whole peak.
    """
    times, signal = _trace((5.0e6, top, 0.05))
    analytes = [_analyte(first, 31.00), _analyte(second, 31.06)]
    true_area = 5.0e6 * 0.05 * np.sqrt(2 * np.pi) * 60.0
    return times, signal, analytes, true_area


def test_a_dominant_peak_displaced_into_a_minor_neighbour_is_given_back():
    """110 is nine percent of Aroclor 1254 and 115 a tenth of a percent.

    Without the rule, the largest peak on the trace is reported as 115 and 110 as absent,
    which is what the archive shows in its most concentrated samples. With it, 110 takes the
    peak, 115 is absent, and the two windows do not overlap, so nothing is counted twice.
    """
    times, signal, analytes, true_area = _displaced_pair("110", "115", top=31.10)
    off = _by_label(times, signal, analytes, NO_RECLAIM)
    on = _by_label(times, signal, analytes, IntegrationParams())

    assert off["110"].flags == ["absent"] and off["115"].detected
    assert off["115"].area > 0.95 * true_area, "the failure: the dominant peak credited to 115"

    winner, loser = on["110"], on["115"]
    assert winner.flags == ["displaced"] and winner.detected
    assert abs(winner.area / true_area - 1.0) < 0.02, f"{winner.area:,.0f} of {true_area:,.0f}"
    assert abs(winner.rt_apex - 31.10) < 0.01
    _assert_auditable(winner)
    assert not loser.detected and "shoulder" not in loser.flags
    assert loser.area < 0.001 * true_area, "the neighbour must give the area back"
    assert loser.i_start >= winner.i_end, "the two windows must not overlap"


def test_two_comparably_abundant_congeners_are_left_alone():
    """87+125 and 86+97+109+119 are both major. The trace cannot say whose peak it is, and
    neither can the profiles, so the rule does not decide."""
    times, signal, analytes, _ = _displaced_pair("87+125", "86+97+109+119", top=31.10)
    assert _by_label(times, signal, analytes, IntegrationParams()) == _by_label(
        times, signal, analytes, NO_RECLAIM
    )


def test_a_neighbour_of_unknown_abundance_keeps_a_peak_at_its_own_retention_time():
    """The profiles give 102 as zero, yet it is a fifth of a percent of a typical archive
    sample, so the table marks its abundance unknown. A peak sitting where 102 elutes stays
    102's. Only a peak well past it is a drifted 99."""
    times, signal, analytes, _ = _displaced_pair("99", "102", top=31.08)
    assert _by_label(times, signal, analytes, IntegrationParams())["102"].detected

    times, signal, analytes, _ = _displaced_pair("99", "102", top=31.13)
    moved = _by_label(times, signal, analytes, IntegrationParams())
    assert moved["99"].flags == ["displaced"] and not moved["102"].detected


def test_a_congener_of_unknown_abundance_never_takes_a_peak():
    """Unknown is not major. On the archive, letting 102 reclaim from 98 produced seventeen
    reclaims no Aroclor profile supports, because every comparison with nan is false."""
    times, signal, analytes, _ = _displaced_pair("102", "98", top=31.13)
    assert _by_label(times, signal, analytes, IntegrationParams()) == _by_label(
        times, signal, analytes, NO_RECLAIM
    )


def test_a_trace_neighbour_gives_back_a_drifted_peak_even_near_its_own_retention_time():
    """112 is genuinely trace in Aroclors and in the archive. On two of the archive's most
    concentrated files the drifted 99 peak topped out 0.02 min from 112's anchor, and 112 was
    reported at 15,800 ng against the analyst's 28."""
    times, signal, analytes, _ = _displaced_pair("99", "112", top=31.08)
    moved = _by_label(times, signal, analytes, IntegrationParams())
    assert moved["99"].flags == ["displaced"] and not moved["112"].detected


def test_congeners_the_profiles_report_together_never_take_each_others_peaks():
    """The reference reports 91 with 88, so the packaged table gives both the combined
    fraction. On the archive, letting 88 take 91's peaks made every one of them worse."""
    times, signal, analytes, _ = _displaced_pair("88", "91", top=31.13)
    assert _by_label(times, signal, analytes, IntegrationParams()) == _by_label(
        times, signal, analytes, NO_RECLAIM
    )


def test_flags_are_always_spelled_in_one_order(cal_chroms, sample_chroms, analytes, anchors, params):
    """Two peaks carrying the same flags must spell them the same way.

    The flags reach the exports as a joined string, so an order that depends on which
    branch built the peak splits one population into two wherever they are counted.
    """
    assert set(FLAG_ORDER) == KNOWN_FLAGS
    for chroms in (cal_chroms, sample_chroms):
        for result in integrate_run(chroms, analytes, anchors, params):
            assert result.flags == sorted(result.flags, key=FLAG_ORDER.index), result.label


def test_the_named_rules_only_ever_add_their_own_peaks(sample_chroms, analytes, anchors, params):
    """On a real sample: a forced peak is one of the named labels, and nothing else moves."""
    named = set(params.forced_riders) | set(params.forced_present)
    assert named, "the shipped defaults name the congeners the rules exist for"
    shipped = {r.label: r for r in integrate_run(sample_chroms, analytes, anchors, params)}
    plain = {r.label: r for r in integrate_run(sample_chroms, analytes, anchors, NO_RULES)}
    assert not any("forced" in r.flags for r in plain.values())
    for label, result in shipped.items():
        if "forced" not in result.flags:
            assert result == plain[label], f"{label}: changed by a rule that does not name it"
            continue
        assert label in named, f"{label}: forced without being named"
        assert result.detected and "absent" in plain[label].flags
        assert not {"shoulder", "clustered"} & set(result.flags)
        if label in params.forced_riders:
            assert result.method == "deconvolved"
        else:
            assert result.method in {"valley_drop", "rt_window"}


# ---------------------------------------------------------------------------
# The pipeline keeps the integrator's verdict
# ---------------------------------------------------------------------------


def test_the_pipeline_keeps_the_integrators_presence_verdict(reference_run: Path):
    """The pipeline reports each peak's presence exactly as the integrator decided it.

    A congener the calibration standard shows is not thereby present in every sample. A
    pass that cleared ``absent`` from every sample peak whose congener the calibration
    detected changed no area and no mass. It only took the detection verdict away: the peak
    stopped appearing in the QA/QC table, its meaningless S/N became readable, and the
    validation report lost the cause it gives for a missed peak. This pins that it is gone.
    """
    from quant209.core.pipeline import run_analysis
    from quant209.core.sequence import find_data_files, scan_inputs

    entries = scan_inputs(find_data_files(reference_run))
    cal = next(e for e in entries if e.file_type == "cal")
    ref = next(e for e in entries if e.file_type == "ref")
    sample = next(e for e in entries if NAMES.sample in e.name)
    result = run_analysis([cal, ref, sample], converted_dir=reference_run / "Converted")

    records = {r.entry.name: r for r in result.files.values()}
    in_cal = {p.label for p in records[cal.name].peaks if p.detected}
    rejected = {p.label for p in records[sample.name].peaks if not p.detected}
    assert rejected & in_cal, "a sample peak the integrator rejected must stay absent"
    for peak in records[sample.name].peaks:
        _assert_auditable(peak)


def test_two_comparable_congeners_merged_by_drift_are_split_by_a_two_peak_fit():
    """86+97+109+119 and 87+125 are both about 3% of Aroclor 1254, so neither can be given
    the other's peak. In concentrated archive samples the drifted pair merged into one
    envelope, credited wholly to 87+125. A rigid two-peak fit (one width, the calibration's
    spacing, one common drift) reproduced the analysts' split within about 10%."""
    times = np.arange(30.0, 32.0, 0.005)
    sigma, drift = 0.045, 0.06
    signal = 400.0 + np.random.default_rng(11).normal(0.0, 60.0, times.size)
    signal += 1.2e6 * np.exp(-0.5 * ((times - (31.00 + drift)) / sigma) ** 2)
    signal += 1.9e6 * np.exp(-0.5 * ((times - (31.06 + drift)) / sigma) ** 2)
    analytes = [_analyte("86+97+109+119", 31.00), _analyte("87+125", 31.06)]
    k = sigma * np.sqrt(2 * np.pi) * 60.0

    off = _by_label(times, signal, analytes, NO_RECLAIM)
    on = _by_label(times, signal, analytes, IntegrationParams())
    assert "absent" in off["86+97+109+119"].flags, "the failure: one congener credited with both"

    first, second = on["86+97+109+119"], on["87+125"]
    assert first.flags == ["apportioned"] and second.flags == ["apportioned"]
    assert first.method == second.method == "deconvolved"
    assert abs(first.area / (1.2e6 * k) - 1) < 0.10, f"{first.area:,.0f} of {1.2e6 * k:,.0f}"
    assert abs(second.area / (1.9e6 * k) - 1) < 0.10, f"{second.area:,.0f} of {1.9e6 * k:,.0f}"
    _assert_auditable(first)
