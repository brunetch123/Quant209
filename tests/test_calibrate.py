"""Tests for ``core.calibrate``: per-run RT alignment, anchors, and the run-fit guard.

Two real acquisitions carry these tests, because the property under test is a relationship
between a retention-time table and a chromatogram and neither half survives being synthesised:

the reference run
    The run the canonical table was extracted from. Its analyst "Cal RT" column *is*
    ``config.canonical_peaks.csv`` to three decimals. Alignment must find nothing here, and
    more strongly must **change nothing**: the anchors it produces have to be bit-for-bit the
    anchors the pre-alignment ``build_anchors`` produced, or every validated number in
    ``docs/validation/`` moves for no reason.

the alternate-method run
    Acquired under a different GC oven program: congener 1 elutes at 8.337 min against a
    canonical 19.563, and the mono/di time segment spans 5.0-22.6 min against 18.0-34.0. This
    is the run alignment exists for, and the run the guard has to catch without it.

Only the calibration injections are read (alignment and the guard are calibration-injection
operations), so these tests stay cheap even though they run on real data. Both run folders
are named in ``tests/dataset.json``. See ``tests/dataset.py``.
"""

from dataclasses import replace
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest

from quant209.core.calibrate import (
    IDENTITY_MAP,
    IDENTITY_TOL,
    MAP_SOURCES,
    MIN_FIT_ANALYTES,
    RtMap,
    _make_map,  # private, but the identity verdict below needs its own direct test
    _settle_identity,
    MIN_STANDARD_ANCHORS,
    align_rt_map,
    assign_in_elution_order,
    build_anchors,
    rt_fit_quality,
    seed_from_standards,
    standards_warp,
    warp_rt,
)
from quant209.core.config import analytes_by_transition, load_canonical_peaks
from quant209.core.mzml import Chromatogram, read_srm_mzml
from quant209.core.peaks import IntegrationParams, local_maxima

from dataset import run_dir

#: The analyst's own calibration retention time for congener 1 in the alternate-method run,
#: read from the "Cal RT" column of that run's quantification workbooks (sheet row 16, column
#: D, the same cells ``tools/bench_areas.workbook_cal_rts`` reads). Hardcoded rather than
#: read at test time so this test states the ground truth it is asserting against instead of
#: deriving it from a file it also has to trust.
ALT_CONGENER_1_CAL_RT = 8.337

#: The same congener's canonical retention time, for contrast: the whole 11-minute error the
#: alignment has to remove.
CANONICAL_CONGENER_1_CAL_RT = 19.563

#: The tetra channel, the alternate-method run's hardest, and the one the cross-channel consensus
#: stage exists for: 32 analytes, whose unconstrained line explains only 22 of them and lands
#: nine minutes off. See :func:`test_consensus_pins_the_tetra_channel`.
TETRA_TRANSITION = (291.9, 222.0)

#: Intercept of that channel's line once the di channel's consistency constraint is applied,
#: and the intercept it takes without it. Measured, and pinned because nothing else in the
#: suite fails when the consensus stage is disabled.
ALT_TETRA_B_CONSTRAINED = -20.29
ALT_TETRA_B_UNCONSTRAINED = -10.58


@lru_cache(maxsize=4)
def _chroms(path_str: str):
    return read_srm_mzml(Path(path_str))


def _cal(run_dir: Path):
    matches = [p for p in sorted((run_dir / "Converted").glob("*.mzML")) if "Cal01" in p.name]
    if not matches:
        pytest.skip(f"no Cal01 mzML in {run_dir}")
    return _chroms(str(matches[0]))


@pytest.fixture
def analytes():
    return load_canonical_peaks()


@pytest.fixture
def params():
    return IntegrationParams()


@pytest.fixture
def alt_method_run() -> Path:
    """A run acquired under a *different* GC method (see the module docstring)."""
    return run_dir("alt_run")


@pytest.fixture
def ref_cal(reference_run: Path):
    return _cal(reference_run)


@pytest.fixture
def alt_cal(alt_method_run: Path):
    return _cal(alt_method_run)


def _native_transitions(analytes) -> list[tuple[float, float]]:
    """Transition keys carrying enough native analytes for a map to be fitted from them."""
    return [
        key
        for key, group in analytes_by_transition(analytes).items()
        if sum(1 for a in group if a.role == "native") >= MIN_FIT_ANALYTES
    ]


# ---------------------------------------------------------------------------
# Same-method run: the alignment must find nothing and change nothing
# ---------------------------------------------------------------------------


def test_same_method_run_fits_the_identity_map(ref_cal, analytes, params):
    """On the run the canonical table came from, every native channel's line is the identity.

    Asserted on the *fitted* ``a``/``b``, not on the ``identity`` verdict. The verdict is
    what the module does with the fit, and testing it instead would pass no matter how far off
    the fit landed.
    """
    alignment = align_rt_map(ref_cal, analytes, params)
    for key in _native_transitions(analytes):
        rt_map = alignment.maps[key]
        assert rt_map.source == "fitted", f"{key}: expected a fit of its own"
        assert abs(rt_map.a - 1.0) < 1e-3, f"{key}: slope {rt_map.a}"
        assert abs(rt_map.b) < IDENTITY_TOL, f"{key}: intercept {rt_map.b}"
        assert rt_map.identity, f"{key}: a near-identity fit must be reported as the identity"
    assert alignment.is_identity
    # Nothing was moved, so nothing may be re-timed either: an identity channel contributes no
    # matched or interpolated retention times (this is what keeps the anchors bit-identical).
    assert alignment.matched == {}
    assert alignment.interpolated == {}


def test_identity_alignment_leaves_anchors_bit_identical(ref_cal, analytes, params):
    """The regression that matters: aligning a same-method run must not move a single anchor.

    Bit equality, not tolerance. ``build_anchors`` reports the *sample index* of an apex, so a
    shifted search window either lands on the same sample or does not. Anything between "all
    185 identical" and "185 differ" would mean the alignment is nudging searches on a run that
    needs no alignment.
    """
    alignment = align_rt_map(ref_cal, analytes, params)
    aligned = build_anchors(ref_cal, analytes, params, alignment=alignment)
    unaligned = build_anchors(ref_cal, analytes, params)
    assert aligned == unaligned


def test_same_method_run_passes_the_fit_guard(ref_cal, analytes, params):
    report = rt_fit_quality(ref_cal, analytes, build_anchors(ref_cal, analytes, params), params)
    assert report.fit_ok
    assert report.reason == ""
    # The calibration standard contains every analyte, so nearly all of them are peaks in it.
    assert report.native_found_frac > 0.9
    assert report.labeled_found_frac == 1.0


# ---------------------------------------------------------------------------
# Different-method run: the alignment must find the shift, and the guard the mismatch
# ---------------------------------------------------------------------------


def test_different_method_run_is_not_the_identity(alt_cal, analytes, params):
    alignment = align_rt_map(alt_cal, analytes, params)
    assert not alignment.is_identity
    mono = alignment.maps[(188.0, 153.1)]
    assert not mono.identity
    # An eleven-minute displacement is not a fit that could be mistaken for column drift.
    assert abs(mono.displacement(CANONICAL_CONGENER_1_CAL_RT)) > 10.0


def test_alignment_recovers_congener_1_on_a_different_gc_method(alt_cal, analytes, params):
    """The headline: congener 1 anchors where the analyst put it, from the cal file alone.

    No workbook is read here: only the alternate-method run's Cal01 mzML. The 8.337 min
    it is checked against is the analyst's own number for this run (see
    :data:`ALT_CONGENER_1_CAL_RT`), which the engine has no access to.
    """
    alignment = align_rt_map(alt_cal, analytes, params)
    anchors = build_anchors(alt_cal, analytes, params, alignment=alignment)
    assert anchors["1"] == pytest.approx(ALT_CONGENER_1_CAL_RT, abs=0.15)
    # Without alignment the same call lands eleven minutes away: the failure being fixed.
    unaligned = build_anchors(alt_cal, analytes, params)
    assert abs(unaligned["1"] - ALT_CONGENER_1_CAL_RT) > 5.0


def test_consensus_pins_the_tetra_channel(alt_cal, analytes, params):
    """The cross-channel consensus stage, which nothing else in this suite exercises.

    The alternate-method run's 32-analyte tetra channel is where a per-channel fit alone goes
    wrong: its unconstrained line explains 22 of its 32 analytes and sits nine minutes off,
    which is a *better-scoring* answer than the truth on that channel considered in isolation.
    Only the di channel next to it (10 of 11 explained, and overlapping in retention time)
    rules it out.

    Pinned on the intercept because that is where the failure shows: disabling the constraint
    moves the channel's line from about -20.3 to about -10.6 while leaving every headline
    number in this suite green (stage 4's pooled table absorbs most of the error), so without
    this assertion the whole stage could be deleted unnoticed.
    """
    alignment = align_rt_map(alt_cal, analytes, params)
    tetra = alignment.maps[TETRA_TRANSITION]
    assert tetra.source in MAP_SOURCES
    assert tetra.b == pytest.approx(ALT_TETRA_B_CONSTRAINED, abs=0.15), (
        f"tetra intercept {tetra.b:.4f}: consensus stage not applied "
        f"(unconstrained lands near {ALT_TETRA_B_UNCONSTRAINED})"
    )
    # Not merely "different from the unconstrained value": far from it, by a margin no
    # amount of fit noise could cross.
    assert abs(tetra.b - ALT_TETRA_B_UNCONSTRAINED) > 5.0

    # And it reaches the anchors: the tetra congeners land on the analyst's own peaks. 81 and
    # 77 are the channel's late pair: canonical 44.058 / 44.636, analyst 36.194 / 37.009 for
    # this run (same workbook column as ALT_CONGENER_1_CAL_RT), where a nine-minute error in
    # the line would be impossible to miss.
    anchors = build_anchors(alt_cal, analytes, params, alignment=alignment)
    assert anchors["81"] == pytest.approx(36.194, abs=0.05)
    assert anchors["77"] == pytest.approx(37.009, abs=0.05)
    tetra_analytes = analytes_by_transition(analytes)[TETRA_TRANSITION]
    assert all(17.0 < anchors[a.label] < 38.0 for a in tetra_analytes)


def test_fit_guard_catches_the_method_mismatch(alt_cal, analytes, params):
    """``fit_ok`` is False without alignment and True with it: the point of the guard.

    The mismatch is otherwise silent: integration returns areas either way, and only the
    calibration injection (where every analyte is known to be present) can say that the
    retention-time table does not describe this run.
    """
    unaligned = rt_fit_quality(alt_cal, analytes, build_anchors(alt_cal, analytes, params), params)
    assert not unaligned.fit_ok
    assert "GC method" in unaligned.reason
    # The labelled channels are the sharp indicator: single peaks on quiet traces, so a
    # misplaced window finds nothing at all.
    assert unaligned.labeled_found_frac < 0.5

    alignment = align_rt_map(alt_cal, analytes, params)
    aligned = rt_fit_quality(
        alt_cal, analytes, build_anchors(alt_cal, analytes, params, alignment=alignment), params
    )
    assert aligned.fit_ok
    assert aligned.reason == ""
    assert aligned.native_found_frac > 0.9
    # Labelled channels carry one analyte each, so since MIN_FIT_ANALYTES they inherit the map
    # from the well-populated channels instead of fitting their own. A channel with nothing
    # to contradict a wrong match can be off by minutes with no signal that anything is wrong.
    # Inheriting places 10 of the 12 here rather than 11, which costs nothing
    # downstream (this run's mass agreement and recovery MAE are unchanged) and is far above
    # the 0.5 the guard actually tests.
    assert aligned.labeled_found_frac > 0.8


# ---------------------------------------------------------------------------
# Synthetic control: a known shift, recovered exactly
# ---------------------------------------------------------------------------


def test_alignment_recovers_a_known_stretch(ref_cal, analytes, params):
    """A known *non-unit* slope, which the shift control cannot see.

    ``times * 1.30 - 5.0`` is what a slower oven ramp does: everything later, and later by
    proportionally more. A pure-offset estimator (or a slope that quietly collapsed to 1)
    would pass every other test in this module and fail here, because on the shift control
    ``a = 1`` is both the truth and the degenerate answer.
    """
    slope, intercept = 1.30, -5.0
    stretched = {
        key: Chromatogram(
            q1=chrom.q1,
            q3=chrom.q3,
            times=chrom.times * slope + intercept,
            intensities=chrom.intensities,
        )
        for key, chrom in ref_cal.items()
    }
    alignment = align_rt_map(stretched, analytes, params)
    assert not alignment.is_identity
    for key in _native_transitions(analytes):
        rt_map = alignment.maps[key]
        assert rt_map.a == pytest.approx(slope, abs=1e-3), f"{key}: slope {rt_map.a}"
        assert rt_map.b == pytest.approx(intercept, abs=0.05), f"{key}: intercept {rt_map.b}"

    # The anchors follow the same map: each lands where the known transform puts it.
    baseline = build_anchors(ref_cal, analytes, params)
    moved = build_anchors(stretched, analytes, params, alignment=alignment)
    residual = np.array(
        [moved[label] - (slope * baseline[label] + intercept) for label in baseline]
    )
    assert np.abs(residual).max() < 0.1


def test_alignment_recovers_a_known_two_minute_shift(ref_cal, analytes, params):
    """Move a real calibration injection +2.0 min in time. The fit must report a=1, b=+2.

    A real trace re-timed by a known amount is the one case where the right answer is known
    exactly, so this pins the estimator itself rather than its agreement with an analyst. The
    shift is applied to the time axis (not by rolling the samples), which is what a slower oven
    ramp does to every peak on every channel at once.
    """
    shift = 2.0
    shifted = {
        key: Chromatogram(
            q1=chrom.q1, q3=chrom.q3, times=chrom.times + shift, intensities=chrom.intensities
        )
        for key, chrom in ref_cal.items()
    }
    alignment = align_rt_map(shifted, analytes, params)
    assert not alignment.is_identity
    for key in _native_transitions(analytes):
        rt_map = alignment.maps[key]
        assert rt_map.a == pytest.approx(1.0, abs=1e-3), f"{key}: slope {rt_map.a}"
        assert rt_map.b == pytest.approx(shift, abs=0.05), f"{key}: intercept {rt_map.b}"

    # And it carries through to the anchors: every one moves by the shift and nothing else.
    baseline = build_anchors(ref_cal, analytes, params)
    moved = build_anchors(shifted, analytes, params, alignment=alignment)
    offsets = np.array([moved[label] - baseline[label] for label in baseline])
    assert np.allclose(offsets, shift, atol=0.05)


# ---------------------------------------------------------------------------
# Degenerate inputs
# ---------------------------------------------------------------------------


def test_alignment_of_an_empty_calibration_file_is_the_identity(analytes, params):
    """Nothing to fit means nothing to change: never a partial or invented map."""
    alignment = align_rt_map({}, analytes, params)
    assert alignment.is_identity
    assert alignment.matched == {}
    assert all(m is IDENTITY_MAP for m in alignment.maps.values())
    assert alignment.expected_rt(analytes[0]) == analytes[0].cal_rt


def test_identity_map_is_a_no_op_on_expected_rt():
    assert IDENTITY_MAP(42.5) == 42.5
    assert IDENTITY_MAP.displacement(42.5) == 0.0


def test_run_identity_verdict_covers_channels_that_borrowed_a_line():
    """A borrower cannot go active on a run whose fitted channels all said "no change".

    The failure this closes is latent rather than live: on the shipped canonical table every
    labelled channel elutes inside or barely past its donor's range, so the arithmetic never
    quite gets there. It is a real edge nonetheless: the identity test is a *displacement*
    judged over each channel's own retention span, and a borrowed line evaluated far from
    where it was fitted displaces further. Constructed directly here, because a run
    exhibiting it does not exist in this dataset and a test that cannot fail is not a test.
    """
    fitted = (256.0, 186.0)
    borrower = (497.7, 427.7)
    # A slope one part in a thousand off identity: nothing across the donor's own 26-37 min
    # span, half a search window by the time a borrower at 66 min evaluates it.
    line = RtMap(a=0.999, b=0.0, identity=True, source="fitted", n_analytes=20, n_matched=20)
    borrowed = _make_map(line.a, line.b, np.array([66.067]), "inherited", 1, 0)
    assert not borrowed.identity, "precondition: the borrower would otherwise be applied"

    settled = _settle_identity({fitted: line, borrower: borrowed}, [fitted])
    assert settled[borrower].identity, "a borrower outlived its donor's identity verdict"
    # The fitted numbers survive the verdict. Provenance is not what is being suppressed.
    assert settled[borrower].a == borrowed.a
    assert settled[borrower].source == "inherited"
    assert settled[borrower](66.067) == 66.067

    # And the rule does not fire when a fitted channel genuinely moved.
    moved = replace(line, identity=False)
    assert not _settle_identity({fitted: moved, borrower: borrowed}, [fitted])[borrower].identity


def test_local_maxima_reports_plateau_left_edges():
    """``peaks.local_maxima``'s contract, which the alignment's peak detection depends on.

    A flat top must report *one* index (its left edge) rather than none or every sample in it:
    a digitised MRM apex is frequently two or three equal samples wide, and a detector that
    dropped those would lose exactly the strong, saturated peaks the alignment relies on.
    """
    assert local_maxima(np.array([0.0, 1.0, 0.0])).tolist() == [1]
    assert local_maxima(np.array([0.0, 1.0, 1.0, 1.0, 0.0])).tolist() == [1]
    assert local_maxima(np.array([0.0, 1.0, 2.0])).tolist() == []  # endpoints are not maxima


def test_build_anchors_default_is_the_unaligned_behaviour(ref_cal, analytes, params):
    """``alignment=None`` must stay the pre-alignment path exactly (callers rely on it)."""
    assert build_anchors(ref_cal, analytes, params, alignment=None) == build_anchors(
        ref_cal, analytes, params
    )


# ---------------------------------------------------------------------------
# Elution-order anchor assignment
# ---------------------------------------------------------------------------


def test_assignment_preserves_order_and_never_double_claims():
    """Two analytes can never be matched to the same peak, and order is a constraint."""
    predicted = [10.0, 10.4, 10.8]
    detected = [10.05, 10.45, 10.85]
    assigned = assign_in_elution_order(predicted, detected, tol=0.5)
    assert assigned == [0, 1, 2]
    placed = [j for j in assigned if j is not None]
    assert placed == sorted(placed), "assignment must be monotone"
    assert len(set(placed)) == len(placed), "a peak was claimed twice"


def test_assignment_recovers_a_window_whose_spacing_changed():
    """The failure a nearest-peak rule cannot handle, and the reason this exists.

    The prediction is uniformly early by more than the peak spacing, so matching each
    analyte to its own nearest peak hands every one of them its neighbour's peak. The
    order constraint recovers the true correspondence instead.
    """
    detected = [52.07, 53.49, 54.53, 54.88, 55.03]
    predicted = [51.6, 53.0, 54.0, 54.4, 54.55]  # every one ~0.45 min early
    nearest = [min(range(len(detected)), key=lambda j: abs(p - detected[j])) for p in predicted]
    assert nearest != [0, 1, 2, 3, 4], "expected nearest-peak matching to fail here"
    assert assign_in_elution_order(predicted, detected, tol=1.2) == [0, 1, 2, 3, 4]


def test_assignment_skips_extra_detected_peaks():
    """Contamination and isotope bleed add peaks. Real analytes must not be forced onto them."""
    detected = [10.0, 10.5, 11.0, 12.9]  # the last is not a real analyte
    assigned = assign_in_elution_order([10.02, 10.52, 11.02], detected, tol=0.5)
    assert assigned == [0, 1, 2]


def test_assignment_leaves_an_unmatchable_analyte_unplaced():
    """An analyte with no peak within tolerance is left for the aligned prediction."""
    assigned = assign_in_elution_order([10.0, 20.0], [10.05], tol=0.5)
    assert assigned[0] == 0
    assert assigned[1] is None


def test_anchors_are_unchanged_on_a_same_method_run(ref_cal, analytes, params):
    """The reference run must not move: its predictions are already right, so the
    assignment has nothing to correct and every anchor stays where it was."""
    alignment = align_rt_map(ref_cal, analytes, params)
    anchors = build_anchors(ref_cal, analytes, params, alignment=alignment)
    # Every anchored analyte still sits within one integration window of its canonical RT.
    for analyte in analytes:
        if analyte.label in anchors and np.isfinite(analyte.cal_rt):
            drift = abs(anchors[analyte.label] - analyte.cal_rt)
            assert drift <= params.max_reach + 1e-9, f"{analyte.label} drifted {drift:.3f}"


# ---------------------------------------------------------------------------
# The run's clock, read off the labelled standards
# ---------------------------------------------------------------------------


def test_standards_are_identified_without_any_retention_time_prior(ref_cal, analytes, params):
    """Each labelled standard sits alone on its own transition, so its peak is unambiguous.

    This is the one measurement in a calibration injection that needs no prior: a native
    channel carries twenty or thirty similar peaks in a row and can be matched several ways,
    but a labelled channel carries one towering peak and nothing to confuse it with.
    """
    pairs = standards_warp(ref_cal, analytes, params)
    assert len(pairs) >= MIN_STANDARD_ANCHORS, "too few standards identified to pin a warp"
    canonicals = [c for c, _ in pairs]
    observed = [o for _, o in pairs]
    assert canonicals == sorted(canonicals) and observed == sorted(observed), "warp must increase"
    # The canonical table came from this run, so the standards sit where it says.
    assert all(abs(o - c) <= 0.05 for c, o in pairs), "same-method run should be near identity"


def test_a_same_method_run_is_not_reseeded(ref_cal, analytes, params):
    """A run whose standards already sit where the table says is returned untouched."""
    assert seed_from_standards(ref_cal, analytes, params) is analytes


def test_a_retimed_run_is_moved_onto_its_own_clock(alt_cal, analytes, params):
    """This run was acquired under a different GC method. The seed must follow the standards."""
    seeded = seed_from_standards(alt_cal, analytes, params)
    assert seeded is not analytes
    by_label = {a.label: a for a in analytes}
    moved = [
        abs(s.cal_rt - by_label[s.label].cal_rt)
        for s in seeded
        if np.isfinite(s.cal_rt) and np.isfinite(by_label[s.label].cal_rt)
    ]
    assert max(moved) > 1.0, "a re-timed run must be moved substantially"
    # order is a chemical fact and cannot be rearranged by re-timing
    order_before = [a.label for a in sorted(
        (a for a in analytes if np.isfinite(a.cal_rt)), key=lambda a: a.cal_rt)]
    order_after = [a.label for a in sorted(
        (a for a in seeded if np.isfinite(a.cal_rt)), key=lambda a: a.cal_rt)]
    assert order_before == order_after


def test_warp_interpolates_and_extends_at_its_own_slope():
    """Between anchors the warp interpolates. Outside them it continues, it does not flatten."""
    pairs = [(20.0, 10.0), (40.0, 34.0), (60.0, 58.0)]
    assert warp_rt(pairs, 30.0) == pytest.approx(22.0)
    assert warp_rt(pairs, 10.0) == pytest.approx(10.0 - 1.2 * 10.0)  # first segment slope 1.2
    assert warp_rt(pairs, 70.0) == pytest.approx(58.0 + 1.2 * 10.0)  # last segment slope 1.2
    assert warp_rt([], 33.3) == 33.3


def test_a_line_shifted_by_one_peak_loses_to_the_composition_of_the_standard():
    """A crowded channel can be matched one peak off and score a match more.

    Sixteen congeners half a minute apart, four of them co-elution groups three times as
    heavy, so three times the peak. The first congener's peak is missing (unresolved) and an
    unrelated peak follows the last. The line shifted by one position puts all sixteen on
    peaks; the right line, fifteen. Only the standard's composition says which is right: the
    big groups on the big peaks.
    """
    from quant209.core.calibrate import _fit_line

    canonical = 10.0 + 0.5 * np.arange(16)
    weights = np.ones(16)
    weights[[2, 6, 10, 13]] = 3.0
    observed = 10.5 + 0.5 * np.arange(16)  # congener k at 10 + 0.5k; 10.0 missing, 18.0 extra
    heights = np.array([weights[int(round((t - 10.0) / 0.5))] * 1000.0 if t < 18.0 else 1000.0
                        for t in observed])

    a, b, n = _fit_line(canonical, observed, 0.1)
    assert abs(b - 0.5) < 0.05 and n == 16, "matches alone pick the shifted line"
    a, b, n = _fit_line(canonical, observed, 0.1, weights=weights, heights=heights)
    assert abs(a - 1.0) < 0.01 and abs(b) < 0.05 and n == 15, "the composition picks the right one"


def test_the_composition_does_not_overrule_a_small_channel():
    """Eight congeners are too few for a correlation of their masses to mean anything."""
    from quant209.core.calibrate import _fit_line

    canonical = np.arange(8, dtype=float) + 10.0
    weights = np.array([1.0, 1.0, 3.0, 1.0, 1.0, 1.0, 3.0, 1.0])
    observed = np.arange(11.0, 19.0)
    heights = np.array([weights[int(t - 10)] * 1000.0 if t < 18 else 1000.0 for t in observed])
    assert _fit_line(canonical, observed, 0.1, weights=weights, heights=heights) ==         _fit_line(canonical, observed, 0.1)
