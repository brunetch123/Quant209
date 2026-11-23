"""Formula-layer tests for the quantification chain, plus one end-to-end pipeline run.

Every number in the unit tests below is hand-computable from the formulas decoded from an
analyst workbook, verified cell by cell against one sample's sheet, row 16: congener 1 area
6353, RRF 2.209448, IS factor 2.48916e-4, mass 3.493954 ng. A regression here is therefore a
broken formula rather than a shifted benchmark.

The end-to-end tests need a real analysis run, chosen by ``tests/dataset.py``. They skip when
no archive is configured.
"""

from __future__ import annotations

import math
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import pytest

from quant209.core.config import Analyte
from quant209.core.peaks import IntegrationParams, PeakResult
from quant209.core.qaqc import (
    QaqcParams,
    censor_below_loq,
    compute_loq,
    integration_flags,
    recovery_flags,
    sample_flags,
)
from quant209.core.quantify import (
    CalibrationModel,
    average_calibrations,
    build_calibration,
    compute_recoveries,
    correct_for_recovery,
    internal_standard_areas,
    is_factor,
    quantify_masses,
)
from quant209.core.standards import IS_SPIKE_NG

from dataset import NAMES, run_dir

# --------------------------------------------------------------------------------------
# Fixtures: a miniature two-analyte + one-IS + one-surrogate canonical table
# --------------------------------------------------------------------------------------


def analyte(
    label: str,
    *,
    role: str = "native",
    homolog: str = "mono",
    is_norm: str = "30D",
    cal_mass: float = 29.06,
    q1: float = 188.0,
    q3: float = 153.1,
    cal_rt: float = 20.0,
    order: int = 1,
) -> Analyte:
    return Analyte(
        label=label,
        q1=q1,
        q3=q3,
        homolog=homolog,
        order=order,
        role=role,
        is_norm=is_norm,
        cal_mass=cal_mass,
        cal_rt=cal_rt,
    )


def peak(label: str, area: float, *, absent: bool = False, flags: list[str] | None = None) -> PeakResult:
    return PeakResult(
        label=label,
        q1=188.0,
        q3=153.1,
        rt_apex=20.0,
        rt_start=19.9,
        rt_end=20.1,
        area=area,
        height=area / 10.0,
        baseline_start=0.0,
        baseline_end=0.0,
        sn=50.0,
        method="valley_drop",
        flags=list(flags or []) + (["absent"] if absent else []),
    )


@pytest.fixture
def mini_analytes() -> list[Analyte]:
    """One native, its 13C surrogate, and the low-Cl internal standard."""
    return [
        analyte("1", cal_mass=29.06),
        analyte("13C 3", role="surrogate_13C", homolog="", cal_mass=61.86),
        analyte("30D", role="IS_d30", homolog="", cal_mass=10.4125),
    ]


# --------------------------------------------------------------------------------------
# RRF
# --------------------------------------------------------------------------------------


def test_rrf_is_hand_computable(mini_analytes):
    """RRF = (cal_mass / cal_area) / (cal_mass_IS / cal_area_IS), exactly."""
    cal = build_calibration(
        [peak("1", 1e6), peak("30D", 5e5), peak("13C 3", 2e5)], mini_analytes
    )
    expected = (29.06 / 1e6) / (10.4125 / 5e5)
    assert cal.rrf["1"] == pytest.approx(expected, rel=1e-12)
    assert cal.rrf["1"] == pytest.approx(1.3954381752701, rel=1e-9)
    # The internal standard is its own reference, so its RRF is exactly 1. The workbook
    # prints a literal 1 in that cell, in the cal sheet's own IS row.
    assert cal.rrf["30D"] == pytest.approx(1.0, rel=1e-12)
    assert cal.is_response["30D"] == pytest.approx(10.4125 / 5e5, rel=1e-12)


def test_rrf_matches_the_analyst_workbook_cell(mini_analytes):
    """The analyst's cal sheet: congener 1 area 115270.09, 30D area 306057.09 -> RRF 2.209448."""
    analytes = [analyte("1", cal_mass=29.0606009064623), analyte("30D", role="IS_d30", cal_mass=34.922610531371)]
    cal = build_calibration([peak("1", 115270.09), peak("30D", 306057.09)], analytes)
    assert cal.rrf["1"] == pytest.approx(2.209448, abs=5e-7)


def test_mis_integrated_cal_peak_is_refused_not_quantified(mini_analytes):
    """A native whose cal peak collapsed gets no RRF at all, refused by the floor guard.

    The calibration standard contains all 209 congeners, so a near-zero *native* cal area is
    a failed integration, not a real absence. RRF = (mass/area)/IS would explode, turning any
    sample signal into an absurd mass, so such analytes are refused rather than reported.
    """
    analytes = [analyte(str(i), cal_mass=29.06) for i in range(1, 11)] + [
        analyte("30D", role="IS_d30", cal_mass=10.4125)
    ]
    peaks = [peak(str(i), 1e6) for i in range(1, 10)] + [peak("10", 5.0), peak("30D", 5e5)]
    cal = build_calibration(peaks, analytes)
    assert "10" not in cal.rrf
    assert "10" in cal.invalid
    assert cal.area_floor == pytest.approx(max(500.0, 0.02 * 1e6))
    assert all(str(i) in cal.rrf for i in range(1, 10))


def test_cal_floor_never_drops_below_its_absolute_minimum(mini_analytes):
    """On a quiet cal injection the floor is the 500 counts·s absolute minimum, not 2% of it."""
    analytes = [analyte(str(i)) for i in range(1, 5)] + [analyte("30D", role="IS_d30", cal_mass=10.4125)]
    peaks = [peak(str(i), 1000.0) for i in range(1, 5)] + [peak("30D", 5e5)]
    cal = build_calibration(peaks, analytes)
    assert cal.area_floor == pytest.approx(500.0)
    assert set(cal.rrf) == {"1", "2", "3", "4", "30D"}


def test_undetected_native_cal_peak_is_refused(mini_analytes):
    """Presence, not S/N, decides whether a cal peak may define a response factor.

    A peak the presence test rejected is not a peak: its area is a bounded estimate of
    whatever sits in the analyte's RT territory. Gating here is on ``PeakResult.detected``,
    never on ``sn``, which an absent peak reports large and meaningless (see ``core.peaks``).
    """
    analytes = [analyte(str(i)) for i in range(1, 11)] + [analyte("30D", role="IS_d30", cal_mass=10.4125)]
    peaks = [peak(str(i), 1e6) for i in range(1, 10)] + [
        peak("10", 9e5, absent=True),  # big area, huge S/N, but not a resolved peak
        peak("30D", 5e5),
    ]
    cal = build_calibration(peaks, analytes)
    assert "10" not in cal.rrf
    assert "not detected" in cal.invalid["10"]


def test_labeled_standards_are_exempt_from_the_native_floor():
    """The floor guard is a *native* rule: a labelled channel legitimately sits low."""
    analytes = [analyte(str(i)) for i in range(1, 11)] + [
        analyte("13C 3", role="surrogate_13C", homolog="", cal_mass=61.86),
        analyte("30D", role="IS_d30", cal_mass=10.4125),
    ]
    peaks = [peak(str(i), 1e6) for i in range(1, 11)] + [peak("13C 3", 100.0), peak("30D", 5e5)]
    cal = build_calibration(peaks, analytes)
    assert "13C 3" in cal.rrf


def test_missing_internal_standard_cal_peak_yields_no_rrfs(mini_analytes):
    """With no IS response there is nothing to normalise against, so nothing is quantifiable."""
    cal = build_calibration([peak("1", 1e6)], mini_analytes)
    assert cal.rrf == {}
    assert "1" in cal.invalid


def test_replicate_cal_injections_are_averaged(mini_analytes):
    """RRF from replicate cal injections is the per-analyte mean."""
    one = build_calibration([peak("1", 1e6), peak("30D", 5e5)], mini_analytes, name="Cal01")
    two = build_calibration([peak("1", 2e6), peak("30D", 5e5)], mini_analytes, name="Cal02")
    avg = average_calibrations([one, two])
    assert avg.rrf["1"] == pytest.approx((one.rrf["1"] + two.rrf["1"]) / 2)
    assert avg.injections == ("Cal01", "Cal02")


def test_an_analyte_valid_in_any_injection_survives_averaging(mini_analytes):
    """Invalid only where *no* replicate produced a usable response."""
    analytes = [analyte(str(i)) for i in range(1, 11)] + [analyte("30D", role="IS_d30", cal_mass=10.4125)]
    good = build_calibration(
        [peak(str(i), 1e6) for i in range(1, 11)] + [peak("30D", 5e5)], analytes, name="Cal01"
    )
    bad = build_calibration(
        [peak(str(i), 1e6) for i in range(1, 10)] + [peak("10", 1.0), peak("30D", 5e5)],
        analytes,
        name="Cal02",
    )
    avg = average_calibrations([good, bad])
    assert "10" in avg.rrf  # one good injection is enough
    assert "10" not in avg.invalid
    assert avg.rrf["10"] == pytest.approx(good.rrf["10"])  # only the usable injection contributes


def test_averaging_a_single_calibration_is_the_identity(mini_analytes):
    one = build_calibration([peak("1", 1e6), peak("30D", 5e5)], mini_analytes, name="Cal01")
    assert average_calibrations([one]) is one


def test_averaging_nothing_is_an_error():
    with pytest.raises(ValueError):
        average_calibrations([])


# --------------------------------------------------------------------------------------
# Mass
# --------------------------------------------------------------------------------------


def test_mass_formula_reproduces_the_workbook_row(mini_analytes):
    """Mass = RRF * Area_sample * (IS_spike / IS_area_sample), to the workbook's own digits."""
    analytes = [
        analyte("1", cal_mass=29.0606009064623),
        analyte("30D", role="IS_d30", cal_mass=34.922610531371),
    ]
    cal = build_calibration([peak("1", 115270.09), peak("30D", 306057.09)], analytes)
    peaks = [peak("1", 6353.0), peak("30D", 41831.21)]
    is_areas = internal_standard_areas(peaks, analytes)
    masses = quantify_masses(peaks, analytes, cal, is_areas)
    assert masses["1"] == pytest.approx(3.493954, abs=1e-5)
    # The internal standard quantifies to its own spike amount, by construction.
    assert masses["30D"] == pytest.approx(IS_SPIKE_NG, rel=1e-12)


def test_is_factor_is_the_drift_correction(mini_analytes):
    """The IS factor is this sample's own IS response, which is what corrects for drift."""
    peaks = [peak("30D", 41831.21)]
    is_areas = internal_standard_areas(peaks, mini_analytes)
    assert is_factor(is_areas, "30D", mini_analytes) == pytest.approx(IS_SPIKE_NG / 41831.21)


def test_undetected_internal_standard_blocks_quantification(mini_analytes):
    """An absent IS peak is not an IS area: without it a sample cannot be quantified.

    Gated on ``PeakResult.detected``. The absent peak below reports S/N 50. Gating on ``sn``
    would have accepted it and scaled the whole sample by a noise measurement.
    """
    cal = build_calibration([peak("1", 1e6), peak("30D", 5e5)], mini_analytes)
    peaks = [peak("1", 1000.0), peak("30D", 41831.21, absent=True)]
    is_areas = internal_standard_areas(peaks, mini_analytes)
    assert is_areas == {}
    assert is_factor(is_areas, "30D", mini_analytes) is None
    assert quantify_masses(peaks, mini_analytes, cal, is_areas) == {}


def test_zero_area_internal_standard_blocks_quantification(mini_analytes):
    peaks = [peak("30D", 0.0)]
    assert internal_standard_areas(peaks, mini_analytes) == {}


def test_analytes_without_an_rrf_are_not_quantified(mini_analytes):
    """An analyte refused at the calibration step stays unreported, never a guessed mass."""
    cal = CalibrationModel(
        rrf={"30D": 1.0},
        cal_area={},
        is_response={"30D": 10.4125 / 5e5},
        invalid={"1": "cal peak below floor"},
        injections=("Cal01",),
        area_floor=500.0,
    )
    peaks = [peak("1", 1e5), peak("30D", 41831.21)]
    masses = quantify_masses(peaks, mini_analytes, cal, internal_standard_areas(peaks, mini_analytes))
    assert "1" not in masses


def test_mass_is_computed_from_the_integrated_area_whatever_the_detection_verdict(mini_analytes):
    """Detection is reported, not applied, at the mass step, matching what the analyst does.

    The analyst types an area for a shoulder they chose to split off a dominant neighbour and
    quantifies it, and so does this chain. ``PeakResult.detected`` travels alongside the mass
    (see ``RunResult.detected``) so the report and the GUI can gate on it, but zeroing here
    would silently disagree with the ground truth the whole benchmark is scored against.
    """
    cal = build_calibration([peak("1", 1e6), peak("30D", 5e5)], mini_analytes)
    peaks = [peak("1", 1e5, absent=True), peak("30D", 5e5)]
    masses = quantify_masses(peaks, mini_analytes, cal, internal_standard_areas(peaks, mini_analytes))
    assert masses["1"] > 0


# --------------------------------------------------------------------------------------
# Recovery
# --------------------------------------------------------------------------------------


def test_recovery_is_surrogate_mass_over_reference_surrogate_mass():
    sample = {"13C 3": 6.355354, "13C 15": 6.0848}
    reference = {"13C 3": 8.07382, "13C 15": 7.79204}
    rec = compute_recoveries(sample, reference)
    assert rec["mono"] == pytest.approx(0.78716, abs=1e-5)  # the workbook's own figure for this sample
    assert rec["di"] == pytest.approx(0.78090, abs=1e-5)


def test_recovery_is_absent_when_the_reference_surrogate_is_missing_or_zero():
    assert compute_recoveries({"13C 3": 5.0}, {}) == {}
    assert compute_recoveries({"13C 3": 5.0}, {"13C 3": 0.0}) == {}
    assert compute_recoveries({}, {"13C 3": 8.0}) == {}


def test_recovery_correction_divides_only_when_recovery_is_below_one():
    """EPA-1668: correct for losses, never deflate for apparent gains."""
    analytes = [analyte("1", homolog="mono"), analyte("2", homolog="di", q1=222.0)]
    masses = {"1": 4.0, "2": 4.0}
    corrected = correct_for_recovery(masses, {"mono": 0.5, "di": 1.25}, analytes)
    assert corrected["1"] == pytest.approx(8.0)
    assert corrected["2"] == pytest.approx(4.0)  # no inflation


def test_recovery_of_exactly_one_leaves_the_mass_alone():
    analytes = [analyte("1", homolog="mono")]
    assert correct_for_recovery({"1": 4.0}, {"mono": 1.0}, analytes)["1"] == pytest.approx(4.0)


def test_non_positive_or_missing_recovery_leaves_the_mass_uncorrected():
    """A recovery outside 0 to 1 is not applied, and neither is a missing one.

    Dividing by a zero or negative recovery would produce an infinity or a sign flip, and a
    homolog whose surrogate was never measured has no evidence to correct by. In both cases
    the uncorrected mass is reported and the sample is flagged, rather than a number that
    looks corrected and is not.
    """
    analytes = [analyte("1", homolog="mono")]
    for recoveries in ({"mono": 0.0}, {"mono": -0.4}, {}, {"hexa": 0.5}):
        assert correct_for_recovery({"1": 4.0}, recoveries, analytes)["1"] == pytest.approx(4.0)


def test_labeled_standards_are_not_recovery_corrected(mini_analytes):
    """Surrogates and internal standards carry no homolog, so no surrogate corrects them."""
    masses = {"13C 3": 6.35, "30D": IS_SPIKE_NG}
    corrected = correct_for_recovery(masses, {"mono": 0.5}, mini_analytes)
    assert corrected["13C 3"] == pytest.approx(6.35)
    assert corrected["30D"] == pytest.approx(IS_SPIKE_NG)


def test_recovery_correction_reproduces_the_workbook_corrected_mass():
    """Congener 1: uncorrected 3.493954, mono recovery 0.78716 -> corrected 4.438709.

    Both inputs are the workbook's *displayed* values, so they carry the workbook's display
    rounding and the quotient can only reproduce its corrected cell to about 2e-5. The
    tolerance is that limit, not a slack allowance: tightening it would be a test of Excel's
    number formatting rather than of this formula.
    """
    analytes = [analyte("1", homolog="mono")]
    corrected = correct_for_recovery({"1": 3.4939539}, {"mono": 0.7871585}, analytes)
    assert corrected["1"] == pytest.approx(4.438709, abs=1e-4)


# --------------------------------------------------------------------------------------
# LOQ and censoring
# --------------------------------------------------------------------------------------


def frame(data: dict[str, dict[str, float]]) -> pd.DataFrame:
    """{sample: {label: mass}} -> the label-by-sample frame the run uses."""
    return pd.DataFrame(data)


def test_loq_is_mean_plus_k_sd_of_the_blanks():
    blanks = frame({"LB1": {"1": 1.0}, "LB2": {"1": 2.0}, "LB3": {"1": 3.0}})
    loq = compute_loq(blanks)
    # sample SD (ddof=1) of 1,2,3 is exactly 1.0
    assert loq["1"] == pytest.approx(2.0 + 3.0 * 1.0)


def test_loq_k_is_parameterised():
    blanks = frame({"LB1": {"1": 1.0}, "LB2": {"1": 2.0}, "LB3": {"1": 3.0}})
    assert compute_loq(blanks, QaqcParams(k_sd=2.0))["1"] == pytest.approx(4.0)


def test_loq_needs_at_least_two_blanks():
    """One blank has no spread to measure, so no LOQ is claimed, per ``min_blanks_for_loq``."""
    assert compute_loq(frame({"LB1": {"1": 1.0}})).empty
    assert compute_loq(pd.DataFrame()).empty


def test_a_single_blank_gives_a_zero_width_loq_when_explicitly_allowed():
    """The ``ddof=1 if size > 1 else 0`` branch, reachable only by lowering the minimum.

    With one blank the population SD is 0, so LOQ is that blank's own mass. This is off by
    default because it is not a limit of quantification in any real sense.
    """
    loq = compute_loq(frame({"LB1": {"1": 1.5}}), QaqcParams(min_blanks_for_loq=1))
    assert loq["1"] == pytest.approx(1.5)


def test_loq_ignores_missing_cells_but_keeps_the_analyte():
    """A blank that never reported an analyte contributes a zero rather than being skipped."""
    blanks = pd.DataFrame({"LB1": {"1": 3.0, "2": 1.0}, "LB2": {"1": 3.0}})
    loq = compute_loq(blanks)
    assert loq["1"] == pytest.approx(3.0)
    assert loq["2"] == pytest.approx(0.5 + 3.0 * math.sqrt(0.5))


def test_mass_exactly_at_the_loq_is_censored():
    """The comparison is ≤: at the limit is not above it."""
    masses = frame({"S1": {"1": 5.0, "2": 5.000001}})
    censored = censor_below_loq(masses, pd.Series({"1": 5.0, "2": 5.0}))
    assert censored.loc["1", "S1"] == 0.0
    assert censored.loc["2", "S1"] == pytest.approx(5.000001)


def test_censoring_leaves_analytes_with_no_loq_alone():
    masses = frame({"S1": {"1": 5.0, "2": 0.4}})
    censored = censor_below_loq(masses, pd.Series({"1": 1.0}))
    assert censored.loc["2", "S1"] == pytest.approx(0.4)


def test_censoring_is_a_no_op_without_an_loq_or_when_switched_off():
    masses = frame({"S1": {"1": 5.0}})
    pd.testing.assert_frame_equal(censor_below_loq(masses, pd.Series(dtype=float)), masses)
    pd.testing.assert_frame_equal(
        censor_below_loq(masses, pd.Series({"1": 10.0}), QaqcParams(censor_below_loq=False)), masses
    )


def test_censoring_does_not_mutate_its_input():
    masses = frame({"S1": {"1": 0.5}})
    censor_below_loq(masses, pd.Series({"1": 1.0}))
    assert masses.loc["1", "S1"] == pytest.approx(0.5)


# --------------------------------------------------------------------------------------
# Flags
# --------------------------------------------------------------------------------------


def test_recovery_flags_fire_outside_fifty_to_one_hundred_fifty_percent():
    flags = recovery_flags({"mono": 0.49, "di": 0.50, "tetra": 1.50, "hexa": 1.51})
    assert any("low_recovery:mono" in f for f in flags)
    assert any("high_recovery:hexa" in f for f in flags)
    assert not any("di" in f or "tetra" in f for f in flags)  # the bounds themselves are in range


def test_low_confidence_count_uses_detection_not_signal_to_noise():
    """An absent peak routinely reports a huge S/N. Counting on ``sn`` would miss every one.

    ``core.peaks`` documents that ``sn`` is meaningless, and frequently inverted, for a peak
    the presence test rejected.
    """
    peaks = [
        peak("1", 100.0, absent=True),  # sn 50, not detected
        peak("2", 100.0, flags=["low_sn"]),
        peak("3", 100.0),
    ]
    flags = integration_flags(peaks)
    assert flags == ["1_undetected_analytes"]


def test_no_integration_flag_when_every_analyte_is_detected():
    assert integration_flags([peak("1", 100.0), peak("2", 100.0)]) == []


def test_sample_flags_combine_recovery_and_integration_verdicts():
    flags = sample_flags({"mono": 0.2}, [peak("1", 10.0, absent=True)])
    assert len(flags) == 2
    assert any(f.startswith("low_recovery") for f in flags)
    assert any(f.endswith("_undetected_analytes") for f in flags)


# --------------------------------------------------------------------------------------
# Reference bracketing (the recovery denominator)
# --------------------------------------------------------------------------------------


def when(hour: int) -> datetime:
    return datetime(2023, 10, 2, hour, tzinfo=timezone.utc)


def test_a_file_is_bracketed_by_the_reference_that_preceded_it():
    """Nearest *preceding* reference, which is what opening a batch with one means."""
    from quant209.core.pipeline import bracket_references

    acquired = {
        "Ref04": when(1), "S1": when(2), "S2": when(3),
        "Ref05": when(4), "S3": when(5), "S4": when(6),
    }
    pairing = bracket_references(acquired, ["Ref04", "Ref05"], ["S1", "S2", "S3", "S4"])
    assert pairing == {"S1": "Ref04", "S2": "Ref04", "S3": "Ref05", "S4": "Ref05"}


def test_bracketing_is_not_nearest_in_time():
    """The case the rule exists for, taken from a real sequence: a blank, ``LB``.

    It was acquired 14 h after ``Ref04`` and 7 h before ``Ref05``, and the analyst measured it
    against ``Ref04``. A bracketing standard opens a batch. It does not retroactively claim the
    batch before it, so proximity is the wrong question.
    """
    from quant209.core.pipeline import bracket_references

    acquired = {"Ref04": when(0), "LB": when(14), "Ref05": when(21)}
    assert bracket_references(acquired, ["Ref04", "Ref05"], ["LB"]) == {"LB": "Ref04"}


def test_a_file_before_every_reference_takes_the_earliest_one():
    from quant209.core.pipeline import bracket_references

    acquired = {"S1": when(1), "Ref04": when(2), "Ref05": when(3)}
    assert bracket_references(acquired, ["Ref05", "Ref04"], ["S1"]) == {"S1": "Ref04"}


def test_a_file_with_no_timestamp_falls_back_to_the_first_reference():
    """An unordered file cannot be bracketed. The pipeline reports the fallback in ``notes``."""
    from quant209.core.pipeline import bracket_references

    acquired = {"Ref04": when(1), "Ref05": when(4), "S1": None}
    assert bracket_references(acquired, ["Ref04", "Ref05"], ["S1"]) == {"S1": "Ref04"}


def test_bracketing_without_any_timestamps_degrades_to_the_legacy_single_reference():
    from quant209.core.pipeline import bracket_references

    acquired = {"Ref04": None, "Ref05": None, "S1": None}
    assert bracket_references(acquired, ["Ref04", "Ref05"], ["S1"]) == {"S1": "Ref04"}


def test_a_run_with_no_reference_injection_pairs_nothing():
    from quant209.core.pipeline import bracket_references

    assert bracket_references({"S1": when(1)}, [], ["S1"]) == {}


def test_an_unknown_reference_policy_is_refused(reference_run):
    from quant209.core.pipeline import run_analysis
    from quant209.core.sequence import find_data_files, scan_inputs

    entries = scan_inputs(find_data_files(reference_run))
    with pytest.raises(ValueError, match="reference policy"):
        run_analysis(entries, reference_run / "Converted", reference="nearest")


# --------------------------------------------------------------------------------------
# End to end
# --------------------------------------------------------------------------------------


def test_run_analysis_on_the_reference_run(reference_run):
    """``run_analysis`` on a real run: the shape and provenance a report needs."""
    from quant209.core.pipeline import run_analysis
    from quant209.core.sequence import find_data_files, scan_inputs

    entries = scan_inputs(find_data_files(reference_run))
    seen: list[tuple[str, int, int, str]] = []
    result = run_analysis(
        entries,
        reference_run / "Converted",
        params=IntegrationParams(),
        progress=lambda phase, i, n, label: seen.append((phase, i, n, label)),
    )

    # Provenance of the run itself.
    assert result.fit.fit_ok
    assert result.status == "OK"
    assert result.alignment.is_identity  # this run is the canonical GC method
    assert "Cal01" in result.cal_files[0]
    assert "Ref04" in result.ref_files[0]
    assert result.params == IntegrationParams()

    # The run is three bracketed batches (Cal01/Ref04, Cal02/Ref05, Cal03/Ref06, from the
    # instrument sequence), and every quantified file must be paired with the reference that
    # opened its own batch. The reported recovery is highly sensitive to that pairing, so it
    # is asserted by name and not merely by "some reference was used".
    assert result.reference_policy == "bracketed"
    assert set(result.reference_for) == set(result.sample_ids) | set(result.blank_ids)
    assert result.reference_for[NAMES.sample] == "Ref04"
    # 7 h nearer Ref05, but Ref04 brackets it:
    assert result.reference_for[NAMES.blank_late] == "Ref04"
    assert result.reference_for[NAMES.sample_second_batch] == "Ref05"
    assert result.reference_for[NAMES.sample_in_workbook_3] == "Ref06"
    assert result.files[NAMES.sample].acquired is not None

    # Shape: 173 native rows, every quantified file kept, hexane blanks excluded.
    assert len(result.sample_ids) >= 14
    assert result.blank_ids
    natives = [a.label for a in result.analytes if a.role == "native"]
    assert len(natives) == 173
    assert set(natives) <= set(result.masses.index)
    for frame_ in (result.areas, result.masses, result.masses_corrected, result.masses_censored):
        assert set(result.sample_ids) <= set(frame_.columns)
    assert not any(e.file_type == "hex" for e in result.quantified_entries)

    # Per-file audit records survive, including the analytes that never got integrated.
    for sample_id in result.sample_ids:
        record = result.files[sample_id]
        assert record.peaks
        assert record.mzml.suffix.lower() == ".mzml"
        assert isinstance(record.skipped_analytes, list)
    assert result.skipped_analytes == {
        name: rec.skipped_analytes for name, rec in result.files.items()
    }

    # Quantification actually happened.
    assert result.calibration.rrf
    assert result.recoveries.notna().to_numpy().any()
    assert (result.recoveries.loc["mono"].dropna() > 0).all()
    assert result.masses.loc[natives].to_numpy(dtype=float).sum() > 0
    assert result.loq.notna().any()

    # The IS quantifies to its own spike amount in every sample, by construction. This is the
    # cheapest end-to-end check that the mass chain is wired the right way round.
    is_masses = result.masses.loc["30D", result.sample_ids].astype(float)
    assert np.allclose(is_masses.dropna(), IS_SPIKE_NG, rtol=1e-9)

    # Native sample masses are censored. The scoping below must not turn censoring off.
    assert (result.masses_censored.loc[natives, result.sample_ids] == 0.0).to_numpy().any()

    # Progress reporting covers every phase and every file.
    assert {p for p, *_ in seen} >= {"calibrate", "integrate", "quantify"}
    integrated = [label for phase, _i, _n, label in seen if phase == "integrate"]
    assert len(integrated) == len(result.quantified_entries)


@pytest.fixture(scope="module")
def real_result():
    """One ``run_analysis`` over the reference run, shared by the four scope tests below.

    Module-scoped because a run is dozens of files and each of those tests interrogates a
    different corner of the same result. Re-integrating per test would quadruple this file's
    runtime for no extra coverage. It resolves the run folder itself rather than taking
    ``conftest``'s ``reference_run``, which is function-scoped and cannot be widened without
    changing the scope of every other test that uses it.
    """
    from quant209.core.pipeline import run_analysis
    from quant209.core.sequence import find_data_files, scan_inputs

    folder = run_dir()
    return run_analysis(scan_inputs(find_data_files(folder)), folder / "Converted")


def test_the_internal_standard_survives_censoring(real_result):
    """``masses_censored`` must still read ``IS_SPIKE_NG`` for 30D in every sample.

    A labelled standard is spiked at a known amount into every file, including the blanks, so
    its blank masses are identical. Their standard deviation is exactly 0, so
    ``LOQ = mean + 3*SD`` comes out at the spike amount itself, making mass ≤ LOQ true in
    every sample. Censoring the whole frame would zero out the one column that demonstrates
    the mass chain is wired the right way round.
    """
    assert real_result.loq["30D"] == pytest.approx(IS_SPIKE_NG)  # LOQ lands on the spike amount
    censored = real_result.masses_censored.loc["30D", real_result.sample_ids].astype(float)
    assert np.allclose(censored.dropna(), IS_SPIKE_NG, rtol=1e-9)
    assert (censored != 0.0).all()


def test_censoring_touches_only_native_rows_and_reported_files(real_result):
    """The LOQ is a question about congeners in samples, and is applied nowhere else.

    Calibration and reference injections are measurements the run is built on, such as
    response factors and the recovery denominator, not results to be reported against a
    blank level. Censoring them would contradict what the pipeline documents about them.
    """
    labeled = [a.label for a in real_result.analytes if a.role != "native"]
    unreported = real_result.cal_ids + real_result.ref_ids
    pd.testing.assert_frame_equal(
        real_result.masses_censored.loc[labeled], real_result.masses_corrected.loc[labeled]
    )
    pd.testing.assert_frame_equal(
        real_result.masses_censored[unreported], real_result.masses_corrected[unreported]
    )


def test_detection_travels_beside_the_mass_as_a_frame(real_result):
    """``RunResult.detected``: the presence verdict per (analyte, file), as a bool mask.

    Quantification does not apply it: the analyst quantifies a shoulder they chose to split
    off a neighbour, and so does this chain. It therefore has to reach the report as data, or
    the gate-on-``detected``-never-on-``sn`` rule would have nothing to gate with.
    """
    detected = real_result.detected
    assert detected.shape == real_result.areas.shape
    assert detected.dtypes.eq(bool).all()
    assert detected.loc["30D", real_result.sample_ids].all()
    assert not detected.to_numpy().all()  # some analyte, somewhere, was rejected
    for sample_id in real_result.sample_ids[:3]:
        expected = {p.label for p in real_result.files[sample_id].peaks if p.detected}
        column = detected[sample_id]
        assert set(column.index[column]) == expected


def test_an_analyte_with_no_peak_reads_undetected_rather_than_null(real_result):
    """``detected`` is a mask, so a missing peak is False and never NaN (three-valued logic
    in a mask invites bugs, and ``skipped_analytes`` keeps the distinction)."""
    assert not real_result.detected.isna().to_numpy().any()
    for sample_id in real_result.sample_ids:
        for label in real_result.files[sample_id].skipped_analytes:
            assert not real_result.detected.at[label, sample_id]


def test_run_analysis_requires_a_calibration_injection(reference_run):
    from quant209.core.pipeline import run_analysis
    from quant209.core.sequence import find_data_files, scan_inputs

    entries = [e for e in scan_inputs(find_data_files(reference_run)) if e.file_type != "cal"]
    with pytest.raises(ValueError, match="calibration"):
        run_analysis(entries, reference_run / "Converted")


# ---------------------------------------------------------------------------
# Calibration bracketing (retention-time drift across a sequence)
# ---------------------------------------------------------------------------


def test_nearest_calibration_pairs_each_file_with_its_closest_injection():
    """Drift is continuous, so the nearest calibration wins, on either side.

    This differs from ``bracket_references``, which takes the nearest *preceding* reference
    because a reference standard opens a batch and defines its recovery denominator. Nothing
    opens a batch here.
    """
    from datetime import datetime, timezone

    from quant209.core.pipeline import nearest_calibration

    def when(hour: int) -> datetime:
        return datetime(2024, 1, 3, hour, tzinfo=timezone.utc)

    acquired = {
        "Cal01": when(0), "Cal02": when(8), "Cal03": when(16),
        "early": when(1), "just_after_cal02": when(9), "just_before_cal03": when(15),
        "untimed": None,
    }
    pairing = nearest_calibration(acquired, ["Cal01", "Cal02", "Cal03"],
                                  ["early", "just_after_cal02", "just_before_cal03", "untimed"])
    assert pairing["early"] == "Cal01"
    assert pairing["just_after_cal02"] == "Cal02"
    # 15:00 is one hour before Cal03 and seven after Cal02: a *preceding*-only rule would
    # wrongly keep Cal02 here, which is the drift this pairing exists to remove.
    assert pairing["just_before_cal03"] == "Cal03"
    assert pairing["untimed"] == "Cal01", "an unordered file falls back to the first"


def test_nearest_calibration_handles_a_run_with_no_timestamps():
    from quant209.core.pipeline import nearest_calibration

    pairing = nearest_calibration({"a": None, "Cal01": None}, ["Cal01"], ["a"])
    assert pairing == {"a": "Cal01"}


def test_run_records_which_calibration_anchored_each_file(real_result):
    """The pairing is persisted as data: a reviewer must be able to see, for any file,
    which calibration injection supplied the retention times it was integrated at."""
    assert real_result.calibration_for, "no calibration pairing recorded"
    named = {rec.entry.name for rec in real_result.files.values()}
    assert set(real_result.calibration_for) <= named
    for file_id in real_result.reported_ids:
        name = real_result.files[file_id].entry.name
        assert name in real_result.calibration_for, f"{name} has no anchoring calibration"
        assert real_result.calibration_for[name] in real_result.cal_files


def test_an_internal_standard_that_disagrees_with_its_run_is_rejected():
    """The standard is spiked identically into every file, so one file disagreeing by an
    order of magnitude is a mis-integration, and it would scale that file's whole profile.

    Measured case: one file's PCB204 integrated to 66 counts·s against a run median of
    65,704, which multiplied its high-chlorination masses by roughly a thousand.
    """
    from quant209.core.pipeline import FileRecord, reject_inconsistent_standards
    from quant209.core.sequence import FileEntry
    from pathlib import Path

    def record(area: float) -> FileRecord:
        rec = FileRecord(entry=FileEntry(path=Path("x.mzML"), file_type="sample"),
                         mzml=Path("x.mzML"))
        rec.is_areas = {"204": area, "30D": 50_000.0}
        return rec

    records = {f"S{i}": record(65_000.0) for i in range(5)}
    records["bad"] = record(66.0)
    notes = reject_inconsistent_standards(records, [*records])

    assert "204" not in records["bad"].is_areas, "the outlier must not be used as a denominator"
    assert "30D" in records["bad"].is_areas, "the consistent standard must be left alone"
    for fid in records:
        if fid != "bad":
            assert "204" in records[fid].is_areas, f"{fid} was wrongly rejected"
    assert any("204" in n and "bad" in n for n in notes), "the rejection must be recorded"


def test_a_diluted_injection_keeps_its_internal_standards():
    """Both standards at an eighth of the run median, in step, is a diluted extract.

    Measured case: six archive samples whose two standards sat at 0.11 to 0.14 of their
    run's median, agreeing within 6%, with areas matching the analyst's. Every mass is an
    analyte-to-standard ratio, so the dilution cancels. Rejecting them reported 137k ng of
    analyst mass as unquantified. A standard that collapses on its own, or both collapsing to
    a thousandth, is still a failed integration and is still dropped.
    """
    from quant209.core.pipeline import FileRecord, reject_inconsistent_standards
    from quant209.core.sequence import FileEntry
    from pathlib import Path

    def record(is204: float, is30: float) -> FileRecord:
        rec = FileRecord(entry=FileEntry(path=Path("x.mzML"), file_type="sample"),
                         mzml=Path("x.mzML"))
        rec.is_areas = {"204": is204, "30D": is30}
        return rec

    records = {f"S{i}": record(65_000.0, 50_000.0) for i in range(5)}
    records["diluted"] = record(65_000.0 * 0.125, 50_000.0 * 0.13)
    records["lopsided"] = record(65_000.0 * 0.125, 50_000.0)
    records["failed"] = record(65.0, 50.0)
    notes = reject_inconsistent_standards(records, [*records])

    assert records["diluted"].is_areas == {"204": 8_125.0, "30D": 6_500.0}
    assert any("diluted" in n and "diluted or concentrated" in n for n in notes)
    assert "204" not in records["lopsided"].is_areas and "30D" in records["lopsided"].is_areas
    assert records["failed"].is_areas == {}, "a thousandfold collapse is not a dilution"


def test_internal_standard_check_needs_a_population():
    """With too few files there is no median to judge against, so nothing is rejected."""
    from quant209.core.pipeline import FileRecord, reject_inconsistent_standards
    from quant209.core.sequence import FileEntry
    from pathlib import Path

    rec = FileRecord(entry=FileEntry(path=Path("x.mzML"), file_type="sample"), mzml=Path("x.mzML"))
    rec.is_areas = {"204": 1.0}
    records = {"only": rec}
    assert reject_inconsistent_standards(records, ["only"]) == []
    assert records["only"].is_areas == {"204": 1.0}


def _record(label_rows, *, is_areas, window_max, usable=True):
    """A stand-in file record: ``label_rows`` is ``[(label, detected, height), ...]``."""
    from types import SimpleNamespace

    peaks = [SimpleNamespace(label=lab, detected=det, height=h) for lab, det, h in label_rows]
    return SimpleNamespace(peaks=peaks, is_areas=is_areas, is_window_max=window_max, usable=usable)


def test_a_missing_internal_standard_is_classified_by_what_its_channel_carries():
    """Measured cases. One archive injection has no internal standard at all: its 30D channel
    tops out at 152 counts within a minute of the expected time, against standard peaks near
    14,000 in the files around it. In another file the standard is there but eluted 0.16 min
    late and was not detected. Both get no mass; the reason says which remedy applies."""
    from quant209.core.pipeline import classify_missing_standards

    normal = {"30D": 60_000.0}
    records = {f"S{i}": _record([("30D", True, 14_000.0)], is_areas=dict(normal), window_max={"30D": 14_000.0})
               for i in range(4)}
    records["empty"] = _record([("30D", False, 150.0)], is_areas={}, window_max={"30D": 152.0})
    records["shifted"] = _record([("30D", False, 900.0)], is_areas={}, window_max={"30D": 2_415.0})
    records["refused"] = _record([("30D", True, 14_000.0)], is_areas={}, window_max={"30D": 14_000.0})
    records["undecoded"] = _record([], is_areas={}, window_max={}, usable=False)

    got = classify_missing_standards(records, list(records), ["30D"])

    assert got == {"empty": {"30D": "absent"}, "shifted": {"30D": "not_found"},
                   "refused": {"30D": "rejected"}}


def test_a_calibration_whose_retention_times_do_not_fit_anchors_nothing():
    """Measured case: a run's first calibration sat 4.3 min late against its other three. The
    four files nearest it were integrated at its retention times and lost both internal
    standards. It is excluded as long as another calibration fits; if none does, nothing is
    excluded and the run fails the fit check as before."""
    from types import SimpleNamespace

    from quant209.core.pipeline import select_calibrations

    fit = lambda ok, why="": SimpleNamespace(fit_ok=ok, reason=why)
    order = ["Cal01", "Cal02", "Cal03", "Cal04"]

    usable, excluded = select_calibrations(
        {"Cal01": fit(False, "only 31% of the native analytes were found"), "Cal02": fit(True),
         "Cal03": fit(True), "Cal04": fit(True)}, order)
    assert usable == ["Cal02", "Cal03", "Cal04"]
    assert excluded == {"Cal01": "only 31% of the native analytes were found"}

    assert select_calibrations({c: fit(True) for c in order}, order) == (order, {})
    assert select_calibrations({c: fit(False, "x") for c in order}, order) == (order, {})


def _shifted_file(shifts: dict[str, float], with_204: bool = False,
                  heights: dict[str, float] | None = None):
    """Synthetic file: one peak per channel, each moved by ``shifts`` from its anchor.

    ``with_204`` adds an internal standard 204, on the channel of native 52 as in the real
    canonical table. ``heights`` overrides a peak's height (default 10,000 counts).
    """
    from types import SimpleNamespace

    from quant209.core.config import Analyte
    from quant209.core.integrate import integrate_run
    from quant209.core.mzml import Chromatogram
    from quant209.core.peaks import IntegrationParams

    spec = [("13C 15", 235.0, 165.0, "surrogate_13C", 17.47), ("30D", 304.0, 234.0, "IS_d30", 16.19),
            ("13C 28", 269.0, 199.0, "surrogate_13C", 20.84), ("18", 256.0, 186.0, "native", 18.00),
            ("52", 290.0, 220.0, "native", 30.00), ("13C 118", 338.0, 268.0, "surrogate_13C", 40.30)]
    if with_204:
        spec.append(("204", 290.0, 220.0, "IS_204", 32.00))
    analytes = [Analyte(l, q1, q3, "tri" if role == "native" else "", 2 if l == "204" else 1, role,
                        "30D", 1.0, rt)
                for l, q1, q3, role, rt in spec]
    times = np.arange(10.0, 45.0, 0.005)
    rng = np.random.default_rng(3)
    traces: dict[tuple[float, float], np.ndarray] = {}
    for label, q1, q3, _role, rt in spec:
        if (q1, q3) not in traces:
            traces[(q1, q3)] = 300.0 + rng.normal(0.0, 30.0, times.size)
        height = (heights or {}).get(label, 1.0e4)
        traces[(q1, q3)] = traces[(q1, q3)] + height * np.exp(
            -0.5 * ((times - (rt + shifts.get(label, 0.0))) / 0.03) ** 2)
    chroms = {key: Chromatogram(q1=key[0], q3=key[1], times=times, intensities=y)
              for key, y in traces.items()}
    anchors = {a.label: a.cal_rt for a in analytes}
    params = IntegrationParams()
    record = SimpleNamespace(peaks=integrate_run(chroms, analytes, anchors, params),
                             entry=SimpleNamespace(name="AN_DEMO05.mzML"), is_areas={}, is_window_max={})
    return chroms, record, analytes, anchors, params


def _realign(file, lost=None):
    """Run the realignment on a synthetic file; ``lost`` defaults to the standards not detected."""
    from quant209.core.pipeline import _realign_with, landmark_labels

    chroms, record, analytes, anchors, params = file
    landmarks = landmark_labels(analytes)
    if lost is None:
        lost = [p.label for p in record.peaks if p.label in landmarks and not p.detected]
    return _realign_with(chroms, record, analytes, anchors, params,
                         {l: 1.0e4 for l in landmarks}, landmarks, lost)


def test_a_file_whose_standards_moved_together_is_searched_where_they_went():
    """Measured case: in one archive file everything from 16 to 24 min eluted 0.2 min late. The
    internal standard 30D and three surrogates were lost, and a hundred congeners with them."""
    from quant209.core.pipeline import REALIGNED_FLAG, landmark_labels

    late = {"13C 15": 0.2, "30D": 0.2, "13C 28": 0.2, "18": 0.2}
    file = _shifted_file(late)
    chroms, record, analytes, anchors, params = file
    before = {p.label: p.detected for p in record.peaks}
    assert not any(before[l] for l in late), "the failure: everything that moved is lost"
    assert landmark_labels(analytes) == {"13C 15", "30D", "13C 28", "13C 118"}

    note = _realign(file)

    assert note and "+0.20 min" in note
    after = {p.label: p for p in record.peaks}
    for label in late:
        assert after[label].detected and REALIGNED_FLAG in after[label].flags, label
        assert abs(after[label].rt_apex - (anchors[label] + 0.2)) < 0.02
    # A standard found where it was expected, and a native that never moved, are still found
    # where they are: the drift curve tapers from the moved standards to the unmoved one, and
    # only the unmoved standard itself, at zero offset, is left unflagged.
    assert after["13C 118"].detected and REALIGNED_FLAG not in after["13C 118"].flags
    assert after["52"].detected and abs(after["52"].rt_apex - anchors["52"]) < 0.02
    assert "30D" in record.is_areas


def test_standards_that_disagree_about_the_shift_leave_the_file_alone():
    file = _shifted_file({"13C 15": 0.2, "30D": -0.2, "13C 28": 0.2})
    record = file[1]
    original = list(record.peaks)
    assert _realign(file) is None
    assert record.peaks == original


def test_a_drift_that_grows_along_the_run_is_followed_by_a_curve():
    """Measured case: in two concentrated archive runs the standards sat 0.05 min late at 23 min
    and 0.33 min late at 67 min. Standards 0.18 min apart in offset do not agree on one shift,
    but they lie on one gentle line, and each native is searched where that line puts it."""
    from quant209.core.pipeline import REALIGNED_FLAG

    drift = {"30D": 0.20, "13C 15": 0.22, "13C 28": 0.27, "13C 118": 0.38}
    natives = {"18": 0.228, "52": 0.322}   # the same line, read at their retention times
    file = _shifted_file({**drift, **natives})
    chroms, record, analytes, anchors, params = file
    assert not any(p.detected for p in record.peaks if p.label in {**drift, **natives}), \
        "the failure: everything that drifted is lost"

    note = _realign(file)

    assert note
    after = {p.label: p for p in record.peaks}
    for label, offset in {**drift, **natives}.items():
        assert after[label].detected and REALIGNED_FLAG in after[label].flags, label
        assert abs(after[label].rt_apex - (anchors[label] + offset)) < 0.02, label


def test_a_shift_the_natives_do_not_share_leaves_the_file_alone():
    """Measured case: two standards sat 0.15 min from their anchors in every file of a run, an
    offset the native congeners did not share. Searching where the standards were found left
    the natives exactly where they had been, and the file was flagged for nothing. A shift is
    kept only if the natives it moves the search for are found in greater number."""
    file = _shifted_file({"13C 15": 0.22, "30D": 0.22, "13C 28": 0.22, "18": 0.08, "52": 0.08})
    record = file[1]
    original = list(record.peaks)
    assert all(p.detected for p in original if p.label in ("18", "52")), "the natives are found"

    assert _realign(file) is None
    assert record.peaks == original


def test_the_standard_search_ignores_a_peak_far_above_its_usual_height():
    """Measured case: a standard's channel held a peak four times its usual height among dozens
    of other maxima, and two such channels agreed on a shift by chance. A standard is spiked at
    one amount, so a peak several times taller is something else, and the standard is left out."""
    from quant209.core.pipeline import _lost_standard_offsets

    chroms, record, analytes, anchors, params = _shifted_file(
        {"13C 15": 0.25, "30D": 0.25, "13C 28": 0.25}, heights={"13C 15": 5.0e4, "13C 28": 1.0e3})
    typical = {"13C 15": 1.0e4, "30D": 1.0e4, "13C 28": 1.0e4}

    found = _lost_standard_offsets(chroms, analytes, anchors, params, typical,
                                   ["13C 15", "30D", "13C 28"])

    # 13C 15 is five times its usual height; 13C 28 is a tenth of it (below the 15% floor).
    assert set(found) == {"30D"}
    assert abs(found["30D"] - 0.25) < 0.01


def test_a_lost_internal_standard_in_a_file_whose_standards_show_a_drift_is_searched_for():
    """Measured case: a whole file sat 0.1 to 0.2 min late. Every labelled standard on its own
    channel was still found, but the internal standard 204, on a crowded channel, was not, and
    nor were the natives near it, so the file got no masses."""
    from quant209.core.pipeline import REALIGNED_FLAG, _internal_standard_misplaced, landmark_labels

    drift = {"30D": 0.12, "13C 15": 0.12, "13C 28": 0.12, "13C 118": 0.12}
    file = _shifted_file({**drift, "18": 0.12, "52": 0.30, "204": 0.30}, with_204=True)
    chroms, record, analytes, anchors, params = file
    before = {p.label: p.detected for p in record.peaks}
    assert all(before[l] for l in drift), "the standards on their own channels are found"
    assert not before["204"] and not before["52"], "the failure: the crowded channel is lost"
    landmarks = landmark_labels(analytes)
    assert _internal_standard_misplaced(record, analytes, anchors, params, landmarks)

    note = _realign(file, lost=[])

    assert note
    after = {p.label: p for p in record.peaks}
    for label, offset in (("204", 0.30), ("52", 0.30)):
        assert after[label].detected and REALIGNED_FLAG in after[label].flags, label
        assert abs(after[label].rt_apex - (anchors[label] + offset)) < 0.02, label
    assert "204" in record.is_areas


def test_a_lost_internal_standard_alone_does_not_start_a_realignment():
    """Without a drift in the standards that were found, a missing 204 says nothing about where
    to look: it may simply not be in the injection."""
    from quant209.core.pipeline import _internal_standard_misplaced, landmark_labels

    chroms, record, analytes, anchors, params = _shifted_file({"204": 0.30}, with_204=True)
    assert not {p.label: p.detected for p in record.peaks}["204"]
    landmarks = landmark_labels(analytes)
    assert not _internal_standard_misplaced(record, analytes, anchors, params, landmarks)


def test_two_standards_with_nothing_else_to_anchor_the_curve_leave_the_file_alone():
    """A drift curve needs three points. Two standards found, with no other standard seen
    anywhere, would otherwise be extended as one shift over the whole chromatogram, and that
    would be accepted here because the natives did move."""
    # 13C 28 and 13C 118 moved beyond the search, so only 30D and 13C 15 are found.
    file = _shifted_file({"30D": 0.2, "13C 15": 0.2, "13C 28": 1.5, "13C 118": 1.5,
                          "18": 0.2, "52": 0.2})
    record = file[1]
    original = list(record.peaks)
    assert not any(p.detected for p in original if p.label in ("30D", "13C 15", "13C 28", "13C 118"))

    assert _realign(file) is None
    assert record.peaks == original


def test_the_note_of_a_sloped_drift_gives_the_range_of_offsets_it_followed():
    drift = {"30D": 0.20, "13C 15": 0.22, "13C 28": 0.27, "13C 118": 0.38}
    file = _shifted_file({**drift, "18": 0.228, "52": 0.322})

    note = _realign(file)

    assert note and "+0.20 to +0.38 min" in note
    assert "analytes from " in note and "were searched where the drift" in note
    assert "Native peaks found there: 0 before, 2 after" in note


def test_a_file_that_cannot_be_read_again_is_left_alone(tmp_path):
    from quant209.core.pipeline import landmark_labels, realign_shifted_file

    chroms, record, analytes, anchors, params = _shifted_file(
        {"13C 15": 0.2, "30D": 0.2, "13C 28": 0.2, "18": 0.2})
    record.mzml = tmp_path / "missing.mzML"
    landmarks = landmark_labels(analytes)
    original = list(record.peaks)

    assert realign_shifted_file(record, analytes, anchors, params,
                                {l: 1.0e4 for l in landmarks}, landmarks) is None
    assert record.peaks == original
