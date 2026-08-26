"""Tests for quant209.report.validation: the validation report generator.

The suite builds one small synthetic ``RunResult``/``RunTruth`` pair (2 analytes, 2 files:
one sample, one blank) rather than reading real run data, so it stays fast and deterministic
while still exercising the whole chain the real report generator drives: area join -> mass
join -> recovery join -> RRF join -> bias diagnosis -> CSVs -> figures -> summary.md.

Every number asserted against ``summary.md`` is recomputed independently from the CSV files
the same call wrote, using plain pandas/numpy rather than by importing the module's own
formatting helpers. The point is to catch summary.md drifting from the CSVs it is supposed
to describe, not to re-run the implementation against itself.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from quant209.core.config import Analyte
from quant209.core.peaks import IntegrationParams, PeakResult
from quant209.core.pipeline import FileRecord, RunResult
from quant209.core.qaqc import QaqcParams
from quant209.core.quantify import CalibrationModel
from quant209.core.sequence import FileEntry
from quant209.report.truth import RunTruth
from quant209.report.validation import (
    BAND,
    DEFAULT_MIN_AREA,
    DEFAULT_MIN_MASS,
    generate_validation_report,
    join_areas,
    join_masses,
)

LABEL_1 = "1"
LABEL_5 = "5"


def _analytes() -> list[Analyte]:
    return [
        Analyte(
            label=LABEL_1, q1=100.0, q3=50.0, homolog="mono", order=1, role="native",
            is_norm="30D", cal_mass=10.0, cal_rt=5.0,
        ),
        Analyte(
            label=LABEL_5, q1=200.0, q3=60.0, homolog="di", order=1, role="native",
            is_norm="30D", cal_mass=10.0, cal_rt=6.0,
        ),
    ]


def _peak(label: str, area: float, flags: list[str] | None = None) -> PeakResult:
    return PeakResult(
        label=label, q1=100.0, q3=50.0, rt_apex=5.0, rt_start=4.8, rt_end=5.2,
        area=area, height=area / 4, baseline_start=4.7, baseline_end=5.3,
        sn=30.0, method="valley_drop", flags=list(flags or []),
    )


@pytest.fixture
def analytes() -> list[Analyte]:
    return _analytes()


@pytest.fixture
def result(analytes: list[Analyte]) -> RunResult:
    """A tiny fake RunResult: 2 native analytes, 1 sample ('S1') + 1 blank ('S2').

    Areas/masses/RRFs are hand-picked so the joined tables land on both sides of the
    detectability threshold and the ±25% band. See the module docstring's numbers,
    reproduced in the tests below rather than hardcoded twice.
    """
    s1_peaks = [_peak(LABEL_1, 10000.0), _peak(LABEL_5, 20000.0, flags=["clustered"])]
    s2_peaks = [_peak(LABEL_1, 500.0), _peak(LABEL_5, 800.0)]

    files = {
        "S1": FileRecord(
            entry=FileEntry(path=Path("S1.mzML"), file_type="sample"),
            mzml=Path("S1.mzML"), peaks=s1_peaks,
        ),
        "S2": FileRecord(
            entry=FileEntry(path=Path("S2.mzML"), file_type="blank"),
            mzml=Path("S2.mzML"), peaks=s2_peaks,
        ),
    }

    areas = pd.DataFrame(
        {"S1": {LABEL_1: 10000.0, LABEL_5: 20000.0}, "S2": {LABEL_1: 500.0, LABEL_5: 800.0}}
    )
    masses = pd.DataFrame(
        {"S1": {LABEL_1: 4.8, LABEL_5: 5.5}, "S2": {LABEL_1: 0.04, LABEL_5: 0.03}}
    )
    detected = pd.DataFrame(
        {"S1": {LABEL_1: True, LABEL_5: True}, "S2": {LABEL_1: True, LABEL_5: True}}
    )
    recoveries = pd.DataFrame(
        {"S1": {"mono": 0.81, "di": 0.93}, "S2": {"mono": 0.93, "di": 0.96}}
    )
    calibration = CalibrationModel(
        rrf={LABEL_1: 2.0, LABEL_5: 1.5},
        cal_area={LABEL_1: 100000.0, LABEL_5: 50000.0},
        is_response={"30D": 1.0},
        invalid={},
        injections=(),
    )

    return RunResult(
        entries=[],
        quantified_entries=[],
        analytes=analytes,
        params=IntegrationParams(),
        qaqc=QaqcParams(),
        cal_files=[],
        ref_files=[],
        alignment=None,  # type: ignore[arg-type] (unused by anything under test)
        anchors={},
        fit=None,  # type: ignore[arg-type] (unused by anything under test)
        calibration=calibration,
        files=files,
        sample_ids=["S1"],
        blank_ids=["S2"],
        cal_ids=[],
        ref_ids=[],
        areas=areas,
        detected=detected,
        masses=masses,
        masses_corrected=masses.copy(),
        masses_censored=masses.copy(),
        recoveries=recoveries,
    )


@pytest.fixture
def truth() -> RunTruth:
    areas = pd.DataFrame(
        {"S1": {LABEL_1: 9800.0, LABEL_5: 25000.0}, "S2": {LABEL_1: 400.0, LABEL_5: 600.0}}
    )
    masses = pd.DataFrame(
        {"S1": {LABEL_1: 5.0, LABEL_5: 8.0}, "S2": {LABEL_1: 0.05, LABEL_5: 0.02}}
    )
    recoveries = pd.DataFrame(
        {"S1": {"mono": 0.80, "di": 0.95}, "S2": {"mono": 0.90, "di": 1.00}}
    )
    cal_rrfs = pd.DataFrame({"wb1": {LABEL_1: 2.05, LABEL_5: 2.4}})
    cal_areas = pd.DataFrame({"wb1": {LABEL_1: 100500.0, LABEL_5: 49500.0}})
    return RunTruth(
        areas=areas, masses=masses, recoveries=recoveries,
        cal_areas=cal_areas, cal_rrfs=cal_rrfs,
    )


@pytest.fixture
def report(result: RunResult, truth: RunTruth, tmp_path: Path):
    return generate_validation_report(result, truth, tmp_path, "TESTRUN")


# --------------------------------------------------------------------------------------
# Files exist
# --------------------------------------------------------------------------------------


def test_writes_every_expected_file(report, tmp_path: Path):
    expected = [
        "summary.md", "area_agreement.csv", "mass_agreement.csv", "recovery_agreement.csv",
        "rrf_agreement.csv", "bias_diagnosis.csv", "rpd_hist.pdf", "regression.pdf",
        "homolog_agreement.pdf",
    ]
    for name in expected:
        path = tmp_path / name
        assert path.exists(), f"{name} was not written"
        assert path.stat().st_size > 0, f"{name} is empty"


def test_figures_are_valid_nonempty_pdfs(report, tmp_path: Path):
    for name in ("rpd_hist.pdf", "regression.pdf", "homolog_agreement.pdf"):
        data = (tmp_path / name).read_bytes()
        assert data[:4] == b"%PDF", f"{name} does not start with a PDF header"
        assert len(data) > 200, f"{name} looks truncated ({len(data)} bytes)"


# --------------------------------------------------------------------------------------
# CSV row counts match the joined table sizes
# --------------------------------------------------------------------------------------


def test_area_csv_row_count_matches_join_areas(report, result: RunResult, truth: RunTruth, tmp_path: Path):
    peaks_by_file = {fid: result.files[fid].peaks for fid in result.reported_ids}
    file_types = {fid: result.files[fid].file_type for fid in result.reported_ids}
    rows = join_areas(peaks_by_file, file_types, result.analytes, truth.areas)
    df = pd.read_csv(tmp_path / "area_agreement.csv")
    assert len(df) == len(rows) == 4  # 2 analytes x 2 files


def test_mass_csv_row_count_matches_join_masses(report, result: RunResult, truth: RunTruth, tmp_path: Path):
    rows, _, _ = join_masses(result, truth.masses)
    df = pd.read_csv(tmp_path / "mass_agreement.csv")
    assert len(df) == len(rows) == 4


def test_recovery_csv_has_one_row_per_sample_homolog_pair(report, tmp_path: Path):
    df = pd.read_csv(tmp_path / "recovery_agreement.csv")
    assert len(df) == 4  # 2 samples x 2 homologs


def test_rrf_csv_has_one_row_per_native_analyte(report, tmp_path: Path):
    df = pd.read_csv(tmp_path / "rrf_agreement.csv")
    assert len(df) == 2


def test_bias_diagnosis_contains_only_the_out_of_band_detectable_row(report, tmp_path: Path):
    """Only S1/label 5 is native, detectable (manual mass 8.0 ≥ 1.0 ng) and outside the ±25%
    band (auto 5.5 vs manual 8.0 -> RPD ~ -37%). S1/label 1 is inside the band, and both S2
    rows are below the detectability threshold, so none of those three belong here."""
    df = pd.read_csv(tmp_path / "bias_diagnosis.csv", dtype={"label": str})
    assert len(df) == 1
    row = df.iloc[0]
    assert row["label"] == LABEL_5
    assert row["sample"] == "S1"
    assert bool(row["clustered"]) is True
    assert row["likely_cause"] == "cal-anomaly"  # this label's RRF itself disagrees > 25%


# --------------------------------------------------------------------------------------
# summary.md reproduces numbers independently recomputed from the CSVs
# --------------------------------------------------------------------------------------


def test_summary_area_headline_matches_csv(report, tmp_path: Path):
    df = pd.read_csv(tmp_path / "area_agreement.csv")
    native = df[df["role"] == "native"]
    detectable = native[native["manual"] >= DEFAULT_MIN_AREA]
    assert len(detectable) == 2  # S1's two rows. S2's are both sub-threshold.
    median_abs_rpd = float(np.median(np.abs(detectable["rpd"])))
    within_band = float(np.mean(np.abs(detectable["rpd"]) <= BAND))

    summary = (tmp_path / "summary.md").read_text(encoding="utf-8")
    assert f"{100 * median_abs_rpd:.1f}%" in summary
    assert f"{100 * within_band:.1f}%" in summary


def test_summary_mass_headline_matches_csv(report, tmp_path: Path):
    df = pd.read_csv(tmp_path / "mass_agreement.csv")
    native = df[df["role"] == "native"]
    detectable = native[native["manual"] >= DEFAULT_MIN_MASS]
    assert len(detectable) == 2  # S1's two rows. S2's are both sub-threshold.
    median_abs_rpd = float(np.median(np.abs(detectable["rpd"])))
    within_band = float(np.mean(np.abs(detectable["rpd"]) <= BAND))

    summary = (tmp_path / "summary.md").read_text(encoding="utf-8")
    assert f"{100 * median_abs_rpd:.1f}%" in summary
    assert f"{100 * within_band:.1f}%" in summary


def test_summary_recovery_headline_matches_csv(report, tmp_path: Path):
    df = pd.read_csv(tmp_path / "recovery_agreement.csv")
    mae = float(np.mean(np.abs(df["error"])))
    summary = (tmp_path / "summary.md").read_text(encoding="utf-8")
    assert f"{mae:.4f}" in summary


def test_summary_rrf_headline_matches_csv(report, tmp_path: Path):
    df = pd.read_csv(tmp_path / "rrf_agreement.csv")
    median_abs_rpd = float(np.median(np.abs(df["rpd"])))
    summary = (tmp_path / "summary.md").read_text(encoding="utf-8")
    assert f"{100 * median_abs_rpd:.1f}%" in summary


def test_summary_reports_a_verdict(report, tmp_path: Path):
    summary = (tmp_path / "summary.md").read_text(encoding="utf-8")
    assert "Verdict" in summary
    assert "PASS" in summary or "FAIL" in summary


def test_report_object_reflects_the_same_bias_row_count(report):
    assert report.n_bias_rows == 1
