"""Tests for quant209.report.excel: the analyst-format Excel workbook writer.

Most of this suite exercises a small, hand-built ``RunResult`` (the pattern used by
``tests/test_validation.py``): seven analytes, five files (one cal, one ref, two samples,
one blank), so structural and value correctness is checked fast and deterministically without
touching the real ~2 GB dataset.

Three tests are different on purpose. They read a real analyst workbook directly with
openpyxl and compare it, cell-for-cell address, against a report generated from a real
``run_analysis`` over the same run folder. The positional-fidelity claim this module makes
can only be checked against the genuine article. That fixture (``real_result``) is
module-scoped, like ``tests/test_quantify.py``'s, and skips cleanly when no instrument-data
archive is configured (see ``tests/dataset.py``).
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
import pytest

from quant209.core.calibrate import RtAlignment, RunFitReport
from quant209.core.config import Analyte
from quant209.core.peaks import IntegrationParams, PeakResult
from quant209.core.pipeline import STATUS_OK, FileRecord, RunResult
from quant209.core.qaqc import QaqcParams
from quant209.core.quantify import CalibrationModel
from quant209.core.sequence import FileEntry
from quant209.core.standards import HOMOLOG_SURROGATE, IS_SPIKE_NG
from quant209.report.excel import (
    CAL_COL_AREA,
    CAL_COL_INVALID,
    CAL_COL_LABEL,
    CAL_COL_RRF,
    CAL_DATA_START_ROW,
    CAL_HEADER_ROW,
    COL_AREA_DUP,
    COL_CAL_RT,
    COL_CORRECTED,
    COL_IS,
    COL_LABEL,
    COL_MASS,
    COL_PEAK,
    COL_RRF,
    SAMPLE_DATA_START_ROW,
    SAMPLE_HEADER_ROW,
    SUMMARY_DATA_START_ROW,
    SUMMARY_HEADER_ROW,
    write_report,
)

from dataset import NAMES, run_dir, workbook

# --------------------------------------------------------------------------------------
# A small, hand-built RunResult: 3 natives (mono/di/hexa, split across both internal
# standards), their 2 surrogates, and the 2 internal standards themselves.
# --------------------------------------------------------------------------------------

ANALYTES = [
    Analyte(label="1", q1=188.0, q3=153.1, homolog="mono", order=1, role="native",
            is_norm="30D", cal_mass=29.06, cal_rt=19.563),
    Analyte(label="13C 3", q1=200.0, q3=165.1, homolog="", order=1, role="surrogate_13C",
            is_norm="30D", cal_mass=61.86, cal_rt=22.3),
    Analyte(label="4", q1=222.0, q3=152.1, homolog="di", order=1, role="native",
            is_norm="30D", cal_mass=29.12, cal_rt=22.633),
    Analyte(label="153", q1=360.0, q3=253.0, homolog="hexa", order=1, role="native",
            is_norm="204", cal_mass=30.5, cal_rt=45.0),
    Analyte(label="13C 153", q1=372.0, q3=265.0, homolog="", order=1, role="surrogate_13C",
            is_norm="204", cal_mass=28.9, cal_rt=44.5),
    Analyte(label="30D", q1=261.0, q3=191.1, homolog="", order=1, role="IS_d30",
            is_norm="30D", cal_mass=34.92, cal_rt=28.3),
    Analyte(label="204", q1=430.0, q3=360.0, homolog="", order=1, role="IS_204",
            is_norm="204", cal_mass=34.92, cal_rt=54.5),
]
LABELS = [a.label for a in ANALYTES]


def _peak(label: str, area: float, rt: float, *, height: float | None = None,
          flags: list[str] | None = None) -> PeakResult:
    return PeakResult(
        label=label, q1=0.0, q3=0.0, rt_apex=rt, rt_start=rt - 0.05, rt_end=rt + 0.05,
        area=area, height=height if height is not None else area / 8.0,
        baseline_start=0.0, baseline_end=0.0, sn=40.0, method="valley_drop",
        flags=list(flags or []),
    )


def _entry(file_type: str, name: str) -> FileEntry:
    return FileEntry(path=Path(f"{name}.mzML"), file_type=file_type)


def _record(file_type: str, name: str, peaks: list[PeakResult],
            skipped: list[str] | None = None) -> FileRecord:
    return FileRecord(
        entry=_entry(file_type, name),
        mzml=Path(f"{name}.mzML"),
        peaks=peaks,
        skipped_analytes=list(skipped or []),
        is_areas={},
        acquired=datetime(2023, 10, 2, 8, tzinfo=timezone.utc),
    )


def _frame(data: dict[str, dict[str, float]], index: list[str], columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {c: pd.Series(data.get(c, {}), dtype=float) for c in columns}
    ).reindex(index=index, columns=columns)


@pytest.fixture(scope="module")
def mini_result() -> RunResult:
    all_labels_peaks = {
        "Cal01": [_peak(l, 1.0e5, 20.0 + i) for i, l in enumerate(LABELS)],
        "Ref01": [_peak(l, 8.0e4, 20.0 + i) for i, l in enumerate(LABELS)],
        "S1": [_peak(l, 5.0e4, 20.0 + i) for i, l in enumerate(LABELS)],
        # S2 has no peak at all for "153": a skipped analyte, exercising blank cells.
        "S2": [_peak(l, 4.5e4, 20.0 + i) for i, l in enumerate(LABELS) if l != "153"],
        "LB1": [_peak(l, 1.0e3, 20.0 + i) for i, l in enumerate(LABELS)],
    }
    files = {
        "Cal01": _record("cal", "Cal01", all_labels_peaks["Cal01"]),
        "Ref01": _record("ref", "Ref01", all_labels_peaks["Ref01"]),
        "S1": _record("sample", "S1", all_labels_peaks["S1"]),
        "S2": _record("sample", "S2", all_labels_peaks["S2"], skipped=["153"]),
        "LB1": _record("blank", "LB1", all_labels_peaks["LB1"]),
    }

    areas = _frame(
        {fid: {p.label: p.area for p in peaks} for fid, peaks in all_labels_peaks.items()},
        LABELS, list(files),
    )
    detected = pd.DataFrame(True, index=LABELS, columns=list(files))
    detected.loc["153", "S2"] = False  # matches the missing peak

    # masses: made up, but internally consistent enough to exercise every column.
    masses = areas * 2.0e-5
    masses.loc[["30D", "204"]] = IS_SPIKE_NG  # the IS quantifies to its own spike, by construction
    masses_corrected = masses.copy()
    masses_corrected.loc["1"] *= 1.2  # pretend a mono recovery < 1 corrected it up
    masses_censored = masses_corrected.copy()
    masses_censored.loc["1", "LB1"] = 0.0  # the one deliberately censored cell

    recoveries = pd.DataFrame(
        index=list(HOMOLOG_SURROGATE), columns=["S1", "S2", "LB1"], dtype=float
    )
    recoveries.loc["mono"] = [0.85, 0.80, 0.90]
    recoveries.loc["hexa"] = [0.55, 0.60, 0.65]

    loq = pd.Series({l: 500.0 for l in LABELS})

    calibration = CalibrationModel(
        rrf={"1": 2.2094, "4": 1.4762, "153": 1.05, "13C 3": 2.1645, "13C 153": 1.0,
             "30D": 1.0, "204": 1.0},
        cal_area={l: 1.0e5 for l in LABELS},
        is_response={"30D": IS_SPIKE_NG / 5.0e5, "204": IS_SPIKE_NG / 2.8e4},
        invalid={},
        injections=("Cal01",),
        area_floor=500.0,
    )

    return RunResult(
        entries=[_entry(f.file_type, fid) for fid, f in files.items()],
        quantified_entries=[_entry(f.file_type, fid) for fid, f in files.items()],
        analytes=list(ANALYTES),
        params=IntegrationParams(),
        qaqc=QaqcParams(),
        cal_files=["Cal01"],
        ref_files=["Ref01"],
        alignment=RtAlignment(),
        anchors={a.label: a.cal_rt for a in ANALYTES},
        fit=RunFitReport(fit_ok=True, native_found_frac=1.0, labeled_found_frac=1.0,
                          n_native=3, n_native_found=3, n_labeled=4, n_labeled_found=4),
        calibration=calibration,
        status=STATUS_OK,
        notes=["one usable blank: LOQ is illustrative only in this fixture"],
        files=files,
        sample_ids=["S1", "S2"],
        blank_ids=["LB1"],
        cal_ids=["Cal01"],
        ref_ids=["Ref01"],
        areas=areas,
        detected=detected,
        masses=masses,
        masses_corrected=masses_corrected,
        masses_censored=masses_censored,
        recoveries=recoveries,
        loq=loq,
        flags={"S1": [], "S2": ["low_recovery:hexa=60%"], "LB1": [], "Cal01": [], "Ref01": []},
        reference_policy="bracketed",
        reference_for={"S1": "Ref01", "S2": "Ref01", "LB1": "Ref01"},
    )


@pytest.fixture(scope="module")
def written_path(tmp_path_factory, mini_result: RunResult) -> Path:
    out_dir = tmp_path_factory.mktemp("excel_report")
    return write_report(mini_result, out_dir / "mini_report.xlsx", run_name="mini-test-run")


@pytest.fixture(scope="module")
def wb(written_path: Path):
    workbook = openpyxl.load_workbook(written_path)
    yield workbook
    workbook.close()


# --------------------------------------------------------------------------------------
# Sheet set
# --------------------------------------------------------------------------------------


def test_sheet_set_mirrors_the_analyst_plus_documented_additions(wb, mini_result: RunResult):
    expected = set(mini_result.reported_ids) | {
        "cal", "ref", "Summary", "SamplesList_Standards", "run_info",
    }
    assert set(wb.sheetnames) == expected


def test_cal_and_ref_injections_get_no_per_file_sheet(wb):
    """Only reported files (samples, blanks) get a per-file sheet. Cal/ref are summarized
    on their own sheets instead, matching the analyst's own workbook, which holds no
    per-sample-shaped sheet for Cal01/Ref01 either."""
    assert "Cal01" not in wb.sheetnames
    assert "Ref01" not in wb.sheetnames


# --------------------------------------------------------------------------------------
# Per-sample sheet
# --------------------------------------------------------------------------------------


def test_per_sample_header_row_matches_the_analyst_texts(wb):
    ws = wb["S1"]
    row = list(ws.iter_rows(min_row=SAMPLE_HEADER_ROW, max_row=SAMPLE_HEADER_ROW, values_only=True))[0]
    assert row[COL_LABEL - 1] is None  # blank, as in the analyst sheet
    assert row[COL_CAL_RT - 1] == "Cal RT"
    assert row[COL_MASS - 1] == "Mass (ng)"
    assert row[COL_CORRECTED - 1] == "Corrected Mass"
    assert row[COL_AREA_DUP - 1] == "Area"


def test_per_sample_rows_are_in_canonical_order_at_the_analyst_addresses(wb, mini_result: RunResult):
    """Row 16 + i is analyte i of ``result.analytes``, exactly as in the real workbook
    (verified there: row 16 = congener '1', row 20 = congener '4')."""
    ws = wb["S1"]
    for i, analyte in enumerate(ANALYTES):
        row = SAMPLE_DATA_START_ROW + i
        assert ws.cell(row, COL_LABEL).value == analyte.label
        assert ws.cell(row, COL_CAL_RT).value == pytest.approx(analyte.cal_rt)


def test_area_and_its_duplicate_column_equal_result_areas(wb, mini_result: RunResult):
    ws = wb["S1"]
    for i, analyte in enumerate(ANALYTES):
        row = SAMPLE_DATA_START_ROW + i
        expected = float(mini_result.areas.at[analyte.label, "S1"])
        assert ws.cell(row, 7).value == pytest.approx(expected)  # G = Area
        assert ws.cell(row, COL_AREA_DUP).value == pytest.approx(expected)  # P = Area (dup)


def test_mass_and_corrected_columns_carry_the_uncensored_masses(wb, mini_result: RunResult):
    """Column K mirrors the analyst's own raw ('Corrected Mass') presentation. The analyst
    workbook's per-sample sheet and its SS_Corrected_Mass pivot were verified, by direct
    inspection of the analyst's own workbook, to carry the *uncensored* corrected mass, and
    this column carries the same uncensored value. Censoring is reported in a separate,
    clearly additive column instead."""
    ws = wb["LB1"]
    row = SAMPLE_DATA_START_ROW  # analyte 0 = "1", the deliberately-censored cell
    assert mini_result.masses_censored.at["1", "LB1"] == 0.0
    assert mini_result.masses_corrected.at["1", "LB1"] != 0.0
    assert ws.cell(row, COL_MASS).value == pytest.approx(mini_result.masses.at["1", "LB1"])
    assert ws.cell(row, COL_CORRECTED).value == pytest.approx(
        mini_result.masses_corrected.at["1", "LB1"]
    )


def test_censored_mass_is_a_separate_additive_column(wb, mini_result: RunResult):
    ws = wb["LB1"]
    header = list(ws.iter_rows(min_row=SAMPLE_HEADER_ROW, max_row=SAMPLE_HEADER_ROW, values_only=True))[0]
    censored_col = next(
        i + 1 for i, v in enumerate(header) if isinstance(v, str) and "Censored" in v
    )
    row = SAMPLE_DATA_START_ROW
    assert ws.cell(row, censored_col).value == pytest.approx(0.0)


def test_missing_peak_leaves_measured_columns_blank_but_keeps_cal_rt(wb, mini_result: RunResult):
    """S2 has no peak for '153' (index 3, row 19). RT/Area/Peak are blank, but Cal RT, a
    property of the analyte rather than of this file, is still populated. Mass/Corrected are
    blank too, since quantify_masses never produced a number for it."""
    ws = wb["S2"]
    row = SAMPLE_DATA_START_ROW + LABELS.index("153")
    assert ws.cell(row, COL_LABEL).value == "153"
    assert ws.cell(row, COL_CAL_RT).value == pytest.approx(45.0)
    assert ws.cell(row, 6).value is None  # F = RT
    assert ws.cell(row, 7).value is None  # G = Area
    assert ws.cell(row, COL_PEAK).value is None
    assert ws.cell(row, COL_MASS).value is None


def test_is_column_carries_the_normalizing_standard_label(wb):
    """The analyst fills 'IS' only on the two internal-standard rows themselves, with a raw
    factor. This report instead puts the *normalizing IS label* on every native/surrogate
    row, which is more useful in an automated report. report.truth never reads this column,
    so there is no ground truth to contradict."""
    ws = wb["S1"]
    by_label = {a.label: SAMPLE_DATA_START_ROW + i for i, a in enumerate(ANALYTES)}
    assert ws.cell(by_label["1"], COL_IS).value == "30D"
    assert ws.cell(by_label["153"], COL_IS).value == "204"
    assert ws.cell(by_label["30D"], COL_IS).value is None
    assert ws.cell(by_label["204"], COL_IS).value is None


def test_rrf_column_matches_the_calibration_model(wb, mini_result: RunResult):
    ws = wb["S1"]
    row = SAMPLE_DATA_START_ROW  # "1"
    assert ws.cell(row, COL_RRF).value == pytest.approx(mini_result.calibration.rrf["1"])


def test_freeze_panes_is_below_the_header_row(wb):
    ws = wb["S1"]
    assert ws.freeze_panes == f"A{SAMPLE_DATA_START_ROW}"


# --------------------------------------------------------------------------------------
# cal sheet
# --------------------------------------------------------------------------------------


def test_cal_sheet_header_and_row_addressing(wb, mini_result: RunResult):
    ws = wb["cal"]
    row = list(ws.iter_rows(min_row=CAL_HEADER_ROW, max_row=CAL_HEADER_ROW, values_only=True))[0]
    assert row[CAL_COL_AREA - 1] == "Area"
    assert row[CAL_COL_RRF - 1] == "RRF"
    for i, analyte in enumerate(ANALYTES):
        r = CAL_DATA_START_ROW + i
        assert ws.cell(r, CAL_COL_LABEL).value == analyte.label


def test_cal_sheet_area_and_rrf_match_the_calibration_model(wb, mini_result: RunResult):
    ws = wb["cal"]
    r = CAL_DATA_START_ROW  # "1"
    assert ws.cell(r, CAL_COL_AREA).value == pytest.approx(mini_result.calibration.cal_area["1"])
    assert ws.cell(r, CAL_COL_RRF).value == pytest.approx(mini_result.calibration.rrf["1"])


def test_cal_sheet_invalid_reason_is_reported_when_present(wb, mini_result: RunResult):
    """An analyte the calibration step refused an RRF for gets no RRF cell, but its reason is
    visible in an additive column of its own, since the analyst sheet has no such concept."""
    from quant209.report.excel import write_report as _wr  # local import: avoid unused at top

    result = mini_result
    broken = CalibrationModel(
        rrf={}, cal_area={"1": 10.0}, is_response={}, invalid={"1": "cal peak below floor"},
        injections=("Cal01",), area_floor=500.0,
    )
    import dataclasses

    broken_result = dataclasses.replace(result, calibration=broken)
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        path = _wr(broken_result, Path(tmp) / "broken.xlsx")
        ws = openpyxl.load_workbook(path)["cal"]
        cell = ws.cell(CAL_DATA_START_ROW, CAL_COL_RRF)
        assert cell.value is None
        assert "below floor" in ws.cell(CAL_DATA_START_ROW, CAL_COL_INVALID).value


# --------------------------------------------------------------------------------------
# ref sheet
# --------------------------------------------------------------------------------------


def test_ref_sheet_surrogate_masses_by_homolog_position(wb, mini_result: RunResult):
    ws = wb["ref"]
    header = list(ws.iter_rows(min_row=2, max_row=2, values_only=True))[0]
    ref_col = header.index("Ref01") + 1
    homologs = list(HOMOLOG_SURROGATE)
    mono_row = 3 + homologs.index("mono")
    hexa_row = 3 + homologs.index("hexa")
    expected_mono = float(mini_result.masses.at["13C 3", "Ref01"])
    expected_hexa = float(mini_result.masses.at["13C 153", "Ref01"])
    assert ws.cell(mono_row, ref_col).value == pytest.approx(expected_mono)
    assert ws.cell(hexa_row, ref_col).value == pytest.approx(expected_hexa)


# --------------------------------------------------------------------------------------
# Summary sheet: the Recoveries block
# --------------------------------------------------------------------------------------


def test_summary_recoveries_block_header_and_homolog_rows(wb, mini_result: RunResult):
    ws = wb["Summary"]
    header = list(
        ws.iter_rows(min_row=SUMMARY_HEADER_ROW, max_row=SUMMARY_HEADER_ROW, values_only=True)
    )[0]
    assert header[0] == "Recoveries"
    s1_col = header.index("S1") + 1
    homologs = list(HOMOLOG_SURROGATE)
    mono_row = SUMMARY_DATA_START_ROW + homologs.index("mono")
    hexa_row = SUMMARY_DATA_START_ROW + homologs.index("hexa")
    assert ws.cell(mono_row, 1).value == "mono"
    assert ws.cell(mono_row, s1_col).value == pytest.approx(mini_result.recoveries.at["mono", "S1"])
    assert ws.cell(hexa_row, s1_col).value == pytest.approx(mini_result.recoveries.at["hexa", "S1"])


def test_summary_recoveries_block_covers_every_reported_file(wb, mini_result: RunResult):
    ws = wb["Summary"]
    header = list(
        ws.iter_rows(min_row=SUMMARY_HEADER_ROW, max_row=SUMMARY_HEADER_ROW, values_only=True)
    )[0]
    assert set(mini_result.reported_ids) <= set(header)


# --------------------------------------------------------------------------------------
# run_info: additive audit-trail sheet
# --------------------------------------------------------------------------------------


def test_run_info_sheet_carries_the_audit_trail(wb, mini_result: RunResult):
    import quant209

    ws = wb["run_info"]
    text = "\n".join(
        str(v) for row in ws.iter_rows(values_only=True) for v in row if v is not None
    )
    assert "mini-test-run" in text
    assert quant209.__version__ in text
    assert mini_result.reference_policy in text
    assert "bracketed" in text.lower() or mini_result.reference_policy in text
    # IntegrationParams field values are visible somewhere.
    assert str(mini_result.params.rt_window) in text
    # Including the two label lists, written out as text: a cell holds one scalar, and a
    # report that cannot say which congeners were integrated by rule is not an audit trail.
    assert ", ".join(mini_result.params.forced_riders) in text
    assert ", ".join(mini_result.params.forced_present) in text
    # The skipped analyte for S2 is disclosed.
    assert "153" in text
    assert "S2" in text


def test_run_info_documents_sheets_deliberately_not_reproduced(wb):
    ws = wb["run_info"]
    text = "\n".join(
        str(v) for row in ws.iter_rows(values_only=True) for v in row if v is not None
    )
    assert "cal_info" in text or "Profile" in text or "QuickView" in text


# --------------------------------------------------------------------------------------
# Round-trip and number formats
# --------------------------------------------------------------------------------------


def test_workbook_opens_cleanly_after_round_trip(written_path: Path):
    reopened = openpyxl.load_workbook(written_path)
    assert "Summary" in reopened.sheetnames
    reopened.close()


def test_numeric_cells_carry_numeric_not_text_number_formats(wb):
    ws = wb["S1"]
    row = SAMPLE_DATA_START_ROW
    for col in (7, COL_AREA_DUP, COL_MASS, COL_CORRECTED, COL_RRF):
        cell = ws.cell(row, col)
        assert cell.value is not None
        assert cell.data_type == "n"
        assert cell.number_format != "@"

    rec_ws = wb["Summary"]
    rec_cell = rec_ws.cell(SUMMARY_DATA_START_ROW, 2)
    assert rec_cell.data_type == "n"
    assert "%" in rec_cell.number_format


# --------------------------------------------------------------------------------------
# Real analyst workbook: positional spot-checks (module-scoped real pipeline run)
# --------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def real_result():
    from quant209.core.pipeline import run_analysis
    from quant209.core.sequence import find_data_files, scan_inputs

    folder = run_dir()
    return run_analysis(scan_inputs(find_data_files(folder)), folder / "Converted")


@pytest.fixture(scope="module")
def real_report_path(tmp_path_factory, real_result) -> Path:
    out_dir = tmp_path_factory.mktemp("excel_report")
    return write_report(real_result, out_dir / "report.xlsx", run_name=NAMES.run)


def test_spot_cell_area_matches_the_real_analyst_workbook_position(real_result, real_report_path):
    """Cell P20 on this report's sample sheet is an Area for the same analyte label as the
    analyst's P20 (congener '4'), and its value equals ``result.areas``."""
    analyst = openpyxl.load_workbook(workbook(), data_only=True, read_only=True)
    try:
        analyst_ws = analyst[NAMES.sample]
        analyst_row = list(
            analyst_ws.iter_rows(min_row=20, max_row=20, max_col=16, values_only=True)
        )[0]
    finally:
        analyst.close()
    assert str(analyst_row[2]) == "4"  # column C = label

    ours = openpyxl.load_workbook(real_report_path)
    try:
        ours_ws = ours[NAMES.sample]
        assert str(ours_ws.cell(20, COL_LABEL).value) == "4"
        assert ours_ws.cell(20, COL_AREA_DUP).value == pytest.approx(
            float(real_result.areas.at["4", NAMES.sample])
        )
    finally:
        ours.close()


def test_spot_cell_cal_rrf_matches_the_real_analyst_workbook_position(real_result, real_report_path):
    """The cal sheet's row 6 is congener '1' in both the analyst's workbook and this report's
    (row header + 1, per report.truth's own documented cal-sheet layout)."""
    analyst = openpyxl.load_workbook(workbook(), data_only=True, read_only=True)
    try:
        analyst_row = list(
            analyst["cal"].iter_rows(min_row=6, max_row=6, max_col=8, values_only=True)
        )[0]
    finally:
        analyst.close()
    assert str(analyst_row[2]) == "1"

    ours = openpyxl.load_workbook(real_report_path)
    try:
        ws = ours["cal"]
        assert str(ws.cell(CAL_DATA_START_ROW, CAL_COL_LABEL).value) == "1"
        assert ws.cell(CAL_DATA_START_ROW, CAL_COL_AREA).value == pytest.approx(
            real_result.calibration.cal_area["1"]
        )
        assert ws.cell(CAL_DATA_START_ROW, CAL_COL_RRF).value == pytest.approx(
            real_result.calibration.rrf["1"]
        )
    finally:
        ours.close()


def test_spot_cell_summary_recovery_matches_the_real_analyst_workbook_row_addressing(
    real_result, real_report_path
):
    """Row 3 is the Recoveries header and row 4 is 'mono' in both workbooks, the row
    addressing report.truth itself relies on, even though the column for a given sample
    differs. This run has no single 'REF' column, since it brackets three references."""
    analyst = openpyxl.load_workbook(workbook(), data_only=True, read_only=True)
    try:
        rows = list(analyst["Summary"].iter_rows(min_row=1, max_row=4, max_col=12, values_only=True))
    finally:
        analyst.close()
    header = rows[2]
    assert header[0] == "Recoveries"
    analyst_col = header.index(NAMES.sample)
    analyst_mono = float(rows[3][analyst_col])

    ours = openpyxl.load_workbook(real_report_path)
    try:
        ws = ours["Summary"]
        our_header = list(
            ws.iter_rows(min_row=SUMMARY_HEADER_ROW, max_row=SUMMARY_HEADER_ROW, values_only=True)
        )[0]
        assert our_header[0] == "Recoveries"
        our_col = our_header.index(NAMES.sample) + 1
        our_mono = ws.cell(SUMMARY_DATA_START_ROW, our_col).value
        assert our_mono == pytest.approx(float(real_result.recoveries.at["mono", NAMES.sample]))
        assert our_mono == pytest.approx(analyst_mono, rel=0.35)  # same ballpark, not identical run
    finally:
        ours.close()


def test_e2e_report_generation_from_a_real_run(real_result, tmp_path):
    """End to end: write a real run's report out and confirm it is a sane xlsx file."""
    out_path = write_report(real_result, tmp_path / "report.xlsx", run_name=NAMES.run)
    assert out_path.exists()
    size = out_path.stat().st_size
    assert size > 10_000  # a real multi-sheet workbook, not an empty shell
    reopened = openpyxl.load_workbook(out_path)
    try:
        assert set(real_result.reported_ids) <= set(reopened.sheetnames)
        assert "run_info" in reopened.sheetnames
    finally:
        reopened.close()
