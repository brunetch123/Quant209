"""Tests for quant209.report.truth, the analyst workbook ground-truth reader.

Every downstream benchmark (area, mass and recovery agreement) joins pipeline output against
the DataFrames read here, so these tests check the reader directly against real workbook
content rather than against a fixture of hand-picked numbers. The run and the names within
it come from ``tests/dataset.py``. Without a configured archive the whole module skips.
"""

from __future__ import annotations

from pathlib import Path

import openpyxl
import pytest

from quant209.core.config import load_canonical_peaks
from quant209.core.standards import HOMOLOG_SURROGATE, internal_standard_label
from quant209.report.truth import RunTruth, read_run_truth

from dataset import NAMES, run_dir


@pytest.fixture(scope="module")
def truth() -> RunTruth:
    """Module-scoped: read_run_truth is a pure function over the on-disk workbooks, so
    every test in this module shares one parse instead of re-reading four ~600KB
    workbooks per test. Computes the run folder path directly (rather than depending on
    conftest.py's function-scoped ``reference_run`` fixture, which a module-scoped fixture
    can't request) so the skip-when-absent behaviour still applies per module import."""
    return read_run_truth(run_dir())


def test_read_run_truth_returns_a_run_truth(truth: RunTruth):
    assert isinstance(truth, RunTruth)
    assert not truth.areas.empty
    assert not truth.masses.empty
    assert not truth.recoveries.empty


def test_at_least_14_quantified_samples(truth: RunTruth):
    """'Quantified' means the analyst actually calculated a mass for at least one analyte
    in that sample column. This excludes the workbooks' blank template sheets (named
    '9'/'10'), whose Mass column is entirely '#DIV/0!' (no injection area behind them)."""
    quantified = truth.masses.columns[truth.masses.notna().any(axis=0)]
    assert len(quantified) >= 14


def test_analyte_labels_align_with_canonical_table(truth: RunTruth):
    """Spot-check the '44+47+65' co-elution group label (see core/config.py's canonical
    table) is present and identically spelled in both the areas and masses index."""
    assert "44+47+65" in truth.areas.index
    assert "44+47+65" in truth.masses.index


def test_areas_and_masses_share_shape_and_orientation(truth: RunTruth):
    """Per the documented convention: index is analyte label, columns are sample id."""
    assert truth.areas.shape == truth.masses.shape
    assert list(truth.areas.index) == list(truth.masses.index)
    assert list(truth.areas.columns) == list(truth.masses.columns)
    assert NAMES.sample in truth.areas.columns


def test_recoveries_indexed_by_homolog_with_fractional_values(truth: RunTruth):
    assert set(truth.recoveries.index) == set(HOMOLOG_SURROGATE)
    values = truth.recoveries.stack()
    assert len(values) > 0
    assert values.between(0.0, 1.5).all()


def test_manual_area_matches_workbook_cell_exactly(reference_run: Path, truth: RunTruth):
    """Derive the expected value by reading the workbook cell directly with openpyxl
    (not a hardcoded magic constant), so this test is self-verifying against the real
    analyst-entered data rather than against the reader's own assumptions."""
    wb = openpyxl.load_workbook(
        reference_run / f"{NAMES.workbook}.xlsx", data_only=True, read_only=True
    )
    try:
        ws = wb[NAMES.sample]
        expected_area = None
        for row in ws.iter_rows(min_row=16, max_col=16, values_only=True):
            if row[2] == "44+47+65":  # column C (index 2) = analyte label
                expected_area = float(row[15])  # column P (index 15) = area
                break
    finally:
        wb.close()

    assert expected_area is not None, "fixture workbook layout changed, cell not found"
    assert truth.areas.loc["44+47+65", NAMES.sample] == expected_area


def test_recoveries_match_workbook_summary_cells_by_row_position(
    reference_run: Path, truth: RunTruth
):
    """Read the Summary sheet's Recoveries block directly. Mono is row 4 (header row 3
    + 1) and deca is row 13 (header row 3 + 10), checking that the by-position homolog
    mapping landed on the right rows.

    A test that only range-checks recoveries (0-1.5) can't catch a row-position
    off-by-one: every homolog's recovery is a similarly-shaped fraction, so a shift
    would silently pair "mono" with the tri-CB row's numbers and still pass a pure
    range check. Pinning two specific cells (first and last row of the block) catches
    that class of bug.
    """
    wb = openpyxl.load_workbook(
        reference_run / f"{NAMES.workbook}.xlsx", data_only=True, read_only=True
    )
    try:
        ws = wb["Summary"]
        rows = list(ws.iter_rows(min_row=1, max_row=13, max_col=12, values_only=True))
    finally:
        wb.close()

    header = rows[2]  # row 3 (1-indexed) -> index 2
    assert header[0] == "Recoveries", "fixture workbook layout changed, header not found"
    col = header.index(NAMES.sample)

    mono_row = rows[3]  # row 4  = header_row(3) + 1 -> mono, per HOMOLOG_SURROGATE order
    deca_row = rows[12]  # row 13 = header_row(3) + 10 -> deca

    assert truth.recoveries.loc["mono", NAMES.sample] == float(mono_row[col])
    assert truth.recoveries.loc["deca", NAMES.sample] == float(deca_row[col])


def test_merges_samples_across_every_workbook_in_the_run_folder(truth: RunTruth):
    """The reference run holds four analyst workbooks, two of them covering the same
    samples. A sample from each distinct one should land in the merged result."""
    assert NAMES.sample in truth.areas.columns
    assert NAMES.sample_in_workbook_2 in truth.areas.columns
    assert NAMES.sample_in_workbook_3 in truth.areas.columns  # the re-issued workbook


def test_sequence_workbook_is_not_treated_as_a_quant_workbook(truth: RunTruth):
    """The run's ``*_sequence.xlsx`` is the instrument sequence log, not analyst
    quantification data, and must not contribute spurious columns/rows."""
    assert "SEQ" not in truth.areas.columns
    assert "Sheet1" not in truth.areas.columns


# --------------------------------------------------------------------------------------
# cal sheet: the calibration injection every RRF in a workbook came from
# --------------------------------------------------------------------------------------


def test_cal_frames_are_analyte_by_workbook(truth: RunTruth):
    """Per the documented convention: index is analyte label, columns are the *workbook*
    that holds the cal sheet, not a sample and not an injection id. See the module
    docstring for why the sheet cannot name its own injection."""
    assert truth.cal_areas.shape == truth.cal_rrfs.shape
    assert list(truth.cal_areas.index) == list(truth.cal_rrfs.index)
    assert list(truth.cal_areas.columns) == list(truth.cal_rrfs.columns)
    assert NAMES.workbook in truth.cal_areas.columns
    # The cal standard contains every congener, so the cal sheet covers the whole analyte
    # table rather than a sample's detectable subset.
    assert set(truth.cal_areas.index) == set(truth.areas.index)


def test_cal_area_and_rrf_match_workbook_cells_exactly(reference_run: Path, truth: RunTruth):
    """Self-verifying against the real sheet: read column E (area) and column H (RRF) of
    the ``cal`` sheet directly rather than asserting hardcoded constants."""
    wb = openpyxl.load_workbook(
        reference_run / f"{NAMES.workbook}.xlsx", data_only=True, read_only=True
    )
    try:
        ws = wb["cal"]
        expected: dict[str, tuple[float, float]] = {}
        for row in ws.iter_rows(min_row=6, max_col=8, values_only=True):
            # Congener labels that are bare numbers come back from openpyxl as ints, not
            # strings. The reader stringifies them, so the comparison has to as well.
            label = None if row[2] is None else str(row[2])
            if label in {"87+125", "1", "13C 153"}:
                expected[label] = (float(row[4]), float(row[7]))  # E = area, H = RRF
    finally:
        wb.close()

    assert set(expected) == {"87+125", "1", "13C 153"}, "fixture cal layout changed"
    for label, (area, rrf) in expected.items():
        assert truth.cal_areas.loc[label, NAMES.workbook] == area
        assert truth.cal_rrfs.loc[label, NAMES.workbook] == rrf


def test_cal_rrf_matches_cal_info_masses(reference_run: Path, truth: RunTruth):
    """The analyst's RRF column is reproducible from their area column and the certified
    calibration masses on the workbook's own ``cal_info`` sheet.

    ``RRF = (cal_mass / cal_area) / response``, where ``response`` is the analyte's
    internal standard's own ``cal_mass / cal_area``. Which standard that is comes from the
    canonical table's ``is_norm`` field rather than from an assumption here. The split is
    not the obvious one (the penta group '87+125' is normalised by d-PCB30, not by PCB204),
    and checking both standards is what proves the area and RRF columns were read from the
    same rows: a one-row offset in either would break the identity for every analyte at
    once rather than for one.
    """
    book = NAMES.workbook
    wb = openpyxl.load_workbook(reference_run / f"{book}.xlsx", data_only=True, read_only=True)
    try:
        cal_mass = {
            str(row[0]): float(row[1])
            for row in wb["cal_info"].iter_rows(min_row=5, max_col=2, values_only=True)
            if row[0] is not None and isinstance(row[1], (int, float))
        }
    finally:
        wb.close()

    analytes = load_canonical_peaks()
    by_label = {a.label: a for a in analytes}

    def response(is_label: str) -> float:
        return cal_mass[is_label] / float(truth.cal_areas.loc[is_label, book])

    checked: set[str] = set()
    for label in ("1", "12+13", "87+125", "206", "13C 153"):
        is_label = internal_standard_label(by_label[label].is_norm, analytes)
        expected = (cal_mass[label] / float(truth.cal_areas.loc[label, book])) / response(is_label)
        assert truth.cal_rrfs.loc[label, book] == pytest.approx(expected, rel=1e-9)
        checked.add(is_label)
    assert checked == {"30D", "204"}, "both internal standards must be exercised"


def test_cal_columns_are_distinct_injections(truth: RunTruth):
    """Each workbook's cal sheet holds a *different* calibration injection's areas. The
    workbooks are not copies of one another. Guards against a reader that silently wrote
    every column from the first workbook it opened."""
    populated = truth.cal_areas.columns[truth.cal_areas.notna().any(axis=0)]
    assert len(populated) >= 3
    congener_1 = truth.cal_areas.loc["1", populated].dropna()
    assert congener_1.nunique() == len(congener_1)


def test_superseded_workbook_contributes_an_unusable_cal_column(truth: RunTruth):
    """One workbook in the run was superseded by a later re-issue, and its cal sheet was
    never filled in: every area cell is a literal 0 and every RRF cell is '#DIV/0!'.

    It must still occupy a column. "The analyst produced no numbers here" is a fact a
    benchmark should be able to see and screen out, not one the reader should hide. The
    areas are reported as the zeros the sheet actually holds. Only the Excel error strings
    become NaN.
    """
    stale = NAMES.workbook_superseded
    assert stale in truth.cal_areas.columns
    assert (truth.cal_areas[stale] == 0.0).all()
    assert truth.cal_rrfs[stale].isna().all()
    assert truth.ref_surrogate_masses[stale].isna().all()


# --------------------------------------------------------------------------------------
# ref sheet: the reference injection every recovery in a workbook was measured against
# --------------------------------------------------------------------------------------


def test_ref_surrogate_masses_are_homolog_by_workbook_in_mono_to_deca_order(truth: RunTruth):
    assert list(truth.ref_surrogate_masses.index) == list(HOMOLOG_SURROGATE)
    assert NAMES.workbook in truth.ref_surrogate_masses.columns


def test_ref_surrogate_mass_matches_workbook_cell_exactly(reference_run: Path, truth: RunTruth):
    """Read the ``ref`` sheet's summary block directly. The first workbook covers the
    batch the pipeline brackets to ``Ref04`` (see ``core.pipeline``'s module docstring for
    that evidence), so this cell is Ref04's 13C 153 surrogate mass: ~11.258 ng."""
    wb = openpyxl.load_workbook(
        reference_run / f"{NAMES.workbook}.xlsx", data_only=True, read_only=True
    )
    try:
        rows = list(wb["ref"].iter_rows(min_row=1, max_row=12, max_col=20, values_only=True))
    finally:
        wb.close()

    header = rows[1]  # row 2 (1-indexed) -> index 1
    assert header[13] == "Recoveries", "fixture ref layout changed, header not found"
    mass_col = header.index("Mass")

    # Mono is the row right under the header. Hexa is five rows further down.
    assert truth.ref_surrogate_masses.loc["mono", NAMES.workbook] == float(
        rows[2][mass_col]
    )
    assert truth.ref_surrogate_masses.loc["hexa", NAMES.workbook] == float(
        rows[7][mass_col]
    )
    assert truth.ref_surrogate_masses.loc["hexa", NAMES.workbook] == pytest.approx(
        11.258, abs=0.001
    )


def test_ref_surrogate_masses_read_by_row_position_not_label(truth: RunTruth):
    """The ref block labels the tri-CB surrogate '13C 31', while the canonical table calls
    it '13C 28'. Homolog identity therefore comes from row position, exactly as in the
    Summary recoveries block, so 'tri' must be populated despite no label ever matching.
    """
    assert HOMOLOG_SURROGATE["tri"] == "13C 28"
    assert truth.ref_surrogate_masses.loc["tri"].notna().any()


def test_ref_surrogate_masses_are_plausible_spike_recoveries(truth: RunTruth):
    """Each surrogate is spiked at 9.84 ng into the reference injection, which skipped the
    extraction the samples went through, so its measured mass sits near that spike. A
    reader that landed on the wrong column (the recovery column, a constant 1.0, sits
    between the label and the mass) would fail this."""
    values = truth.ref_surrogate_masses.stack().dropna()
    assert len(values) >= 30  # three usable workbooks x ten homologs
    assert values.between(3.0, 20.0).all()


def test_workbook_samples_maps_each_workbook_to_its_own_sample_sheets(truth: RunTruth):
    """Needed to resolve a cal/ref column to an injection: a caller intersects these
    samples with the run's own reference bracketing."""
    assert NAMES.sample in truth.workbook_samples[NAMES.workbook]
    assert NAMES.sample_in_workbook_2 in truth.workbook_samples[NAMES.workbook_2]
    assert NAMES.sample_in_workbook_3 in truth.workbook_samples[NAMES.workbook_3]
    # Every real sample column in ``areas`` traces back to exactly one workbook. The
    # numerically-named sheets ('8', '9', '10') are the unused blank slots of the workbook
    # template and repeat across every workbook, so they are excluded from that claim.
    # See ``test_at_least_14_quantified_samples`` for what they are.
    owners = [
        s
        for samples in truth.workbook_samples.values()
        for s in samples
        if not s.isdigit()
    ]
    assert len(owners) == len(set(owners))
    assert set(truth.areas.columns) <= set(owners) | {"8", "9", "10"}
