"""Tests for the canonical peak table loader (quant209.core.config) and the
internal-standard/surrogate mapping built on top of it (quant209.core.standards)."""

from __future__ import annotations

from pathlib import Path

import pytest

from quant209.core.config import Analyte, analytes_by_transition, load_canonical_peaks
from quant209.core.mzml import read_srm_mzml
from quant209.core.standards import HOMOLOG_SURROGATE, internal_standard_label

from dataset import NAMES


def test_load_canonical_peaks_yields_185_analytes():
    analytes = load_canonical_peaks()
    assert len(analytes) == 185
    assert all(isinstance(a, Analyte) for a in analytes)


def test_mono_transition_contains_natives_ordered_by_cal_rt():
    grouped = analytes_by_transition(load_canonical_peaks())
    mono = grouped[(188.0, 153.1)]

    assert [a.role for a in mono] == ["native", "native", "native"]
    cal_rts = [a.cal_rt for a in mono]
    assert cal_rts == sorted(cal_rts)


def test_roles_include_native_surrogate_and_internal_standard_rows():
    roles = {a.role for a in load_canonical_peaks()}
    assert "native" in roles
    assert "surrogate_13C" in roles
    assert {"IS_d30", "IS_204"} <= roles


def test_every_analyte_has_a_valid_is_norm():
    for a in load_canonical_peaks():
        assert a.is_norm in ("30D", "204")


def test_every_analyte_has_positive_cal_mass_and_cal_rt():
    for a in load_canonical_peaks():
        assert a.cal_mass > 0
        assert a.cal_rt > 0


def test_homolog_is_blank_only_for_labeled_rows():
    for a in load_canonical_peaks():
        if a.homolog == "":
            assert a.role != "native"


def test_every_analyte_has_a_label():
    for a in load_canonical_peaks():
        assert a.label


def test_coelution_label_exists():
    """44+47+65 is an inseparable co-elution group reported as a single combined peak."""
    labels = {a.label for a in load_canonical_peaks()}
    assert "44+47+65" in labels


@pytest.mark.parametrize("is_norm,expected_label", [("30D", "30D"), ("204", "204")])
def test_internal_standard_label_resolves_from_csv(is_norm, expected_label):
    analytes = load_canonical_peaks()
    assert internal_standard_label(is_norm, analytes) == expected_label


def test_internal_standard_label_reads_the_row_not_the_is_norm_tag():
    """Synthetic rows where the IS label differs from its own is_norm tag: proves the
    lookup actually reads label off the matching row, rather than e.g. returning
    is_norm verbatim (which would pass the CSV-backed test above by coincidence, since
    in the real table the IS rows' labels happen to equal their is_norm tags)."""
    synthetic = [
        Analyte(
            label="d30-IS", q1=261.0, q3=191.0, homolog="", order=1,
            role="IS_d30", is_norm="30D", cal_mass=1.0, cal_rt=1.0,
        ),
        Analyte(
            label="204-IS", q1=429.8, q3=359.8, homolog="octa", order=3,
            role="IS_204", is_norm="204", cal_mass=1.0, cal_rt=1.0,
        ),
    ]

    assert internal_standard_label("30D", synthetic) == "d30-IS"
    assert internal_standard_label("204", synthetic) == "204-IS"


def test_homolog_surrogate_labels_exist_among_surrogate_13c_rows():
    surrogate_labels = {a.label for a in load_canonical_peaks() if a.role == "surrogate_13C"}
    assert set(HOMOLOG_SURROGATE.values()) <= surrogate_labels


def test_native_transition_keys_join_with_real_mzml_reader_keys(reference_run: Path):
    """The rounded (Q1, Q3) key config.py hands the integration engine must match what
    core.mzml.read_srm_mzml actually produces for a real acquisition file, for every
    transition that carries a native/co-elution analyte. These are what the
    integration engine walks peak-by-peak.

    Not asserted here: full equality of config's key set with the reader's.
    canonical_peaks.csv derives each labelled-standard's (surrogate, IS) Q1/Q3 as the
    native homolog's transition plus the isotope's *nominal* mass shift (13C12 -> +12,
    d5 -> +5), not by reading the true acquired isolation-window target. Checked against
    a real acquisition, that nominal-shift approximation lands exactly on 2 of 11
    labelled-standard channels and is 0.1 Da off on the other 9 (e.g. 13C-3's config key
    (200.0, 165.1) against the file's actual (200.1, 165.1)). It is a characteristic of
    the canonical table, which is why labelled-standard channels are looked up by nearest
    transition within a tolerance rather than by exact key.
    """
    converted = reference_run / "Converted"
    candidates = sorted(p for p in converted.glob("*.mzML") if NAMES.sample in p.name)
    assert candidates, f"no sample mzML found under {converted}"

    reader_keys = set(read_srm_mzml(candidates[0]).keys())
    grouped = analytes_by_transition(load_canonical_peaks())
    native_keys = {key for key, analytes in grouped.items() if any(a.role == "native" for a in analytes)}

    assert len(native_keys) == 10  # one MRM transition per homolog, mono..deca
    assert native_keys <= reader_keys
