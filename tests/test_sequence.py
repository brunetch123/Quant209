"""Tests for quant209.core.sequence: filename type detection and run-folder discovery."""

from __future__ import annotations

from pathlib import Path

import pytest

from quant209.core.sequence import FileEntry, detect_type, find_data_files


@pytest.mark.parametrize(
    "stem,expected",
    [
        ("AN_Cal01_20231002", "cal"),
        ("AN_Ref04_20231002", "ref"),
        ("AN_Hex03_20231002", "hex"),
        ("AN_LB72_20231002", "blank"),  # lab blank
        ("AN_SMPB3_20231002", "blank"),  # batch/field blank
        ("AN_SMP16_20231002", "sample"),
        ("AB01_Cal-1_something", "cal"),
        # Regression: a site code ending in 'B' must not be misread as
        # a batch blank just because it contains "b" + a digit somewhere.
        ("AN2_ARB01_20231218", "sample"),
        ("AN2_ARB11_20231218", "sample"),
        ("AN2_ARBB1_20231218", "blank"),  # genuine batch blank for the ARB run
    ],
)
def test_detect_type(stem, expected):
    assert detect_type(stem) == expected


def test_file_entry_sample_id_strips_analyst_prefix_and_trailing_date():
    entry = FileEntry(path=Path("AN_SMP16_20231002.mzML"), file_type="sample")
    assert entry.sample_id == "SMP16"


def test_file_entry_name_and_defaults():
    entry = FileEntry(path=Path("AN_SMP16_20231002.mzML"), file_type="sample")
    assert entry.name == "AN_SMP16_20231002.mzML"
    assert entry.auto is True


def test_find_data_files_prefers_converted_mzml_over_d_dirs(reference_run: Path):
    files = find_data_files(reference_run)

    assert files, "expected a non-empty set of Converted/*.mzML files"
    assert all(f.suffix == ".mzML" for f in files)
    assert all(f.parent.name == "Converted" for f in files)


def test_find_data_files_falls_back_to_d_dirs_when_mzml_not_preferred(reference_run: Path):
    files = find_data_files(reference_run, prefer_mzml=False)

    assert files, "expected a non-empty set of .D acquisition directories"
    assert all(f.suffix == ".D" for f in files)


@pytest.mark.parametrize(
    "stem,expected",
    [
        ("AN_SMP16_20231002", "SMP16"),   # the common case: prefix and date on the file
        ("AN2_PQ201", "PQ201"),          # folder carries the date. Prefix still comes off
        ("AN2_LB139", "LB139"),
        ("PQ10_20231008", "PQ10"),        # two parts, but the second is a date, not a sample
        ("09_AN2_PQ201_20250911", "PQ201"),  # acquisition-order prefix as well
        ("ARB01_11", "ARB01_11"),        # not an analyst prefix: both halves kept
        ("AB01_QR05_20220901", "QR05"),   # numeric analyst code, still stripped by the date rule
        ("SomethingLong_X", "SomethingLong_X"),  # not a short analyst code: left alone
    ],
)
def test_sample_id_strips_the_analyst_prefix_not_the_sample(stem, expected):
    """The id has to match the name the analyst's workbook gives the sample, or the file
    cannot be scored against its own ground truth."""
    assert FileEntry(path=Path(f"{stem}.mzML"), file_type="sample").sample_id == expected
