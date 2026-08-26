"""Detect a data file's role (calibration, reference, blank, sample, ...) from its
Agilent-instrument filename, and locate the .mzML/.D files that make up a run folder.

Filenames follow an analyst-prefixed convention observed across every run in the archive
(analyst initials or instrument id, then a type-and-number token, then a run date). Written
with placeholder site codes, the shapes are:

    AN_Cal01_...    AB01_Cal-1_...  -> calibration standard
    AN_Ref04_...                    -> reference standard
    AN_Hex01_...                    -> hexane blank (carryover check, not quantified)
    AN_LB72_...                     -> lab blank (quantified, used for LOQ)
    AN_SMPB3_...                    -> field/batch blank (quantified, used for LOQ)
    AN_SMP16_...    AB01_SMP08_...  -> sample
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

FileType = str  # one of: cal, ref, hex, blank, sample

# Ordered (type, regex) rules, most specific first, tested against the filename stem.
_DEFAULT_RULES: list[tuple[str, re.Pattern]] = [
    ("cal", re.compile(r"cal[-_ ]?\d", re.IGNORECASE)),
    ("ref", re.compile(r"ref[-_ ]?\d", re.IGNORECASE)),
    ("hex", re.compile(r"hex[-_ ]?\d", re.IGNORECASE)),
    ("blank", re.compile(r"(^|_)lb\d", re.IGNORECASE)),  # lab blank: LB##
    # batch/field blank: a site code, then 'b', then a single digit. The trailing (?!\d)
    # requires *exactly one* digit after the 'b', so a site code that itself ends in B
    # followed by a two-digit sample number is not misread as a blank. A genuine
    # single-digit sample code ending in 'b' would still false-positive here. The GUI lets
    # the analyst override the detected type for that case.
    ("blank", re.compile(r"[A-Za-z]{2,}b\d(?!\d)", re.IGNORECASE)),
]


@dataclass
class FileEntry:
    """A data file paired with its detected (or user-overridden) type."""

    path: Path
    file_type: FileType
    auto: bool = True  # True if file_type came from detect_type(). False once user-overridden

    @property
    def name(self) -> str:
        return self.path.name

    @property
    def sample_id(self) -> str:
        """The name the analyst's workbook gives this sample, recovered from the filename.

        Acquisition filenames in this archive carry up to three things that are not the
        sample: an acquisition-order number, the analyst's initials, and the run date. For
        example, ``09_AN_PQ201_20250911`` is sample ``PQ201``. Each is stripped only when it
        is unambiguously itself, and the sample id is whatever remains:

        * a leading run of digits is the position in the sequence, never a sample.
        * a trailing 6-8 digit group is the run date (``AN_SMP16_20231002`` -> ``SMP16``).
        * a leading short token of letters (optionally with a digit or two, so that instrument
          codes like ``AB01`` are covered) is the analyst, but only when what follows still
          starts with a letter, so ``ARB01_11`` keeps both halves rather than being read as
          sample ``11``.

        This matters beyond display: the id is what joins a file to its column in the
        analyst's workbook, so a prefix left on is a sample with no ground truth at all.
        """
        parts = self.path.stem.split("_")
        if len(parts) > 1 and re.fullmatch(r"\d{1,3}", parts[0]):
            parts = parts[1:]
        if len(parts) > 1 and re.fullmatch(r"\d{6,8}", parts[-1]):
            parts = parts[:-1]
        if (
            len(parts) > 1
            and re.fullmatch(r"[A-Za-z]{2,4}\d{0,2}", parts[0])
            and parts[1][:1].isalpha()
        ):
            parts = parts[1:]
        return "_".join(parts)


def detect_type(path: str | Path, rules: list[tuple[str, re.Pattern]] | None = None) -> FileType:
    """Classify a file by its name against the ordered rules. Unmatched files are samples."""
    stem = Path(path).stem
    for file_type, pattern in rules or _DEFAULT_RULES:
        if pattern.search(stem):
            return file_type
    return "sample"


def scan_inputs(
    paths: list[str | Path],
    rules: list[tuple[str, re.Pattern]] | None = None,
) -> list[FileEntry]:
    """Build a typed FileEntry for each of a list of .D / .mzML paths."""
    return [FileEntry(path=Path(p), file_type=detect_type(p, rules)) for p in paths]


def find_data_files(folder: str | Path, prefer_mzml: bool = True) -> list[Path]:
    """Find this run's data files: msconvert output (Converted/ or the run folder itself)
    when present and preferred, else the raw Agilent .D acquisition directories."""
    folder = Path(folder)
    mzml = sorted((folder / "Converted").glob("*.mzML")) + sorted(folder.glob("*.mzML"))
    if prefer_mzml and mzml:
        return mzml
    d_files = sorted(folder.glob("*.D"))
    return d_files or mzml
