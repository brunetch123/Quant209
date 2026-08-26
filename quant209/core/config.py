"""Load the canonical peak/transition table shipped with the package.

Each row describes one reported analyte: a single congener, an inseparable co-elution
group (e.g. ``44+47+65``, chromatographically inseparable congeners reported as one
combined peak), or a labelled standard (a ``13C``-surrogate or internal-standard channel).
Every row ties its analyte to the MRM transition it's found in, its elution order there,
the homolog class, which internal standard normalises it, and its certified calibration-
standard mass. The table is method-/column-specific and re-derivable from any run's
calibration sheet (see ``tools/extract_config.py``).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources

_CSV_PACKAGE = "quant209.config"
_CSV_NAME = "canonical_peaks.csv"


@dataclass(frozen=True)
class Analyte:
    """One reported analyte and how to find/quantify it."""

    label: str  # congener number, co-elution group ("44+47+65"), or standard label ("13C 3", "30D")
    q1: float  # precursor (Q1) m/z of this analyte's MRM transition
    q3: float  # product (Q3) m/z of this analyte's MRM transition
    homolog: str  # mono..deca, or "" for labelled standards that don't belong to a native homolog
    order: int  # 1-based elution order among analytes sharing this transition
    role: str  # native | surrogate_13C | IS_d30 | IS_204
    is_norm: str  # internal standard that normalises this analyte: "30D" or "204"
    cal_mass: float  # certified mass of this analyte in the calibration standard
    cal_rt: float  # analyst calibration retention time (anchor / cross-check)


@lru_cache()
def load_canonical_peaks() -> list[Analyte]:
    """Load the packaged canonical analyte table (185 analytes), in file order.

    Cached: this is called once per file processed downstream, and the packaged CSV
    never changes within a process lifetime. Callers must treat the returned list as
    read-only. Mutating it would corrupt every subsequent caller's view.
    """
    csv_path = resources.files(_CSV_PACKAGE).joinpath(_CSV_NAME)
    analytes: list[Analyte] = []
    with csv_path.open("r", newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            try:
                analytes.append(
                    Analyte(
                        label=row["label"],
                        q1=float(row["transition_q1"]),
                        q3=float(row["transition_q3"]),
                        homolog=row["homolog"],
                        order=int(row["order_in_transition"]),
                        role=row["role"],
                        is_norm=row["is_norm"],
                        cal_mass=float(row["cal_mass"]),
                        cal_rt=float(row["cal_rt"]),
                    )
                )
            except ValueError as exc:
                raise ValueError(
                    f"{_CSV_NAME}: malformed row for label {row['label']!r}: {exc}"
                ) from exc
    return analytes


def analytes_by_transition(analytes: list[Analyte]) -> dict[tuple[float, float], list[Analyte]]:
    """Group analytes by (Q1, Q3) transition, each list sorted by calibration RT.

    Elution order within a transition is what the integration engine walks peak-by-peak,
    so it must reflect actual retention time rather than the CSV's row order. Keys are
    rounded to 1 decimal place, the same convention ``core.mzml.read_srm_mzml`` uses for
    its own (Q1, Q3) keys, to absorb acquisition jitter. This makes every native/
    co-elution transition key join exactly against a real reader's keys (verified in
    ``tests/test_config.py``). Labelled-standard (surrogate/IS) keys are chemically
    *derived* in the source CSV rather than read off an instrument, and can differ from
    a real acquisition's key by ~0.1 on this same rounding grid. Callers that need to
    locate a labelled standard's own chromatogram should fall back to a tolerant nearest-
    transition lookup rather than assuming an exact match.
    """
    grouped: dict[tuple[float, float], list[Analyte]] = {}
    for analyte in analytes:
        key = (round(analyte.q1, 1), round(analyte.q3, 1))
        grouped.setdefault(key, []).append(analyte)
    for key in grouped:
        grouped[key].sort(key=lambda a: a.cal_rt)
    return grouped
