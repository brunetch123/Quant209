"""Load one pooled comparison, as written by ``tools/master_bench.py``, for the figures to share."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path


def _number(text: str) -> float:
    return float(text) if text else 0.0


@dataclass
class Dataset:
    """Every table the figures read, parsed once.

    ``areas`` rows carry ``manual_f`` and ``auto_f``. ``masses`` rows carry ``_m``, ``_a`` and
    ``_det``. Both keep the CSV's own columns as strings.
    """

    name: str
    areas: list[dict]
    masses: list[dict]
    recoveries: list[dict]
    runs: list[dict]
    #: Aroclor names in the reference file's column order, and each one's congener profile.
    aroclor_names: list[str] = field(default_factory=list)
    aroclors: dict[str, dict[str, float]] = field(default_factory=dict)
    #: Per run, how many files of each type the run folder holds. Empty when not scanned.
    file_types: dict[str, dict[str, int]] = field(default_factory=dict)
    #: Draw the legends, axis titles and notes added since the figures were first drawn.
    extras: bool = True
    _memo: dict = field(default_factory=dict, repr=False)

    def memo(self, key: str, build):
        """Return ``build()``, computed once for this dataset and remembered under ``key``."""
        if key not in self._memo:
            self._memo[key] = build()
        return self._memo[key]

    @property
    def sample_natives(self) -> list[dict]:
        """Native analytes in samples: the rows every per-sample figure is built from."""
        return self.memo("sample_natives", lambda: [
            r for r in self.masses if r["role"] == "native" and r["file_type"] == "sample"])


def _read(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def read_aroclors(path: Path) -> tuple[list[str], dict[str, dict[str, float]]]:
    """Read the congener-level Aroclor reference: one row per congener, one column per Aroclor."""
    with open(path, encoding="utf-8") as fh:
        reader = csv.reader(fh)
        names = next(reader)[1:]
        profiles: dict[str, dict[str, float]] = {name: {} for name in names}
        for row in reader:
            for name, value in zip(names, row[1:]):
                profiles[name][row[0]] = float(value)
    return names, profiles


def scan_file_types(root: Path, runs: list[str]) -> dict[str, dict[str, int]]:
    """Count the files of each type in every run folder that is present under ``root``."""
    from quant209.core.sequence import find_data_files, scan_inputs

    counts: dict[str, dict[str, int]] = {}
    for run in runs:
        folder = root / run
        if not folder.is_dir():
            continue
        tally = {"cal": 0, "ref": 0, "blank": 0, "hex": 0, "sample": 0}
        for entry in scan_inputs(find_data_files(folder)):
            tally[entry.file_type] = tally.get(entry.file_type, 0) + 1
        counts[run] = tally
    return counts


def load(data_dir: Path, aroclor: Path | None = None, root: Path | None = None) -> Dataset:
    """Load the comparison in ``data_dir``.

    ``aroclor`` is the Aroclor reference table; without it the fingerprint figures are skipped.
    ``root`` is the archive of run folders; without it the run table leaves its file counts blank.
    """
    areas = _read(data_dir / "areas.csv")
    for r in areas:
        r["manual_f"] = _number(r["manual"])
        r["auto_f"] = _number(r["auto"])
    masses = _read(data_dir / "masses.csv")
    for r in masses:
        r["_m"] = _number(r["manual"])
        r["_a"] = _number(r["auto"])
        r["_det"] = r["detected"] == "True"
    ds = Dataset(
        name=data_dir.name,
        areas=areas,
        masses=masses,
        recoveries=_read(data_dir / "recoveries.csv"),
        runs=_read(data_dir / "runs.csv"),
    )
    if aroclor is not None and aroclor.exists():
        ds.aroclor_names, ds.aroclors = read_aroclors(aroclor)
    if root is not None:
        ds.file_types = scan_file_types(root, [r["run"] for r in ds.runs if r["status"] == "ok"])
    return ds
