"""Synthesise a complete, fictional analysis run.

The laboratory's own data is not distributed with this project, which leaves anyone who clones
it with nothing to point the application at. This writes a run folder of valid SRM mzML files
(a calibration injection, a reference injection, a lab blank, a hexane rinse and six samples)
built from the canonical peak table rather than from any acquisition. Every number in it is
invented. It is what the screenshots in the README were taken against, and it is enough to
click through all six steps of the application.

    python tools/make_demo_run.py --out .

Writes ``DEMO01-06_20240115/Converted/*.mzML`` under ``--out``.
"""

from __future__ import annotations

import argparse
import base64
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quant209.core.config import analytes_by_transition, load_canonical_peaks
from quant209.core.standards import HOMOLOG_SURROGATE

RUN_NAME = "DEMO01-06_20240115"

#: Peak area in counts·s per ng of analyte. Chosen so the calibration injection's areas land in
#: the same 10^4-10^6 range a real one does.
AREA_PER_NG = 4000.0

#: Gaussian half-width in minutes. Real PCB peaks on this method are ~0.07 min at half height.
PEAK_SIGMA = 0.03

#: Samples per minute of the synthesised time grid.
POINTS_PER_MINUTE = 200.0

#: Minutes of quiet grid kept either side of a transition's own analytes.
GRID_MARGIN = 2.0

#: ng of each internal standard spiked into every quantified file.
IS_SPIKE_NG = 10.4125

#: ng of each 13C surrogate spiked into every quantified file and the reference injection.
SURROGATE_SPIKE_NG = 9.84

BASELINE = 40.0
NOISE = 12.0


def homolog_profile(rng: np.random.Generator) -> dict[str, float]:
    """A per-homolog concentration multiplier, so each sample looks like its own mixture."""
    return {
        h: float(np.clip(rng.lognormal(mean=-0.4, sigma=0.8), 0.02, 6.0))
        for h in ("mono", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona", "deca")
    }


def native_ng(analyte, profile: dict[str, float], rng: np.random.Generator) -> float:
    """How many ng of one native congener this file contains."""
    if rng.random() < 0.35:  # most congeners sit below the reporting limit in a real sample
        return float(rng.uniform(0.0, 0.25))
    return float(profile.get(analyte.homolog, 1.0) * rng.lognormal(-0.2, 0.7))


def masses_for(kind: str, analytes, seed: int) -> dict[str, float]:
    """ng of every analyte in one file, keyed by analyte label."""
    rng = np.random.default_rng(seed)
    profile = homolog_profile(rng)
    surrogates = set(HOMOLOG_SURROGATE.values())
    masses: dict[str, float] = {}

    for a in analytes:
        if a.role in ("IS_d30", "IS_204"):
            masses[a.label] = 0.0 if kind == "hex" else IS_SPIKE_NG
        elif a.role == "surrogate_13C" or a.label in surrogates:
            if kind == "hex":
                masses[a.label] = 0.0
            elif kind == "cal":
                masses[a.label] = a.cal_mass
            elif kind == "ref":
                masses[a.label] = SURROGATE_SPIKE_NG
            else:
                # Extraction losses: a sample recovers 65-95% of what was spiked into it.
                masses[a.label] = SURROGATE_SPIKE_NG * float(rng.uniform(0.65, 0.95))
        elif kind == "cal":
            masses[a.label] = a.cal_mass
        elif kind == "sample":
            masses[a.label] = native_ng(a, profile, rng)
        elif kind == "blank":
            masses[a.label] = float(rng.uniform(0.0, 0.35))
        elif kind == "ref":
            masses[a.label] = float(rng.uniform(0.0, 0.05))
        else:
            masses[a.label] = 0.0
    return masses


def chromatogram(group, masses: dict[str, float], rng: np.random.Generator):
    """One transition's (time, intensity) arrays: a Gaussian per analyte on a noisy baseline."""
    rts = [a.cal_rt for a in group]
    start, stop = min(rts) - GRID_MARGIN, max(rts) + GRID_MARGIN
    n = int((stop - start) * POINTS_PER_MINUTE) + 1
    time = np.linspace(start, stop, n)

    intensity = np.full(n, BASELINE) + rng.normal(0.0, NOISE, n)
    for a in group:
        area = masses.get(a.label, 0.0) * AREA_PER_NG
        if area <= 0.0:
            continue
        # A Gaussian's area is amplitude * sigma * sqrt(2*pi). Areas are in counts·s while the
        # grid is in minutes, hence the 60.
        amplitude = area / (PEAK_SIGMA * np.sqrt(2.0 * np.pi) * 60.0)
        intensity += amplitude * np.exp(-0.5 * ((time - a.cal_rt) / PEAK_SIGMA) ** 2)
    return time, np.clip(intensity, 0.0, None)


def encode(values: np.ndarray) -> str:
    return base64.b64encode(np.asarray(values, dtype="<f4").tobytes()).decode("ascii")


CHROMATOGRAM_TEMPLATE = """        <chromatogram index="{index}" id="SRM SIC Q1={q1} Q3={q3}" defaultArrayLength="{n}">
          <cvParam cvRef="MS" accession="MS:1001473" name="selected reaction monitoring chromatogram" value=""/>
          <precursor><isolationWindow>
            <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{q1}"/>
          </isolationWindow></precursor>
          <product><isolationWindow>
            <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{q3}"/>
          </isolationWindow></product>
          <binaryDataArrayList count="2">
            <binaryDataArray>
              <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" value=""/>
              <cvParam cvRef="MS" accession="MS:1000595" name="time array" value="" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>
              <binary>{time}</binary>
            </binaryDataArray>
            <binaryDataArray>
              <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" value=""/>
              <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/>
              <binary>{intensity}</binary>
            </binaryDataArray>
          </binaryDataArrayList>
        </chromatogram>"""

MZML_TEMPLATE = """<?xml version="1.0" encoding="utf-8"?>
<indexedmzML xmlns="http://psi.hupo.org/ms/mzml">
  <mzML xmlns="http://psi.hupo.org/ms/mzml" id="{run_name}" version="1.1.0">
    <run id="demo" startTimeStamp="{stamp}">
      <chromatogramList count="{count}">
{chromatograms}
      </chromatogramList>
    </run>
  </mzML>
</indexedmzML>
"""


def build_mzml(chroms, acquired: datetime) -> bytes:
    blocks = [
        CHROMATOGRAM_TEMPLATE.format(
            index=index, q1=q1, q3=q3, n=len(time),
            time=encode(time), intensity=encode(intensity),
        )
        for index, ((q1, q3), time, intensity) in enumerate(chroms)
    ]
    xml = MZML_TEMPLATE.format(
        run_name=RUN_NAME,
        stamp=acquired.strftime("%Y-%m-%dT%H:%M:%SZ"),
        count=len(chroms),
        chromatograms="\n".join(blocks),
    )
    return xml.encode("utf-8")


#: The sequence, in acquisition order: (filename stem, kind).
SEQUENCE = [
    ("AN_Cal01_20240115", "cal"),
    ("AN_Ref01_20240115", "ref"),
    ("AN_LB01_20240115", "blank"),
    ("AN_DEMO01_20240115", "sample"),
    ("AN_DEMO02_20240115", "sample"),
    ("AN_DEMO03_20240115", "sample"),
    ("AN_Hex01_20240115", "hex"),
    ("AN_DEMO04_20240115", "sample"),
    ("AN_DEMO05_20240115", "sample"),
    ("AN_DEMO06_20240115", "sample"),
]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--out", default=".", help="directory to create the run folder in")
    args = ap.parse_args(argv)

    analytes = load_canonical_peaks()
    groups = analytes_by_transition(analytes)
    converted = Path(args.out) / RUN_NAME / "Converted"
    converted.mkdir(parents=True, exist_ok=True)

    start = datetime(2024, 1, 15, 8, 0, tzinfo=timezone.utc)
    for i, (stem, kind) in enumerate(SEQUENCE):
        masses = masses_for(kind, analytes, seed=1000 + i)
        rng = np.random.default_rng(2000 + i)
        chroms = [(key, *chromatogram(group, masses, rng)) for key, group in groups.items()]
        path = converted / f"{stem}.mzML"
        path.write_bytes(build_mzml(chroms, start + timedelta(minutes=52 * i)))
        print(f"  {path.name}  {path.stat().st_size / 1e6:.1f} MB  ({kind})")

    print(f"\n{RUN_NAME} written to {Path(args.out).resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
