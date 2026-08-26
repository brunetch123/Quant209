"""Tests for quant209.core.mzml: the SRM mzML reader."""

import base64
from pathlib import Path

import numpy as np
import pytest

from quant209.core.mzml import Chromatogram, read_srm_mzml

from dataset import NAMES

# The acquisition method defines 22 distinct SRM transitions: 10 native homolog
# transitions plus 12 labelled-standard transitions (13C surrogates and internal
# standards). That includes the deuterated PCB65-d5 surrogate at 297.0 -> 227.0, which
# is acquired but is *not* one of the 21 transitions in the downstream quantification
# table (quant209/config/canonical_peaks.csv).
EXPECTED_TRANSITION_COUNT = 22

NATIVE_HOMOLOG_TRANSITIONS = [
    (188.0, 153.1),  # mono
    (222.0, 152.1),  # di
    (256.0, 186.0),  # tri
    (291.9, 222.0),  # tetra
    (325.9, 255.9),  # penta
    (359.8, 289.9),  # hexa
    (393.8, 323.9),  # hepta
    (429.8, 359.8),  # octa
    (463.7, 393.8),  # nona
    (497.7, 427.7),  # deca
]


def _find_sample_mzml(reference_run: Path) -> Path:
    """Locate one fully-converted sample mzML in this run's Converted folder."""
    converted = reference_run / "Converted"
    candidates = sorted(p for p in converted.glob("*.mzML") if NAMES.sample in p.name)
    assert candidates, f"no sample mzML found under {converted}"
    return candidates[0]


def test_reads_21_transitions(reference_run: Path):
    """A real acquisition decodes into 22 distinct (Q1, Q3) transitions.

    (Named for the 21 canonical quantification transitions. The raw acquisition carries
    one extra labelled-standard channel not used downstream. See the comment above.)
    """
    path = _find_sample_mzml(reference_run)
    chroms = read_srm_mzml(path)

    assert len(chroms) == EXPECTED_TRANSITION_COUNT
    for key, chrom in chroms.items():
        assert isinstance(chrom, Chromatogram)
        # Keys must be the rounded-to-1-dp (Q1, Q3) pair, consistent with the stored fields.
        assert key == (round(chrom.q1, 1), round(chrom.q3, 1))
        assert key == (round(key[0], 1), round(key[1], 1))
        assert chrom.times.size == chrom.intensities.size
        assert chrom.times.size > 0
        assert np.all(np.diff(chrom.times) > 0), f"times not strictly increasing for {key}"
        assert np.all(chrom.intensities >= 0), f"negative intensity found for {key}"

    for key in NATIVE_HOMOLOG_TRANSITIONS:
        assert key in chroms, f"missing native transition {key}"
        assert chroms[key].times.size > 1000


def test_truncated_file_returns_empty(tmp_path: Path):
    """A header-only, truncated-mid-document mzML must not raise. It must just yield {}."""
    stub = tmp_path / "truncated.mzML"
    stub.write_text(
        '<?xml version="1.0" encoding="utf-8"?>\n'
        '<indexedmzML xmlns="http://psi.hupo.org/ms/mzml">\n'
        '  <mzML xmlns="http://psi.hupo.org/ms/mzml" id="stub" version="1.1.0">\n'
        '    <run id="stub_run">\n',
        encoding="utf-8",
    )

    assert read_srm_mzml(stub) == {}


def _minimal_srm_mzml(time_b64: str, intensity_b64: str, q1: float, q3: float, dtype_bits: int = 64) -> bytes:
    """Build a small, structurally-complete mzML with exactly one SRM chromatogram.

    Unlike ``test_truncated_file_returns_empty`` (which cuts off the XML document itself),
    this produces a fully well-formed document so callers can control corruption at the
    binary-payload level instead (see ``test_mid_array_truncation_is_tolerated``).
    """
    float_cv_name = "32-bit float" if dtype_bits == 32 else "64-bit float"
    xml = f"""<?xml version="1.0" encoding="utf-8"?>
<indexedmzML xmlns="http://psi.hupo.org/ms/mzml">
  <mzML xmlns="http://psi.hupo.org/ms/mzml" id="stub" version="1.1.0">
    <run id="stub_run">
      <chromatogramList count="1">
        <chromatogram index="0" id="SRM SIC Q1={q1} Q3={q3}" defaultArrayLength="0">
          <cvParam cvRef="MS" accession="MS:1001473" name="selected reaction monitoring chromatogram" value=""/>
          <precursor>
            <isolationWindow>
              <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{q1}"/>
            </isolationWindow>
          </precursor>
          <product>
            <isolationWindow>
              <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{q3}"/>
            </isolationWindow>
          </product>
          <binaryDataArrayList count="2">
            <binaryDataArray>
              <cvParam cvRef="MS" accession="MS:1000523" name="{float_cv_name}" value=""/>
              <cvParam cvRef="MS" accession="MS:1000595" name="time array" value="" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>
              <binary>{time_b64}</binary>
            </binaryDataArray>
            <binaryDataArray>
              <cvParam cvRef="MS" accession="MS:1000523" name="{float_cv_name}" value=""/>
              <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/>
              <binary>{intensity_b64}</binary>
            </binaryDataArray>
          </binaryDataArrayList>
        </chromatogram>
      </chromatogramList>
    </run>
  </mzML>
</indexedmzML>
"""
    return xml.encode("utf-8")


def _encode_floats(values: np.ndarray, dtype_bits: int = 64) -> str:
    raw = np.asarray(values, dtype=f"<f{dtype_bits // 8}").tobytes()
    return base64.b64encode(raw).decode("ascii")


def test_mid_array_truncation_is_tolerated(tmp_path: Path):
    """A structurally-complete mzML whose intensity payload is a few bytes short of a
    whole float64 must not raise.

    The byte buffer is decoded with ``struct.unpack`` after first being trimmed to a whole
    multiple of the element width. A buffer that is valid base64 (and so raises nothing at
    the base64/zlib stage) but not an exact multiple of 8 bytes long would otherwise raise
    ``struct.error``, a type the truncation tolerance would not catch, crashing the whole
    read.

    Behaviour under truncation: the corrupt array is *partially decoded*, its trailing
    partial float dropped. Since the two arrays of this transition then disagree in length,
    both are trimmed to their shared length rather than the transition being dropped
    entirely.
    """
    time_vals = np.linspace(0.0, 1.0, 10)
    time_b64 = _encode_floats(time_vals)

    full_intensity_raw = np.linspace(10.0, 100.0, 10).astype("<f8").tobytes()
    truncated_intensity_raw = full_intensity_raw[:-3]  # a few bytes short of the 10th float64
    truncated_intensity_b64 = base64.b64encode(truncated_intensity_raw).decode("ascii")

    stub = tmp_path / "mid_array_truncated.mzML"
    stub.write_bytes(_minimal_srm_mzml(time_b64, truncated_intensity_b64, q1=100.0, q3=50.0))

    chroms = read_srm_mzml(stub)  # must not raise

    assert (100.0, 50.0) in chroms
    chrom = chroms[(100.0, 50.0)]
    # 77 bytes // 8 = 9 whole float64s decoded from the truncated array. The time array
    # (fully intact at 10 points) is trimmed down to match.
    assert chrom.times.size == chrom.intensities.size == 9
    assert np.all(chrom.intensities >= 0)


def test_32bit_float_arrays_decode_correctly(tmp_path: Path):
    """32-bit float binaryDataArrays (dtype_bits=32) round-trip correctly to float64."""
    time_vals = np.array([0.0, 0.25, 0.5, 0.75, 1.0], dtype="<f4")
    intensity_vals = np.array([1.0, 2.5, 3.25, 4.125, 5.0], dtype="<f4")

    stub = tmp_path / "float32.mzML"
    stub.write_bytes(
        _minimal_srm_mzml(
            _encode_floats(time_vals, dtype_bits=32),
            _encode_floats(intensity_vals, dtype_bits=32),
            q1=200.0,
            q3=99.0,
            dtype_bits=32,
        )
    )

    chroms = read_srm_mzml(stub)

    chrom = chroms[(200.0, 99.0)]
    np.testing.assert_allclose(chrom.times, time_vals.astype(np.float64))
    np.testing.assert_allclose(chrom.intensities, intensity_vals.astype(np.float64))
