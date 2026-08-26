"""Read Agilent MRM/SRM chromatograms from msconvert-produced ``.mzML`` files.

The acquisition is targeted MRM: each ``.mzML`` holds a TIC plus a set of SRM
chromatograms, one per (precursor Q1 -> product Q3) transition, gated to retention-time
*segments*. Because Agilent splits a transition across time segments, the same nominal
(Q1, Q3) transition can appear in several ``<chromatogram>`` elements covering different
RT ranges. These are concatenated, in time order, into one continuous trace per transition.

Binary data arrays are decoded by hand (base64 -> optional zlib -> float32/float64) per
the ``cvParam`` accessions on each ``<binaryDataArray>`` rather than via a heavier
dependency, keeping the parse path transparent and dependency-free.
"""

from __future__ import annotations

import base64
import zlib
from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

_MZML_NS = "{http://psi.hupo.org/ms/mzml}"

_TIME_UNIT_SCALE_TO_MINUTES = {
    "minute": 1.0,
    "second": 1.0 / 60.0,
}


@dataclass(frozen=True, eq=False)
class Chromatogram:
    """A single MRM transition trace, concatenated across all its time segments.

    ``eq=False``: the default dataclass ``__eq__``/``__hash__`` compare fields with ``==``,
    which raises ``ValueError`` (ambiguous truth value) on the ``np.ndarray`` fields here.
    Identity-based equality (the ``object`` default) is required instead.
    """

    q1: float  # precursor (Q1) isolation-window target m/z, rounded to 1 dp
    q3: float  # product (Q3) isolation-window target m/z, rounded to 1 dp
    times: np.ndarray  # retention time, minutes, strictly increasing
    intensities: np.ndarray  # detector counts, same length as ``times``


def get_transition(
    chromatograms: dict[tuple[float, float], Chromatogram],
    q1: float,
    q3: float,
    tol: float = 0.3,
) -> Chromatogram | None:
    """Look up a transition by (Q1, Q3): exact rounded key first, then nearest within ``tol``.

    The canonical table's *labelled-standard* transitions are derived chemically (native
    precursor/product plus the nominal ``13C``/``d`` isotope shift) rather than read off the
    instrument, so they land ~0.1 Da from the acquired isolation-window targets: for example,
    the ``13C 3`` row's (200.0, 165.1) against an acquired (200.1, 165.1). All ten native
    homolog transitions match exactly. Only the labelled channels need the fallback, and
    without it every surrogate/IS channel would silently resolve to no chromatogram at all.

    Distance is the Manhattan sum ``|dQ1| + |dQ3|`` against the reader's rounded keys, and
    ties resolve to the last candidate in file order. ``tol`` is small enough that a real
    neighbouring transition (homologs are ~34 Da apart) can never be captured.
    """
    exact = chromatograms.get((round(q1, 1), round(q3, 1)))
    if exact is not None:
        return exact
    best: Chromatogram | None = None
    best_distance = tol
    for chrom in chromatograms.values():
        distance = abs(chrom.q1 - q1) + abs(chrom.q3 - q3)
        if distance <= best_distance:
            best, best_distance = chrom, distance
    return best


def _decode_binary(b64_text: str, compressed: bool, dtype_bits: int) -> np.ndarray:
    """Decode an mzML ``<binary>`` payload (base64 [+ zlib], float32 or float64) to float64.

    Decodes with a zero-copy ``np.frombuffer`` read (no boxed-float round trip through
    ``struct``) since this host is memory-constrained. If the decoded byte buffer is a few
    bytes short of a whole element (for example, an interrupted conversion truncated
    mid-float), the trailing partial bytes are dropped rather than raising.
    """
    raw = base64.b64decode(b64_text)
    if compressed:
        raw = zlib.decompress(raw)
    if not raw:
        return np.empty(0, dtype=np.float64)
    width = 4 if dtype_bits == 32 else 8
    n = len(raw) // width
    raw = raw[: n * width]  # drop a trailing partial element rather than misaligning the read
    return np.frombuffer(raw, dtype=f"<f{width}").astype(np.float64, copy=False)


def _isolation_target(node: ET.Element | None) -> float | None:
    """Extract the isolation-window target m/z from a ``<precursor>``/``<product>`` node."""
    if node is None:
        return None
    iso = node.find(f"{_MZML_NS}isolationWindow")
    if iso is None:
        return None
    for cv in iso.findall(f"{_MZML_NS}cvParam"):
        if cv.get("name") == "isolation window target m/z":
            value = cv.get("value")
            return float(value) if value else None
    return None


def _parse_arrays(chrom_el: ET.Element) -> dict[str, np.ndarray]:
    """Decode the time + intensity ``binaryDataArray`` children of one ``<chromatogram>``."""
    arrays: dict[str, np.ndarray] = {}
    time_scale = 1.0
    for bda in chrom_el.iter(f"{_MZML_NS}binaryDataArray"):
        compressed = False
        dtype_bits = 64
        kind: str | None = None
        for cv in bda.findall(f"{_MZML_NS}cvParam"):
            name = cv.get("name")
            if name == "zlib compression":
                compressed = True
            elif name == "32-bit float":
                dtype_bits = 32
            elif name == "64-bit float":
                dtype_bits = 64
            elif name == "time array":
                kind = "time"
                time_scale = _TIME_UNIT_SCALE_TO_MINUTES.get(cv.get("unitName", "minute"), 1.0)
            elif name == "intensity array":
                kind = "intensity"
        binary = bda.find(f"{_MZML_NS}binary")
        if kind is None or binary is None or not binary.text:
            continue
        arrays[kind] = _decode_binary(binary.text, compressed, dtype_bits)
    if "time" in arrays and time_scale != 1.0:
        arrays["time"] = arrays["time"] * time_scale
    return arrays


def read_srm_mzml(path: Path) -> dict[tuple[float, float], Chromatogram]:
    """Parse an mzML file into per-transition concatenated SRM chromatograms.

    Chromatograms without both a precursor and a product isolation window (the TIC, and
    any other non-SRM chromatogram) are skipped. Time-segmented chromatograms sharing the
    same rounded (Q1, Q3) are concatenated in time order, with duplicate timestamps at
    segment seams dropped (keeping the first occurrence), so each transition ends up as
    one strictly-increasing trace. Q1/Q3 are rounded to 1 decimal place because the raw
    isolation-window target m/z jitters slightly around the nominal transition value
    (acquisition noise, differing between segments and between runs). Rounding lets
    same-transition segments merge and lets the resulting key match the nominal values in the
    canonical peak table.

    Truncated or corrupt files (e.g. an interrupted msconvert run) are tolerated: whatever
    chromatograms parsed cleanly before the failure are returned rather than raising. This
    covers both a document truncated mid-XML (caught below) and one that is structurally
    complete but whose binary payload was truncated mid-element (handled in
    ``_decode_binary``, plus the length-mismatch guard below when only one of a
    chromatogram's two arrays was cut short).
    """
    path = Path(path)
    segments: dict[tuple[float, float], list[tuple[np.ndarray, np.ndarray]]] = {}

    try:
        for _event, el in ET.iterparse(str(path), events=("end",)):
            if el.tag != f"{_MZML_NS}chromatogram":
                continue
            q1 = _isolation_target(el.find(f"{_MZML_NS}precursor"))
            q3 = _isolation_target(el.find(f"{_MZML_NS}product"))
            if q1 is None or q3 is None:
                # TIC or other non-SRM chromatogram: nothing to quantify against.
                el.clear()
                continue
            arrays = _parse_arrays(el)
            t = arrays.get("time", np.empty(0))
            y = arrays.get("intensity", np.empty(0))
            if t.size != y.size:
                # One of the pair was truncated independently of the other (e.g. an
                # interrupted conversion cut off mid-element). Keep the shared, fully
                # paired prefix, a partial decode, rather than raising downstream.
                n = min(t.size, y.size)
                t, y = t[:n], y[:n]
            key = (round(q1, 1), round(q3, 1))
            segments.setdefault(key, []).append((t, y))
            el.clear()
    except (ET.ParseError, zlib.error, ValueError):
        # Interrupted conversion: keep whatever chromatograms parsed fully before this point.
        pass

    chromatograms: dict[tuple[float, float], Chromatogram] = {}
    for key, segs in segments.items():
        times = np.concatenate([s[0] for s in segs])
        intensities = np.concatenate([s[1] for s in segs])
        order = np.argsort(times, kind="mergesort")
        times, intensities = times[order], intensities[order]
        if times.size:
            keep = np.concatenate(([True], np.diff(times) > 0))
            times, intensities = times[keep], intensities[keep]
        chromatograms[key] = Chromatogram(q1=key[0], q3=key[1], times=times, intensities=intensities)

    return chromatograms
