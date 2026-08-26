"""Wrap the MSConvert (ProteoWizard) CLI, reproducing a fixed set of conversion settings.

Settings match a saved MSConvertGUI screenshot (``MSconvert/Msconvertsettings.png``): output
mzML, 64-bit binary encoding, zlib compression, write index, TPP compatibility, the default
``titleMaker`` filter, SRM kept as chromatograms (not spectra). These map to the CLI flags
below. MSConvert is third party and invoked as a subprocess.

:func:`is_truncated` and :func:`plan_conversion` answer the question a conversion step always
has to ask before it runs anything: *which* ``.D`` acquisitions in this run actually need
MSConvert this time. A run re-opened after a partial conversion (the host killed mid-write, or
a prior run just stopped early) should not have to reconvert everything, but it must never
trust a truncated output as if it were done. An interrupted MSConvert run can leave several
files in one folder as ~8 KB header-only mzML, which a size-blind check would accept.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree as ET

from quant209.core.sequence import FileEntry

# Flags reproducing a saved MSConvertGUI "Generic Defaults" + checkboxes.
LAB_MSCONVERT_SETTINGS: tuple[str, ...] = (
    "--mzML",
    "--64",  # 64-bit binary encoding precision
    "--zlib",  # zlib compression
    "--filter",
    "titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState> "
    'File:"<SourcePath>", NativeID:"<Id>"',  # default titleMaker. SRM stays as chromatograms
)

#: Below this size, an mzML is certainly not a complete SRM conversion. A real converted file
#: in this method runs well over 1 MB. This constant matches ``convert_file``'s own success
#: check, so "converted successfully" and "not truncated" agree.
_TRUNCATED_SIZE_FLOOR = 50_000

#: How many header bytes ``is_truncated`` probes for a parseable mzML/indexedmzML root. This
#: is generous over any real header (well before the chromatogram list) and cheap either way.
_HEADER_PROBE_BYTES = 65_536


@dataclass
class ConversionResult:
    source: Path
    output: Path | None
    ok: bool
    message: str


@dataclass
class ConversionPlan:
    """Which of a run's ``.D`` entries need MSConvert run, and why.

    Three disjoint buckets over the entries handed to :func:`plan_conversion`, in input
    order within each bucket:

    * ``to_convert``: no output at all yet in ``out_dir``.
    * ``reconvert``: an output exists but :func:`is_truncated` flags it (an interrupted prior
      conversion). It must be run again, not trusted.
    * ``skip``: a fresh, complete output already exists (or the entry is already ``.mzML``
      and needs no conversion at all).
    """

    to_convert: list[FileEntry] = field(default_factory=list)
    reconvert: list[FileEntry] = field(default_factory=list)
    skip: list[FileEntry] = field(default_factory=list)

    @property
    def pending(self) -> list[FileEntry]:
        """Every entry :func:`convert_many` should actually run, in ``to_convert`` then
        ``reconvert`` order, new work before repairs."""
        return self.to_convert + self.reconvert


def _output_path(source: Path, out_dir: Path) -> Path:
    return out_dir / f"{source.stem}.mzML"


def _header_parses(path: Path, probe_bytes: int = _HEADER_PROBE_BYTES) -> bool:
    """Whether the first ``probe_bytes`` of ``path`` look like a well-formed mzML header.

    Feeds only a bounded prefix into an incremental XML parser and never calls ``close()``,
    so a real file's inevitably-unterminated prefix (the probe stops short of the closing
    tags) parses cleanly. What would raise is a document that is not just incomplete but
    actually malformed (binary garbage, a corrupted write), which is exactly what
    "unparseable header" means to catch. Returns ``False`` for anything unreadable or
    genuinely invalid, and also if the probe never reaches an ``mzML``/``indexedmzML`` root
    start tag at all.
    """
    try:
        with path.open("rb") as fh:
            head = fh.read(probe_bytes)
    except OSError:
        return False
    parser = ET.XMLPullParser(events=("start",))
    try:
        parser.feed(head)
        saw_root = False
        for _event, elem in parser.read_events():
            tag = elem.tag.rsplit("}", 1)[-1]
            if tag in ("mzML", "indexedmzML"):
                saw_root = True
            elem.clear()
        return saw_root
    except ET.ParseError:
        return False


def is_truncated(path: str | Path) -> bool:
    """Whether ``path`` is a truncated/incomplete mzML rather than a finished conversion.

    True when the file is missing, smaller than a real converted file could plausibly be
    (:data:`_TRUNCATED_SIZE_FLOOR`), or its header is not parseable XML at all. A file that
    clears the size floor and parses is treated as complete. This is a cheap, filesystem-level
    check for planning purposes, not a full decode. ``core.pipeline.run_analysis`` still
    tolerates (and reports) a file that parses but decoded to no chromatogram data.
    """
    path = Path(path)
    try:
        size = path.stat().st_size
    except OSError:
        return True
    if size < _TRUNCATED_SIZE_FLOOR:
        return True
    return not _header_parses(path)


def plan_conversion(entries: list[FileEntry], out_dir: str | Path) -> ConversionPlan:
    """Classify a run's entries into what :func:`convert_many` actually needs to do.

    An entry already carrying a ``.mzML`` path needs no conversion and is filed under
    ``skip`` without touching ``out_dir`` at all. Every other entry (typically an Agilent
    ``.D`` acquisition directory) is checked against its would-be output in ``out_dir``
    (``core.pipeline``'s own ``_output_path`` convention: ``<stem>.mzML``): missing ->
    ``to_convert``, present but :func:`is_truncated` -> ``reconvert``, present and fresh ->
    ``skip``.

    Does not filter by :attr:`FileEntry.file_type`. A hexane-blank entry the caller (GUI or
    CLI) explicitly selected is planned exactly like a sample. Only quantification later
    declines to use it (see ``core.pipeline.QUANTIFIED_TYPES``).
    """
    out_dir = Path(out_dir)
    plan = ConversionPlan()
    for entry in entries:
        if entry.path.suffix.lower() == ".mzml":
            plan.skip.append(entry)
            continue
        output = _output_path(entry.path, out_dir)
        if not output.exists():
            plan.to_convert.append(entry)
        elif is_truncated(output):
            plan.reconvert.append(entry)
        else:
            plan.skip.append(entry)
    return plan


def convert_file(
    source: str | Path,
    msconvert_exe: str | Path,
    out_dir: str | Path,
    extra_args: tuple[str, ...] = (),
    timeout: float = 600.0,
) -> ConversionResult:
    """Convert one Agilent ``.D`` directory to ``.mzML`` in ``out_dir``."""
    source = Path(source)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    exe = Path(msconvert_exe)
    if not exe.exists():
        return ConversionResult(source, None, False, f"msconvert not found: {exe}")
    cmd = [str(exe), str(source), *LAB_MSCONVERT_SETTINGS, *extra_args, "-o", str(out_dir)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return ConversionResult(source, None, False, "conversion timed out")
    out = _output_path(source, out_dir)
    if proc.returncode == 0 and out.exists() and out.stat().st_size > 50_000:
        return ConversionResult(source, out, True, "ok")
    msg = (proc.stderr or proc.stdout or "unknown error").strip().splitlines()[-1:] or ["failed"]
    return ConversionResult(source, out if out.exists() else None, False, msg[0])


def convert_many(
    sources: list[str | Path],
    msconvert_exe: str | Path,
    out_dir: str | Path,
    extra_args: tuple[str, ...] = (),
    progress=None,
) -> list[ConversionResult]:
    """Convert several ``.D`` directories, optionally reporting progress via a callback."""
    results: list[ConversionResult] = []
    for i, src in enumerate(sources):
        if progress:
            progress(i, len(sources), Path(src).name)
        results.append(convert_file(src, msconvert_exe, out_dir, extra_args))
    return results
