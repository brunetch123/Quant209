"""End-to-end run orchestration: align -> anchor -> integrate -> quantify -> correct -> QA/QC.

This is the headless entry point the CLI and the Qt worker both call. It takes typed file
entries (``core.sequence.FileEntry``) and returns one :class:`RunResult` carrying every number
a report, a dashboard or an audit could ask for, including the evidence behind each one.

The order of operations, and why
--------------------------------

1. **The calibration injection is read first**, because everything downstream is measured
   against it. It supplies the run's retention-time map (:func:`core.calibrate.align_rt_map`),
   the anchors every file is then integrated at (:func:`core.calibrate.build_anchors`), and the
   answer to the run-fit guard (:func:`core.calibrate.rt_fit_quality`).

2. **Every quantified file is integrated one at a time**, its chromatograms released before the
   next is read. The host this runs on has ~2 GB of RAM and a run is ~30 files. Holding them
   all would not fit, and holding two is not necessary. Integration is a separate pass from
   quantification so that no result depends on file *order*. The calibration injections and
   the reference are just files in the list, and their peaks are already in hand by the time
   anything is computed from them.

3. **Response factors, then masses, then recoveries, then correction, then LOQ, then
   censoring.** Each stage's output is kept, not overwritten: ``masses``,
   ``masses_corrected`` and ``masses_censored`` are three separate frames, so an analyst
   asking "what did the recovery correction do to this congener?" can subtract two columns
   instead of re-running the pipeline.

Orchestration rules
--------------------
A few choices below are not obvious from the code, so they are stated here. The files used
for each are recorded on the result (``cal_files``, ``ref_files``) so a report can name them.

* **Retention-time anchors come from the FIRST calibration injection only** (``Cal01``), never
  averaged across replicates. An anchor is a position. Averaging positions across injections
  separated by a whole sequence would blur exactly the drift the anchors exist to track. The
  same file supplies the alignment and answers the fit guard.

* **Response factors are averaged over ALL calibration injections**, every one of
  them integrated at the anchors from the first. An RRF is a ratio of responses, and averaging
  four measurements of it is simply a better measurement. A congener mis-integrated in one
  injection is still quantified from the others (see
  :func:`core.quantify.average_calibrations`).

* **Recovery is measured against a single reference injection per sample**, never an average.
  The reference *defines* 100% recovery, so it is a denominator rather than a measurement to
  be improved. Averaging three would make each sample's recovery depend on injections it was
  never compared to. A run with no reference file is quantified without recovery correction
  and says so (``notes``). *Which* reference is chosen is covered below.

* **The LOQ is computed from the recovery-corrected blank masses**, and censoring is applied
  to the corrected masses. The limit has to live in the same units as the numbers it censors.

Bracketed reference injections
-------------------------------
**Each sample is paired with the reference injection that most recently preceded it in
acquisition order** (:func:`bracket_references`), which is what a bracketed sequence means and
what the analyst actually did. Using a single reference injection for the whole run is still
selectable as ``reference="first"``.

The evidence, from a three-batch run in the archive:

* The instrument sequence log runs ``Cal01, Ref04, [8 files], Cal02, Ref05, [7 files], Cal03,
  Ref06, [7 files], Cal04``: three analytical batches, each opened by its own calibration
  and reference injection.
* The analyst kept one workbook per batch, and each workbook's ``ref`` sheet holds that
  batch's own reference: the first workbook's surrogate masses are ``Ref04``'s, the second
  workbook's are ``Ref05``'s, the third workbook's are ``Ref06``'s.
* Nearest-preceding bracketing reproduces that pairing for all 22 quantified files, and moves
  the run's recovery error against the analyst from a mean absolute **0.066** to **0.023**.

Nearest-*preceding* rather than nearest-in-time: in that run a lab blank was acquired 14 h
after ``Ref04`` and 7 h before ``Ref05``, so nearest-in-time would pair it with ``Ref05``, and
the analyst paired it with ``Ref04``. A bracketing standard opens a batch. It does not
retroactively claim the batch before it. A file with no reference before it at all (a sequence
that opens with samples) falls back to the nearest one after it, which is the only reference
that can be said to bracket it.

A few more orchestration choices, all of them additive:

* The calibration and reference injections are **integrated and kept** like any other file.
  They are the two files a disputed result is most likely to be traced back to, so their peaks
  belong in the audit trail. Hexane carryover blanks are still excluded. The analyst never
  quantifies them.
* Every analyte an integration could not place is kept **as data** (``skipped_analytes``).
  ``integrate_run`` warns about them, and a warning is not an audit trail.
* The run-fit guard is surfaced as a run **status** rather than being advisory-and-ignored.
"""

from __future__ import annotations

import re
import warnings
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from .calibrate import (
    RtAlignment,
    RunFitReport,
    align_rt_map,
    build_anchors,
    reconcile_anchors,
    seed_from_standards,
    rt_fit_quality,
)
from .config import Analyte, load_canonical_peaks
from .integrate import integrate_run
from .mzml import read_srm_mzml
from .peaks import IntegrationParams, PeakResult
from .qaqc import QaqcParams, censor_below_loq, compute_loq, sample_flags
from .quantify import (
    CalibrationModel,
    average_calibrations,
    build_calibration,
    compute_recoveries,
    correct_for_recovery,
    internal_standard_areas,
    quantify_masses,
)
from .sequence import FileEntry
from .standards import HOMOLOG_SURROGATE

#: File types this pipeline integrates and quantifies. ``hex`` is deliberately absent: hexane
#: injections are carryover checks between samples and the analyst never quantifies them.
QUANTIFIED_TYPES = ("cal", "ref", "blank", "sample")

#: File types that are *results*: recovery-corrected against their own bracketing reference,
#: and censored against the LOQ. The calibration and reference injections are quantified.
#: Their masses are what the response factors and the recovery denominators are made of, but
#: they are not results, so no recovery is measured for them and their masses are reported
#: exactly as measured. This tuple is what scopes the censoring columns. See
#: :data:`CENSORED_ROLES` for the row scope.
REPORTED_TYPES = ("sample", "blank")

#: Analyte roles the LOQ is *applied* to. Censoring answers "is this congener present in this
#: sample above the laboratory background", which is a question about natives only.
#:
#: This is a named constant rather than an inline filter. A labelled standard is spiked at a
#: known amount into every file including the blanks, so its blank masses are identical by
#: construction. The internal standard quantifies to exactly ``IS_SPIKE_NG`` in every file,
#: giving the blanks zero spread, an LOQ equal to the spike itself, and therefore
#: ``mass <= LOQ`` true in every sample. Censoring the full frame zeroed 177 of 180
#: labelled-standard cells on one run, deleting the one column a reviewer uses to check that the
#: mass chain is wired the right way round.
CENSORED_ROLES = ("native",)

#: Run status values. ``FAILED-QC`` does not mean the run stopped. Every number was still
#: computed and is on the result. It means a run-level check failed and the numbers must not be
#: reported without a human looking at ``status_reason`` first.
STATUS_OK = "OK"
STATUS_FAILED_QC = "FAILED-QC"

#: How each file's recovery denominator is chosen. ``bracketed`` pairs a file with the
#: reference injection that most recently preceded it in acquisition order. This is the
#: default, and what a bracketed sequence means. ``first`` uses one reference for the whole
#: run and remains selectable. See the module docstring for the evidence.
REFERENCE_BRACKETED = "bracketed"
REFERENCE_FIRST = "first"
REFERENCE_POLICIES = (REFERENCE_BRACKETED, REFERENCE_FIRST)

#: Bytes of an mzML read when looking for its acquisition timestamp. ``startTimeStamp`` is an
#: attribute of the ``<run>`` element, which opens after the instrument/software header and
#: before the (very large) chromatogram list. On this dataset it lands inside the first few
#: kilobytes, and this is an order of magnitude of margin over that.
ACQUIRED_SCAN_BYTES = 262_144

_START_TIME = re.compile(rb'startTimeStamp="([^"]+)"')

#: ``progress(phase, i, n, label)``. ``i`` is 0-based and ``n`` is the count of steps in that
#: phase, so ``i / n`` is a fraction of the phase rather than of the run.
ProgressFn = Callable[[str, int, int, str], None]

#: Steps reported in the ``calibrate`` and ``quantify`` phases, in order.
_CALIBRATE_STEPS = ("read calibration injection", "retention-time alignment", "run-fit guard")
_QUANTIFY_STEPS = (
    "response factors",
    "congener masses",
    "surrogate recoveries",
    "recovery correction",
    "limit of quantification",
)


@dataclass(eq=False)
class FileRecord:
    """One data file's integration, kept whole so any result can be traced back to a trace.

    ``eq=False``: identity equality, matching the other result containers in this package.
    """

    entry: FileEntry
    #: The mzML actually read. Not derivable from ``entry.path``, which may be an Agilent
    #: ``.D`` directory that was converted into ``converted_dir``.
    mzml: Path
    peaks: list[PeakResult] = field(default_factory=list)
    #: Canonical analytes this file produced no peak for: no retention-time anchor, or the
    #: transition missing from a partially converted file. Kept as data because a run's output
    #: silently shrinking against the canonical table is the kind of loss that surfaces as a
    #: blank column in a report weeks later.
    skipped_analytes: list[str] = field(default_factory=list)
    #: This file's usable internal-standard areas (see
    #: :func:`core.quantify.internal_standard_areas`).
    is_areas: dict[str, float] = field(default_factory=dict)
    #: When the instrument acquired this file, from the mzML's ``startTimeStamp``. ``None``
    #: when the conversion carries no timestamp. It is what orders the run for reference
    #: bracketing, so its absence is reported rather than guessed around.
    acquired: datetime | None = None

    @property
    def sample_id(self) -> str:
        return self.entry.sample_id

    @property
    def file_type(self) -> str:
        return self.entry.file_type

    @property
    def usable(self) -> bool:
        """Whether this file decoded to anything integrable at all."""
        return bool(self.peaks)


@dataclass(eq=False)
class RunResult:
    """Everything a report, a dashboard or an audit needs from one analysis run.

    The frames are all **analyte-label by file-id** (``report.truth``'s orientation), so one
    congener's behaviour across the run is a single ``.loc[label]`` and the automated result
    joins to the analyst's ground truth without a transpose. Every frame carries the same row
    index, the full canonical analyte table, so an analyte that could not be quantified in
    a file reads as ``NaN`` there rather than vanishing: "not measurable" and "measured, and
    below the limit" (0.0, after censoring) are different report entries and must not be
    conflated.

    ``eq=False``: the generated ``__eq__`` would compare DataFrame fields with ``==``, which
    returns a DataFrame rather than a bool and raises. Identity equality is what is wanted, the
    same defect class as ``core.mzml.Chromatogram`` and ``report.truth.RunTruth``.
    """

    # ---- what was run ------------------------------------------------------------------
    #: Every entry handed in, including the hexane blanks that were not quantified.
    entries: list[FileEntry]
    #: The entries actually integrated, in the order they were processed.
    quantified_entries: list[FileEntry]
    analytes: list[Analyte]
    params: IntegrationParams
    qaqc: QaqcParams

    # ---- how the run was calibrated ------------------------------------------------------
    #: Names of the calibration injections, first-listed first. ``cal_files[0]`` supplied the
    #: alignment and the anchors. All of them contributed response factors.
    cal_files: list[str]
    #: Names of the reference injections, in the order they were listed.
    ref_files: list[str]
    alignment: RtAlignment
    anchors: dict[str, float]
    fit: RunFitReport
    calibration: CalibrationModel

    # ---- verdict --------------------------------------------------------------------------
    status: str = STATUS_OK
    #: Analyst-readable reason when ``status`` is not OK, empty otherwise.
    status_reason: str = ""
    #: Run-level observations that are not failures (a missing reference injection, a file
    #: that decoded to nothing). Written for the person reading the report.
    notes: list[str] = field(default_factory=list)

    # ---- per-file evidence ---------------------------------------------------------------
    #: File id -> its whole integration. Keyed by ``FileEntry.sample_id``.
    files: dict[str, FileRecord] = field(default_factory=dict)
    sample_ids: list[str] = field(default_factory=list)
    blank_ids: list[str] = field(default_factory=list)
    cal_ids: list[str] = field(default_factory=list)
    ref_ids: list[str] = field(default_factory=list)

    # ---- results --------------------------------------------------------------------------
    #: Integrated peak areas (counts·s).
    areas: pd.DataFrame = field(default_factory=pd.DataFrame)
    #: The presence test's verdict per (analyte, file), plain ``bool``, same shape as
    #: ``areas``. This is what a report, the GUI or a QA/QC rule gates on, never ``sn``, which
    #: an absent peak reports large and meaningless (see ``core.peaks.PeakResult.sn``).
    #: Quantification does *not* apply it (see ``core.quantify.quantify_masses``): the analyst
    #: integrates and quantifies a shoulder they chose to split off a dominant neighbour, and
    #: so does this chain, so the verdict travels beside the mass instead of silently zeroing
    #: it. ``False`` also covers "this file produced no peak for the analyte at all".
    #: ``skipped_analytes`` is where that distinction is kept.
    detected: pd.DataFrame = field(default_factory=pd.DataFrame)
    #: ``RRF * area * (IS_spike / IS_area)``, before recovery correction. This is the frame
    #: that joins to the analyst workbook's "Mass (ng)" column.
    masses: pd.DataFrame = field(default_factory=pd.DataFrame)
    #: Masses after one-sided surrogate-recovery correction (unchanged for cal/ref files).
    masses_corrected: pd.DataFrame = field(default_factory=pd.DataFrame)
    #: Corrected masses with at-or-below-LOQ cells set to 0. **The reported numbers.**
    masses_censored: pd.DataFrame = field(default_factory=pd.DataFrame)
    #: Homolog (mono..deca) by file id. Only the reported files carry recoveries.
    recoveries: pd.DataFrame = field(default_factory=pd.DataFrame)
    #: Per-analyte limit of quantification (ng), from this run's usable blanks. Reported for
    #: every analyte, but **applied** only to native analytes in reported files (see
    #: :data:`CENSORED_ROLES` / :data:`REPORTED_TYPES`). A labelled standard's entry here is
    #: arithmetic rather than a limit. It is spiked at a known amount into every blank, so
    #: its blank spread is zero and its "LOQ" is its own spike amount.
    loq: pd.Series = field(default_factory=lambda: pd.Series(dtype=float))
    #: File id -> QA/QC flags (see ``core.qaqc.sample_flags``).
    flags: dict[str, list[str]] = field(default_factory=dict)
    #: How recovery denominators were chosen: one of :data:`REFERENCE_POLICIES`.
    reference_policy: str = REFERENCE_BRACKETED
    #: File id -> the reference-injection file id its recoveries were measured against. Every
    #: recovery in the run traces to one row of this mapping, which is the point of keeping it.
    reference_for: dict[str, str] = field(default_factory=dict)
    #: Which calibration injection supplied the retention-time anchors each file was
    #: integrated at, keyed by the file's name on disk. Retention times drift across a
    #: sequence, so this is not always the first injection. See
    #: :func:`nearest_calibration`. Persisted because an integration cannot be checked
    #: without knowing the retention times it was performed at.
    calibration_for: dict[str, str] = field(default_factory=dict)

    @property
    def skipped_analytes(self) -> dict[str, list[str]]:
        """File id -> the analytes that file's integration could not place."""
        return {name: record.skipped_analytes for name, record in self.files.items()}

    @property
    def reported_ids(self) -> list[str]:
        """Every file that is a *result*: samples then blanks, in processing order."""
        return self.sample_ids + self.blank_ids

    @property
    def ok(self) -> bool:
        return self.status == STATUS_OK


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------


def _mzml_path(entry: FileEntry, converted_dir: Path) -> Path:
    """The mzML to read for an entry: itself, or its msconvert output in ``converted_dir``.

    Conversion is not performed here. The GUI and the CLI convert before calling (see
    ``quant209.convert``), and a pipeline that silently shelled out to msconvert mid-run would
    make "why is this run slow" and "why did this run fail" the same question.
    """
    if entry.path.suffix.lower() == ".mzml":
        return entry.path
    return converted_dir / f"{entry.path.stem}.mzML"


def _unique_id(entry: FileEntry, taken: set[str]) -> str:
    """``entry.sample_id``, or the file stem if two files claim the same id.

    Sample ids are derived from filenames and are unique in every observed sequence, but they
    are not unique *by construction*. A re-injection saved under a different date would
    collide. Silently overwriting one file's results with another's is the worst available
    outcome, so the second file falls back to its full stem.
    """
    sample_id = entry.sample_id
    if sample_id not in taken:
        return sample_id
    stem = entry.path.stem
    if stem not in taken:
        return stem
    n = 2
    while f"{stem}#{n}" in taken:
        n += 1
    return f"{stem}#{n}"


def read_acquired(path: Path) -> datetime | None:
    """The instrument acquisition time from an mzML's ``startTimeStamp``, or None.

    A bounded header read (:data:`ACQUIRED_SCAN_BYTES`) rather than a parse: the whole point
    is to order a run's files without paying to decode any of them, and the attribute sits in
    the ``<run>`` element well before the chromatogram data. Anything unreadable, absent or
    unparseable returns ``None``, which the caller reports rather than guesses around.
    """
    try:
        with path.open("rb") as handle:
            head = handle.read(ACQUIRED_SCAN_BYTES)
    except OSError:
        return None
    match = _START_TIME.search(head)
    if match is None:
        return None
    try:
        # mzML writes UTC as a trailing 'Z', which fromisoformat accepts from 3.11.
        stamp = datetime.fromisoformat(match.group(1).decode("ascii").replace("Z", "+00:00"))
    except (UnicodeDecodeError, ValueError):
        return None
    # A converter that omitted the offset yields a naive datetime, and one naive stamp in a run
    # of aware ones makes every ordering comparison in :func:`bracket_references` raise. A
    # single odd file would abort a whole run rather than degrade. mzML's startTimeStamp is
    # specified as UTC, so a missing offset is treated as UTC rather than guessed.
    return stamp if stamp.tzinfo is not None else stamp.replace(tzinfo=timezone.utc)



def nearest_calibration(
    acquired: Mapping[str, datetime | None],
    cal_names: Sequence[str],
    target_names: Iterable[str],
) -> dict[str, str]:
    """Pair each file with the calibration injection nearest it in acquisition time.

    Retention times drift as a sequence runs. The column conditions, the inlet and the oven
    are not identical eight hours in, which is why the lab injects the calibration standard
    several times across a run rather than once at the start. On one 2024 run the drift from
    the first calibration to the last is a median of 0.29 min and as much as 1.03 min for
    individual congeners, against an integration territory of ±0.20 min: anchoring the whole
    run to the first injection therefore puts the later files' peaks outside the window they
    are searched in, and they are reported absent. Ten of that run's twenty-four files lost
    the PCB204 internal standard that way, and with it every high-chlorination mass in those
    files.

    *Nearest* rather than *preceding*: unlike :func:`bracket_references`, where the reference
    standard opens a batch and defines its recovery denominator. Nothing opens a batch here:
    drift is continuous, so the best estimate of a file's retention times is the calibration
    injection closest to it in time, whichever side it falls on.

    Files or calibrations with no acquisition timestamp fall back to ``cal_names[0]``, and
    the caller records the pairing so a reviewer can see which calibration produced the
    anchors any given file was integrated at.
    """
    if not cal_names:
        return {}
    fallback = cal_names[0]
    timed = sorted(
        ((acquired[c], c) for c in cal_names if acquired.get(c) is not None),
        key=lambda pair: pair[0],
    )
    if not timed:
        return {target: fallback for target in target_names}
    pairing: dict[str, str] = {}
    for target in target_names:
        when = acquired.get(target)
        if when is None:
            pairing[target] = fallback
            continue
        pairing[target] = min(timed, key=lambda pair: abs(pair[0] - when))[1]
    return pairing



#: How far one file's internal-standard area may sit from the run's own median, as a
#: multiplying factor either way, before it is treated as a mis-integration.
#:
#: The internal standard is spiked *identically* into every sample and blank in a run, so its
#: area is the one number in the file that is known in advance to be the same as its
#: neighbours'. A file that disagrees by more than this is not a sample with unusual
#: chemistry. It is a peak that was integrated in the wrong place. That matters more than
#: any other single number, because every congener sharing that standard is scaled by
#: ``IS_spike / IS_area``: on a 2025 run one file's PCB204 integrated to 66
#: counts·s against a run median of 65,704, which would multiply that file's whole
#: high-chlorination profile by a thousand.
#:
#: The range is wide. Real between-sample variation in an internal standard is tens of
#: percent. This is looking for order-of-magnitude nonsense, and a threshold tight enough to
#: argue about would start discarding good data.
IS_CONSISTENCY_RANGE = (0.2, 5.0)


def reject_inconsistent_standards(
    records: Mapping[str, "FileRecord"], reported_ids: Sequence[str]
) -> list[str]:
    """Drop internal-standard areas that disagree with the run's own median. Returns notes.

    Compared across the *reported* files only: samples and blanks. Calibration and
    reference injections carry the standard at a different level by design (on
    one run the calibration injections' internal standards run about seven
    times the samples'), so including them would move the median toward a population the
    samples are not part of.

    Dropping the area, rather than correcting it, is the whole point: a congener normalised
    by a standard that was never measured is not a quantity, and the pipeline reports it as
    unquantified instead of as a plausible wrong number.
    """
    notes: list[str] = []
    if len(reported_ids) < 3:
        return notes  # too few files for a median to mean anything
    labels = {label for fid in reported_ids for label in records[fid].is_areas}
    low, high = IS_CONSISTENCY_RANGE
    for label in sorted(labels):
        areas = {
            fid: records[fid].is_areas[label]
            for fid in reported_ids
            if label in records[fid].is_areas
        }
        if len(areas) < 3:
            continue
        median = float(np.median(list(areas.values())))
        if median <= 0:
            continue
        for fid, area in areas.items():
            if area < low * median or area > high * median:
                del records[fid].is_areas[label]
                notes.append(
                    f"{fid}: internal standard {label} integrated to {area:,.0f} against a run "
                    f"median of {median:,.0f}, so it is not used. Congeners normalised by it are "
                    f"reported unquantified"
                )
    return notes


def bracket_references(
    acquired: Mapping[str, datetime | None],
    ref_ids: Sequence[str],
    target_ids: Iterable[str],
) -> dict[str, str]:
    """Pair each target file with the reference injection that brackets it.

    "Brackets" means the nearest reference injection *preceding* the file in acquisition
    order, which is what opening a batch with a reference standard is for. A file with no
    reference before it falls back to the earliest one after it, the only reference that can
    be said to bracket a sequence that opens with samples.

    Files or references with no acquisition time are paired with ``ref_ids[0]``: an unordered
    file cannot be bracketed, and the caller reports the fallback rather than letting it pass
    as a bracketing decision. Returns ``{}`` when the run has no reference injections at all.

    The function is pure, depending only on timestamps. The choice of denominator is the
    single orchestration decision that most changes a reported recovery, so it can be checked
    against a sequence log without reading a chromatogram.
    """
    if not ref_ids:
        return {}
    fallback = ref_ids[0]
    timed = sorted(
        ((acquired[r], r) for r in ref_ids if acquired.get(r) is not None),
        key=lambda pair: pair[0],
    )
    if not timed:
        return {target: fallback for target in target_ids}

    pairing: dict[str, str] = {}
    for target in target_ids:
        when = acquired.get(target)
        if when is None:
            pairing[target] = fallback
            continue
        preceding = [ref for stamp, ref in timed if stamp <= when]
        pairing[target] = preceding[-1] if preceding else timed[0][1]
    return pairing


def _frame(
    per_file: dict[str, dict[str, float]], index: Sequence[str], columns: Sequence[str]
) -> pd.DataFrame:
    """``{file id: {label: value}}`` -> a label-by-file frame on the full analyte index."""
    frame = pd.DataFrame(
        {name: pd.Series(values, dtype=float) for name, values in per_file.items()}
    )
    return frame.reindex(index=list(index), columns=list(columns))


def _bool_frame(
    per_file: dict[str, dict[str, bool]], index: Sequence[str], columns: Sequence[str]
) -> pd.DataFrame:
    """``{file id: {label: verdict}}`` -> a label-by-file frame of plain ``bool``.

    Missing cells become ``False`` rather than NaN: this frame's whole purpose
    is to be used as a mask, and an analyte with no peak at all was certainly not detected as
    one. A nullable dtype here would make ``frame[mask]`` a source of three-valued-logic bugs
    in every caller for the sake of a distinction that ``skipped_analytes`` already records.
    """
    frame = pd.DataFrame(
        {name: pd.Series(values, dtype=object) for name, values in per_file.items()}
    )
    return frame.reindex(index=list(index), columns=list(columns)).fillna(False).astype(bool)


def _integrate_one(
    record: FileRecord, analytes: list[Analyte], anchors: dict[str, float], params: IntegrationParams
) -> None:
    """Read, integrate and release one file's chromatograms, filling ``record`` in place.

    Chromatograms are dropped before returning: a run is ~30 files on a host with ~2 GB, and
    the peaks are all that is needed afterwards. ``integrate_run``'s warning about unplaced
    analytes is suppressed and recomputed as data on the record: the same fact, kept
    somewhere a report can read it.
    """
    if not record.mzml.exists():
        record.skipped_analytes = [a.label for a in analytes]
        return
    record.acquired = read_acquired(record.mzml)
    chroms = read_srm_mzml(record.mzml)
    try:
        if not any(c.times.size > 1 for c in chroms.values()):
            # A truncated or header-only conversion. ``read_srm_mzml`` tolerates it and
            # returns whatever prefix it decoded, which can be nothing at all.
            record.skipped_analytes = [a.label for a in analytes]
            return
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            record.peaks = integrate_run(chroms, analytes, anchors, params)
    finally:
        del chroms
    placed = {p.label for p in record.peaks}
    record.skipped_analytes = [a.label for a in analytes if a.label not in placed]
    record.is_areas = internal_standard_areas(record.peaks, analytes)


# --------------------------------------------------------------------------------------
# The run
# --------------------------------------------------------------------------------------


def run_analysis(
    entries: list[FileEntry],
    converted_dir: str | Path,
    params: IntegrationParams | None = None,
    progress: ProgressFn | None = None,
    qaqc: QaqcParams | None = None,
    analytes: list[Analyte] | None = None,
    reference: str = REFERENCE_BRACKETED,
) -> RunResult:
    """Process one analysis run from typed file entries to QA/QC'd, auditable results.

    ``entries`` is the run's files with their detected types (``core.sequence.scan_inputs``).
    ``converted_dir`` is where an Agilent ``.D`` entry's msconvert output lives (ignored for
    entries that are already mzML). ``progress(phase, i, n, label)`` is called once per file
    in the ``integrate`` phase and once per step in the ``calibrate`` and ``quantify`` phases.
    ``reference`` selects the recovery denominator policy (:data:`REFERENCE_POLICIES`). See
    the module docstring for why bracketing is the default.

    Raises ``ValueError`` if the run carries no calibration injection: without one there are no
    retention-time anchors and no response factors, so there is nothing to compute and no
    partial answer worth returning.

    **A failing run-fit guard does not stop the run.** If the retention-time table does not
    describe this run's calibration injection, every number is still computed. An analyst
    looking at a flagged run is exactly who can tell a GC method change from a bad injection,
    and needs the numbers to tell them apart. However, ``status`` becomes ``FAILED-QC`` and
    ``status_reason`` carries the guard's own analyst-readable explanation for the GUI to
    display. Silently continuing, or silently refusing, would both be worse.
    """
    params = params or IntegrationParams()
    qaqc = qaqc or QaqcParams()
    analytes = analytes if analytes is not None else load_canonical_peaks()
    converted_dir = Path(converted_dir)
    labels = [a.label for a in analytes]
    if reference not in REFERENCE_POLICIES:
        raise ValueError(
            f"unknown reference policy {reference!r}; expected one of {REFERENCE_POLICIES}"
        )

    def report(phase: str, i: int, n: int, label: str) -> None:
        if progress is not None:
            progress(phase, i, n, label)

    cal_entries = [e for e in entries if e.file_type == "cal"]
    if not cal_entries:
        raise ValueError(
            "this run has no calibration injection: retention-time anchors and response "
            "factors both come from one, so nothing can be quantified"
        )
    ref_entries = [e for e in entries if e.file_type == "ref"]

    # ---- 1) Calibration injection: alignment, anchors, fit guard -----------------------
    # All three come from the FIRST cal injection (see the module docstring). Its
    # chromatograms are released here and the file is re-read in the integration loop like
    # any other: one extra read of one file, in exchange for a flat memory profile and a
    # single code path that integrates every file identically.
    first_cal = cal_entries[0]
    cal_mzml = _mzml_path(first_cal, converted_dir)
    report("calibrate", 0, len(_CALIBRATE_STEPS), f"{_CALIBRATE_STEPS[0]}: {first_cal.name}")
    if not cal_mzml.exists():
        raise ValueError(f"calibration injection not found: {cal_mzml}")
    cal_chroms = read_srm_mzml(cal_mzml)
    try:
        if not any(c.times.size > 1 for c in cal_chroms.values()):
            raise ValueError(
                f"calibration injection {cal_mzml.name} decoded to no chromatograms "
                "(truncated or header-only conversion)"
            )
        report("calibrate", 1, len(_CALIBRATE_STEPS), _CALIBRATE_STEPS[1])
        # Put the analyte table on this run's own clock first, using the labelled standards,
        # so the per-channel alignment refines in the right neighbourhood instead of
        # searching for one (see core.calibrate.seed_from_standards).
        seeded = seed_from_standards(cal_chroms, analytes, params)
        alignment = align_rt_map(cal_chroms, seeded, params)
        anchors = build_anchors(cal_chroms, seeded, params, alignment=alignment)
        report("calibrate", 2, len(_CALIBRATE_STEPS), _CALIBRATE_STEPS[2])
        fit = rt_fit_quality(cal_chroms, seeded, anchors, params)
    finally:
        del cal_chroms

    # Every calibration injection gets its own anchors, so that a file can be integrated at
    # the retention times of the calibration nearest it rather than at the first one's. See
    # :func:`nearest_calibration` for the drift this corrects. Each injection is read,
    # anchored and released in turn: only the anchor tables are kept.
    anchors_by_cal: dict[str, dict[str, float]] = {first_cal.name: anchors}
    for extra in cal_entries[1:]:
        extra_mzml = _mzml_path(extra, converted_dir)
        if not extra_mzml.exists():
            continue
        extra_chroms = read_srm_mzml(extra_mzml)
        try:
            if not any(c.times.size > 1 for c in extra_chroms.values()):
                continue
            extra_seeded = seed_from_standards(extra_chroms, analytes, params)
            extra_alignment = align_rt_map(extra_chroms, extra_seeded, params)
            anchors_by_cal[extra.name] = build_anchors(
                extra_chroms, extra_seeded, params, alignment=extra_alignment
            )
        finally:
            del extra_chroms

    # Each calibration was anchored on its own, so each could independently slip a window
    # by one peak. Comparing them with each other catches that. See
    # :func:`core.calibrate.reconcile_anchors`. Every repair is recorded rather than applied
    # silently.
    anchors_by_cal, anchor_repairs = reconcile_anchors(
        anchors_by_cal, [e.name for e in cal_entries if e.name in anchors_by_cal]
    )
    anchors = anchors_by_cal.get(first_cal.name, anchors)

    notes: list[str] = list(anchor_repairs)

    # ---- 2) Integrate every quantified file, one at a time ------------------------------
    quantified = [e for e in entries if e.file_type in QUANTIFIED_TYPES]
    records: dict[str, FileRecord] = {}
    ids_by_type: dict[str, list[str]] = {t: [] for t in QUANTIFIED_TYPES}

    # Which calibration's anchors each file is integrated at. Keyed by entry *name* (the
    # file on disk) rather than by file id, because the pairing is a fact about acquisition
    # order and is decided before any file is read.
    acquired_by_name = {
        e.name: read_acquired(_mzml_path(e, converted_dir))
        for e in quantified
    }
    calibration_for = nearest_calibration(
        acquired_by_name,
        [e.name for e in cal_entries if e.name in anchors_by_cal],
        [e.name for e in quantified],
    )

    for i, entry in enumerate(quantified):
        report("integrate", i, len(quantified), entry.sample_id)
        file_id = _unique_id(entry, set(records))
        record = FileRecord(entry=entry, mzml=_mzml_path(entry, converted_dir))
        file_anchors = anchors_by_cal.get(calibration_for.get(entry.name, ""), anchors)
        _integrate_one(record, analytes, file_anchors, params)
        if not record.usable:
            notes.append(
                f"{entry.name}: decoded to no chromatograms (truncated or missing "
                "conversion). Integrated nothing, quantified nothing"
            )
        records[file_id] = record
        ids_by_type[entry.file_type].append(file_id)

    cal_ids, ref_ids = ids_by_type["cal"], ids_by_type["ref"]
    sample_ids, blank_ids = ids_by_type["sample"], ids_by_type["blank"]
    file_ids = list(records)

    # An internal standard is spiked identically into every reported file, so one that
    # disagrees with its own run is a mis-integration, and it would scale every congener
    # that shares it. Rejected here, before any mass is computed.
    notes.extend(reject_inconsistent_standards(records, sample_ids + blank_ids))

    # ---- 3) Response factors, averaged over every usable calibration injection -----------
    report("quantify", 0, len(_QUANTIFY_STEPS), _QUANTIFY_STEPS[0])
    models = [
        build_calibration(records[cid].peaks, analytes, name=records[cid].entry.name)
        for cid in cal_ids
        if records[cid].usable
    ]
    if not models:
        raise ValueError(
            "no calibration injection in this run produced integrable chromatograms, "
            "so no response factors could be derived"
        )
    calibration = average_calibrations(models)

    # ---- 4) Masses, for every file (cal and ref included: they are what 5 is made of) ----
    report("quantify", 1, len(_QUANTIFY_STEPS), _QUANTIFY_STEPS[1])
    areas: dict[str, dict[str, float]] = {}
    detected: dict[str, dict[str, bool]] = {}
    masses: dict[str, dict[str, float]] = {}
    for file_id, record in records.items():
        areas[file_id] = {p.label: float(p.area) for p in record.peaks}
        detected[file_id] = {p.label: bool(p.detected) for p in record.peaks}
        masses[file_id] = quantify_masses(
            record.peaks, analytes, calibration, record.is_areas
        )

    # ---- 5) Recoveries against each file's own bracketing reference injection ------------
    report("quantify", 2, len(_QUANTIFY_STEPS), _QUANTIFY_STEPS[2])
    reported_ids = sample_ids + blank_ids
    usable_refs = [rid for rid in ref_ids if records[rid].usable]
    reference_for: dict[str, str] = {}
    if not usable_refs:
        notes.append(
            (
                f"reference injection(s) {', '.join(e.name for e in ref_entries)} produced no "
                "usable masses"
                if ref_entries
                else "this run has no reference injection"
            )
            + ": surrogate recovery cannot be measured, so masses are reported uncorrected"
        )
    elif reference == REFERENCE_BRACKETED:
        acquired = {fid: rec.acquired for fid, rec in records.items()}
        reference_for = bracket_references(acquired, usable_refs, reported_ids)
        # Two different losses, reported separately: a file with no timestamp cannot be
        # placed in the sequence, and a *reference* with no timestamp cannot bracket
        # anything and drops out of the candidate set for every file.
        untimed_refs = sorted(rid for rid in usable_refs if records[rid].acquired is None)
        untimed_files = sorted(fid for fid in reported_ids if records[fid].acquired is None)
        if untimed_refs:
            notes.append(
                f"{len(untimed_refs)} reference injection(s) carry no mzML acquisition "
                "timestamp and were excluded from bracketing: " + ", ".join(untimed_refs)
            )
        if untimed_files:
            notes.append(
                f"{len(untimed_files)} file(s) carry no mzML acquisition timestamp and could "
                f"not be bracketed; they fall back to {records[usable_refs[0]].entry.name}: "
                + ", ".join(untimed_files)
            )
    else:
        reference_for = {fid: usable_refs[0] for fid in reported_ids}

    recoveries = {
        file_id: compute_recoveries(
            masses[file_id],
            masses[reference_for[file_id]] if file_id in reference_for else {},
            (
                {p.label: p.detected for p in records[reference_for[file_id]].peaks}
                if file_id in reference_for
                else None
            ),
        )
        for file_id in reported_ids
    }

    # ---- 6) One-sided recovery correction -------------------------------------------------
    report("quantify", 3, len(_QUANTIFY_STEPS), _QUANTIFY_STEPS[3])
    corrected = {
        file_id: correct_for_recovery(
            per_file, recoveries.get(file_id, {}), analytes
        )
        for file_id, per_file in masses.items()
    }

    areas_frame = _frame(areas, labels, file_ids)
    masses_frame = _frame(masses, labels, file_ids)
    corrected_frame = _frame(corrected, labels, file_ids)
    detected_frame = _bool_frame(detected, labels, file_ids)

    # ---- 7) LOQ from the blanks, then censoring ------------------------------------------
    # Only blanks that actually decoded count towards the minimum: a file that produced no
    # chromatograms contributes a column of zeros, which is not a measurement of the
    # laboratory background and must not be what lets a run clear ``min_blanks_for_loq``.
    report("quantify", 4, len(_QUANTIFY_STEPS), _QUANTIFY_STEPS[4])
    loq_blank_ids = [bid for bid in blank_ids if records[bid].usable]
    loq = compute_loq(corrected_frame[loq_blank_ids], qaqc)
    if loq.empty and len(loq_blank_ids) < qaqc.min_blanks_for_loq:
        notes.append(
            f"{len(loq_blank_ids)} usable blank file(s): fewer than the "
            f"{qaqc.min_blanks_for_loq} needed to measure a blank spread, so no LOQ was "
            "computed and nothing was censored"
        )

    # Censoring is scoped to native analytes in reported files (see CENSORED_ROLES and
    # REPORTED_TYPES) and the untouched rest of the frame is carried through unchanged, so
    # ``masses_censored`` stays the same shape as every other frame and differs from
    # ``masses_corrected`` only where a censoring decision was actually taken.
    censor_rows = [a.label for a in analytes if a.role in CENSORED_ROLES]
    censor_cols = [fid for fid in file_ids if records[fid].file_type in REPORTED_TYPES]
    censored_frame = corrected_frame.copy()
    if censor_rows and censor_cols:
        censored_frame.loc[censor_rows, censor_cols] = censor_below_loq(
            corrected_frame.loc[censor_rows, censor_cols], loq, qaqc
        )

    flags = {
        file_id: sample_flags(recoveries.get(file_id, {}), records[file_id].peaks, qaqc)
        for file_id in file_ids
    }

    return RunResult(
        entries=list(entries),
        quantified_entries=list(quantified),
        analytes=list(analytes),
        params=params,
        qaqc=qaqc,
        cal_files=[e.name for e in cal_entries],
        ref_files=[e.name for e in ref_entries],
        alignment=alignment,
        anchors=anchors,
        fit=fit,
        calibration=calibration,
        status=STATUS_OK if fit.fit_ok else STATUS_FAILED_QC,
        status_reason="" if fit.fit_ok else f"{fit.reason} ({fit.summary})",
        notes=notes,
        files=records,
        sample_ids=sample_ids,
        blank_ids=blank_ids,
        cal_ids=cal_ids,
        ref_ids=ref_ids,
        areas=areas_frame,
        detected=detected_frame,
        masses=masses_frame,
        masses_corrected=corrected_frame,
        masses_censored=censored_frame,
        recoveries=_frame(recoveries, list(HOMOLOG_SURROGATE), reported_ids),
        loq=loq,
        flags=flags,
        reference_policy=reference,
        reference_for=reference_for,
        calibration_for=calibration_for,
    )
