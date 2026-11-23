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
from .deconvolve import fit_unresolved_pair
from .ghosts import SCALE_MIN_R, Ghost, GhostScale, integrate_with_ghosts, is_fraction, run_scale
from .integrate import integrate_run, unresolved_pairs
from .mzml import get_transition, read_srm_mzml
from .peaks import IntegrationParams, PeakResult, local_maxima, smooth
from .qaqc import QaqcParams, censor_below_loq, compute_loq, sample_flags
from .quantify import (
    CalibrationModel,
    IS_ROLES,
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
    #: For each internal standard, the tallest point on its channel within
    #: :data:`IS_PRESENCE_WINDOW` of where it was searched for, smoothed. Read only for a
    #: standard that was not detected, to tell one that is not in the injection at all from
    #: one that is there but was not found where expected (:func:`classify_missing_standards`).
    is_window_max: dict[str, float] = field(default_factory=dict)
    #: When the instrument acquired this file, from the mzML's ``startTimeStamp``. ``None``
    #: when the conversion carries no timestamp. It is what orders the run for reference
    #: bracketing, so its absence is reported rather than guessed around.
    acquired: datetime | None = None
    #: Every ghost peak found in this file and kept out of its integration
    #: (:mod:`core.ghosts`): whose copy it is, where, how large, and whether it was fenced off
    #: or subtracted. Empty when the correction is off.
    ghosts: list[Ghost] = field(default_factory=list)

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
    #: Calibration injections whose retention-time table did not fit, keyed by name, with the
    #: fit check's reason. They anchor no other file and contribute no response factor (see
    #: :func:`select_calibrations`). Empty when every injection fits, or when none does.
    excluded_calibrations: dict[str, str] = field(default_factory=dict)
    #: Neighbouring congeners this run's calibration injections could not separate, as
    #: ``(earlier, later)`` label pairs (see :func:`core.integrate.unresolved_pairs`). Wherever
    #: such a pair shares one envelope it is split by a two-peak fit, in every file. Empty when
    #: the calibration resolves every pair.
    unresolved_pairs: tuple[tuple[str, str], ...] = ()
    #: The internal standard 204's ghost-peak fraction this run's correction was scaled by
    #: (:class:`core.ghosts.GhostScale`). ``None`` when ``params.ghost_correction`` is off.
    ghost_scale: GhostScale | None = None

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

#: Two internal standards that leave ``IS_CONSISTENCY_RANGE`` together, within this factor of
#: each other, moved because the whole injection did: a diluted or concentrated extract, or a
#: different injection volume. Both standards and every analyte scale alike, so the
#: analyte-to-standard ratio that every mass is computed from is unchanged, and the standards
#: are used as measured. A mis-integrated standard moves alone, by an order of magnitude or more.
#: On the archive the diluted samples' two standards agree within 6%.
IS_DILUTION_AGREEMENT = 1.5

#: The furthest a file's internal standards may move together and still be read as a dilution
#: or concentration. Beyond it the standards have not been measured: on the archive the two
#: files whose standards came out at a thousandth of the run median were integration failures.
IS_DILUTION_RANGE = (0.01, 100.0)


def reject_inconsistent_standards(
    records: Mapping[str, "FileRecord"], reported_ids: Sequence[str]
) -> list[str]:
    """Drop internal-standard areas that disagree with the run's own median. Returns notes.

    Compared across the *reported* files only: samples and blanks. Calibration and
    reference injections carry the standard at a different level by design (on
    one run the calibration injections' internal standards run about seven
    times the samples'), so including them would move the median toward a population the
    samples are not part of.

    A file whose standards all left the band together, in step, is a diluted or concentrated
    injection rather than a mis-integration (``IS_DILUTION_AGREEMENT``). Its standards are
    kept and the reading is noted. Only a standard that moved on its own is dropped.

    Dropping the area, rather than correcting it, is the whole point: a congener normalised
    by a standard that was never measured is not a quantity, and the pipeline reports it as
    unquantified instead of as a plausible wrong number.
    """
    notes: list[str] = []
    if len(reported_ids) < 3:
        return notes  # too few files for a median to mean anything
    low, high = IS_CONSISTENCY_RANGE
    medians: dict[str, float] = {}
    for label in sorted({label for fid in reported_ids for label in records[fid].is_areas}):
        values = [records[fid].is_areas[label] for fid in reported_ids if label in records[fid].is_areas]
        if len(values) >= 3:
            median = float(np.median(values))
            if median > 0:
                medians[label] = median
    for fid in reported_ids:
        ratios = {
            label: records[fid].is_areas[label] / median
            for label, median in medians.items()
            if label in records[fid].is_areas
        }
        outside = [label for label, ratio in ratios.items() if not low <= ratio <= high]
        if not outside:
            continue
        smallest, largest = min(ratios.values()), max(ratios.values())
        in_step = (
            len(ratios) >= 2
            and smallest > 0
            and largest / smallest <= IS_DILUTION_AGREEMENT
            and IS_DILUTION_RANGE[0] <= smallest
            and largest <= IS_DILUTION_RANGE[1]
        )
        if in_step:
            factor = float(np.exp(np.mean(np.log(list(ratios.values())))))
            notes.append(
                f"{fid}: every internal standard at {factor:.2f} times the run median, in step, "
                f"read as a diluted or concentrated injection. The standards are used as measured"
            )
            continue
        for label in outside:
            area = records[fid].is_areas.pop(label)
            notes.append(
                f"{fid}: internal standard {label} integrated to {area:,.0f} against a run "
                f"median of {medians[label]:,.0f}, so it is not used. Congeners normalised by it "
                f"are reported unquantified"
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
    record: FileRecord, analytes: list[Analyte], anchors: dict[str, float], params: IntegrationParams,
    unresolved: frozenset[tuple[str, str]] = frozenset(),
    ghost_scale: GhostScale | None = None,
) -> None:
    """Read, integrate and release one file's chromatograms, filling ``record`` in place.

    Chromatograms are dropped before returning: a run is ~30 files on a host with ~2 GB, and
    the peaks are all that is needed afterwards. ``integrate_run``'s warning about unplaced
    analytes is suppressed and recomputed as data on the record: the same fact, kept
    somewhere a report can read it. With a ``ghost_scale`` the file's ghost peaks are kept out
    of its integration (:func:`core.ghosts.integrate_with_ghosts`) and recorded on it.
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
            record.peaks, record.ghosts = integrate_with_ghosts(
                chroms, analytes, anchors, params, unresolved, ghost_scale
            )
        record.is_window_max = _standard_window_max(chroms, analytes, anchors, params)
    finally:
        del chroms
    placed = {p.label for p in record.peaks}
    record.skipped_analytes = [a.label for a in analytes if a.label not in placed]
    record.is_areas = internal_standard_areas(record.peaks, analytes)


#: How far either side of its expected retention time an internal standard's channel is read
#: to decide whether the standard is in the injection at all (minutes). Wider than the drift
#: any calibration corrects, narrower than the gap between the two standards.
IS_PRESENCE_WINDOW = 1.0

#: A channel whose tallest point in that window stays below this share of the standard's
#: typical peak height in the run's samples and blanks carries no standard. Measured case: the
#: one archive injection without an internal standard tops out at 152 counts on the 30D
#: channel and 223 on the 204 channel, against standard peaks of 9,000 to 14,000 in the files
#: around it.
IS_PRESENCE_FRACTION = 0.05

#: Sample-level flag for an internal standard a reported file has no usable area for, as
#: ``missing_internal_standard:<label>=<reason>``. ``reason`` is ``absent`` (the standard is not
#: in the injection), ``not_found`` (its channel carries a peak near where it was expected, but
#: no standard was detected there) or ``rejected`` (detected, and refused by
#: :func:`reject_inconsistent_standards`). Every congener the standard normalises is then
#: reported without a mass. A number scaled by a standard that was never measured would be a
#: guess, and no number is the honest result.
MISSING_IS_FLAG = "missing_internal_standard"


def _standard_window_max(
    chroms, analytes: list[Analyte], anchors: dict[str, float], params: IntegrationParams
) -> dict[str, float]:
    """Each internal standard's tallest smoothed point within :data:`IS_PRESENCE_WINDOW`."""
    out: dict[str, float] = {}
    for analyte in analytes:
        if analyte.role not in IS_ROLES or analyte.label not in anchors:
            continue
        chrom = get_transition(chroms, analyte.q1, analyte.q3)
        if chrom is None or chrom.times.size < 3:
            continue
        window = np.abs(chrom.times - anchors[analyte.label]) <= IS_PRESENCE_WINDOW
        if window.any():
            trace = smooth(chrom.intensities, params.smooth_window, params.smooth_poly)
            out[analyte.label] = float(np.max(trace[window]))
    return out


def classify_missing_standards(
    records: Mapping[str, "FileRecord"], reported_ids: Sequence[str], is_labels: Sequence[str]
) -> dict[str, dict[str, str]]:
    """File id -> {internal-standard label: reason} for every standard a reported file lacks.

    Reasons as in :data:`MISSING_IS_FLAG`. "Absent" is judged against the standard's own
    typical height in this run's reported files, so it needs no absolute threshold: a channel
    that stays below :data:`IS_PRESENCE_FRACTION` of it near the expected time has no standard.
    A file that decoded to nothing is reported as such elsewhere and is not classified here.
    """
    typical: dict[str, float] = {}
    for label in is_labels:
        heights = [
            float(peak.height)
            for fid in reported_ids
            for peak in records[fid].peaks
            if peak.label == label and peak.detected and peak.height > 0
        ]
        typical[label] = float(np.median(heights)) if heights else float("nan")
    out: dict[str, dict[str, str]] = {}
    for fid in reported_ids:
        record = records[fid]
        if not record.usable:
            continue
        for label in is_labels:
            if label in record.is_areas:
                continue
            peak = next((p for p in record.peaks if p.label == label), None)
            if peak is not None and peak.detected:
                reason = "rejected"
            else:
                height = record.is_window_max.get(label, 0.0)
                expected = typical.get(label, float("nan"))
                absent = not np.isfinite(expected) or height < IS_PRESENCE_FRACTION * expected
                reason = "absent" if absent else "not_found"
            out.setdefault(fid, {})[label] = reason
    return out


def _missing_standard_note(name: str, missing: dict[str, str], is_labels: Sequence[str],
                           expected_rt: Mapping[str, float]) -> str | None:
    """The run note for a file's missing standards. ``rejected`` ones are noted by the guard."""
    own = {label: why for label, why in missing.items() if why != "rejected"}
    if not own:
        return None
    parts = []
    for label, why in own.items():
        rt = expected_rt.get(label)
        where = f" at {rt:.2f} min" if rt is not None else ""
        if why == "absent":
            parts.append(f"internal standard {label} is not in this injection (no peak within "
                         f"{IS_PRESENCE_WINDOW:g} min of its expected time{where})")
        else:
            parts.append(f"internal standard {label} was not found at its expected time{where}, "
                         "although its channel carries a peak nearby. Check the retention times")
    scope = ("this file" if set(missing) >= set(is_labels)
             else "the congeners " + " and ".join(missing) + " normalises")
    return f"{name}: " + "; ".join(parts) + f". No mass is reported for {scope}"


#: How far either side of its expected time a lost labelled standard is looked for (minutes).
LANDMARK_SEARCH = 0.4
#: The peak found there must reach this share of the standard's typical height in the run ...
LANDMARK_FRACTION = 0.15
#: ... and stay below this multiple of it. A standard is spiked at one amount, so a peak several
#: times its usual height on that channel is something else eluting there.
LANDMARK_CEILING = 3.0
#: The offsets of the standards found must lie on one straight line, none further from it than
#: this (minutes) ...
LANDMARK_AGREEMENT = 0.1
#: ... and the line may rise or fall no faster than this (minutes of offset per minute of
#: retention time). Real drift is far gentler.
LANDMARK_SLOPE = 0.05
#: Peak flag for an analyte searched at a retention time moved by :func:`realign_shifted_file`.
REALIGNED_FLAG = "realigned"


def landmark_labels(analytes: Sequence[Analyte]) -> set[str]:
    """Labelled standards on channels that carry no native congener.

    On such a channel the tallest peak near the expected time can only be the standard, so its
    position is evidence of where this file's chromatography put that part of the run.
    """
    native = {(round(a.q1, 1), round(a.q3, 1)) for a in analytes if a.role == "native"}
    return {a.label for a in analytes
            if a.role != "native" and (round(a.q1, 1), round(a.q3, 1)) not in native}


def typical_heights(
    records: Mapping[str, "FileRecord"], ids: Sequence[str], labels: set[str]
) -> dict[str, float]:
    """Median height of each label's detected peaks across ``ids``."""
    out: dict[str, float] = {}
    for label in labels:
        heights = [float(p.height) for fid in ids for p in records[fid].peaks
                   if p.label == label and p.detected and p.height > 0]
        out[label] = float(np.median(heights)) if heights else float("nan")
    return out


def _internal_standard_misplaced(
    record: "FileRecord", analytes: Sequence[Analyte], anchors: Mapping[str, float],
    params: IntegrationParams, landmarks: set[str],
) -> bool:
    """Whether an internal standard is lost in a file whose other standards show a drift.

    Standards found where they were expected still say where the chromatography went. Three of
    them with a median offset beyond half the apex search window, in a file whose internal
    standard was not detected, is a file searched in the wrong place even though fewer than two
    standards were lost outright.
    """
    by_label = {p.label: p for p in record.peaks}
    lost_internal = [a.label for a in analytes
                     if a.role in IS_ROLES and a.label in anchors and a.label in by_label
                     and not by_label[a.label].detected]
    offsets = [float(by_label[label].rt_apex) - anchors[label] for label in landmarks
               if label in anchors and label in by_label and by_label[label].detected]
    return bool(lost_internal) and len(offsets) >= 3 and (
        abs(float(np.median(offsets))) > params.rt_window / 2)


def _lost_standard_offsets(
    chroms, analytes: Sequence[Analyte], anchors: Mapping[str, float], params: IntegrationParams,
    typical: Mapping[str, float], lost: Sequence[str],
) -> dict[str, float]:
    """Where each lost labelled standard went, as its offset from its anchor (minutes).

    The standard is the tallest peak within :data:`LANDMARK_SEARCH` of its anchor, if that peak
    reaches :data:`LANDMARK_FRACTION` of the standard's typical height and stays under
    :data:`LANDMARK_CEILING` times it. A standard with no such peak is left out.
    """
    analyte_by = {a.label: a for a in analytes}
    found: dict[str, float] = {}
    for label in lost:
        analyte = analyte_by[label]
        chrom = get_transition(chroms, analyte.q1, analyte.q3)
        if chrom is None or chrom.times.size < 5:
            continue
        trace = smooth(chrom.intensities, params.smooth_window, params.smooth_poly)
        near = np.abs(chrom.times - anchors[label]) <= LANDMARK_SEARCH
        maxima = [int(i) for i in local_maxima(trace) if near[int(i)]]
        if not maxima:
            continue
        top = max(maxima, key=lambda i: trace[i])
        if LANDMARK_FRACTION * typical[label] <= trace[top] <= LANDMARK_CEILING * typical[label]:
            found[label] = float(chrom.times[top]) - anchors[label]
    return found


def realign_shifted_file(
    record: "FileRecord", analytes: list[Analyte], anchors: dict[str, float],
    params: IntegrationParams, typical: Mapping[str, float], landmarks: set[str],
    unresolved: frozenset[tuple[str, str]] = frozenset(),
    ghost_scale: GhostScale | None = None,
) -> str | None:
    """Re-integrate a file whose chromatography drifted from its calibration's, where it went.

    Each file is integrated at the retention times of its nearest calibration injection. If
    part of one file's chromatogram elutes late, the analytes there are searched in the wrong
    place: reported absent, or credited with a neighbour's peak, internal standard included.
    Measured cases: in one archive file everything from 16 to 24 min eluted 0.2 min late; its
    internal standard 30D and three surrogates were lost, and with them the masses of a hundred
    congeners. In two concentrated archive runs the drift was a slope rather than a step, the
    standards sitting 0.05 min late at 23 min and 0.33 min late at 67 min, which one shift fits
    only roughly.

    The labelled standards are spiked into every file, so they show where this file's
    chromatography went. The correction acts when at least two of them (on channels with no
    native congener, :func:`landmark_labels`) were lost here, or when an internal standard was
    lost while three standards found where expected show an offset
    (:func:`_internal_standard_misplaced`). The lost standards must be found within
    :data:`LANDMARK_SEARCH` at :data:`LANDMARK_FRACTION` to :data:`LANDMARK_CEILING` times their
    usual height, and lie on one gentle line (:data:`LANDMARK_AGREEMENT`,
    :data:`LANDMARK_SLOPE`). The file's drift curve then runs through every standard it shows,
    those found where they went and those found where expected, and each analyte whose offset
    on that curve exceeds half the apex search window is searched where the curve puts it.

    The new integration is kept only if it recovers the standards that were lost, loses none
    that was found before, and finds *more* native peaks where it moved the search than before.
    A shift the native congeners do not share leaves their count where it was, so a file with no
    native peak to corroborate it is left as it was. Otherwise the file is left exactly as
    integrated. Every moved peak carries the ``realigned`` flag. Returns the run note, or
    ``None``. The re-integration keeps ghost peaks out exactly as the first did
    (``ghost_scale``, see :func:`_integrate_one`).
    """
    by_label = {p.label: p for p in record.peaks}
    lost = [label for label in landmarks
            if label in anchors and label in by_label and not by_label[label].detected
            and np.isfinite(typical.get(label, np.nan))]
    if not record.mzml.exists():
        return None
    if len(lost) < 2 and not _internal_standard_misplaced(record, analytes, anchors, params,
                                                           landmarks):
        return None
    chroms = read_srm_mzml(record.mzml)
    try:
        return _realign_with(chroms, record, analytes, anchors, params, typical, landmarks, lost,
                             unresolved, ghost_scale)
    finally:
        del chroms


def _realign_with(chroms, record, analytes, anchors, params, typical, landmarks, lost,
                  unresolved=frozenset(), ghost_scale=None):
    by_label = {p.label: p for p in record.peaks}
    internal = [a.label for a in analytes if a.role in IS_ROLES and a.label in by_label]
    found = _lost_standard_offsets(chroms, analytes, anchors, params, typical, lost)
    kept = {label: float(by_label[label].rt_apex) - anchors[label] for label in landmarks
            if label in anchors and label in by_label and by_label[label].detected}
    points = {**kept, **found}
    if len(points) < 3 or (len(lost) >= 2 and len(found) < 2):
        return None

    # The standards that must agree: the lost ones found, or, when fewer than two were lost
    # and an internal standard was, every standard that shows an offset.
    steady = found if len(found) >= 2 else points
    if abs(float(np.median(list(steady.values())))) <= params.rt_window / 2:
        return None
    at = np.array([anchors[label] for label in steady])
    off = np.array(list(steady.values()))
    slope, level = np.polyfit(at, off, 1) if np.ptp(at) > 0 else (0.0, float(off.mean()))
    if abs(slope) > LANDMARK_SLOPE or float(np.max(np.abs(off - (level + slope * at)))) > (
            LANDMARK_AGREEMENT):
        return None

    # The drift curve runs through every standard the file shows. An analyte is searched
    # where it puts it, once that is more than half the apex search window from the anchor.
    by_anchor = {anchors[label]: offset for label, offset in points.items()}
    xs = np.array(sorted(by_anchor))
    ys = np.array([by_anchor[x] for x in xs])
    drift = {label: float(np.interp(rt, xs, ys)) for label, rt in anchors.items()}
    moved_labels = {label for label, d in drift.items() if abs(d) >= params.rt_window / 2}
    if not moved_labels:
        return None
    moved = {label: rt + drift[label] if label in moved_labels else rt
             for label, rt in anchors.items()}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        peaks, ghosts = integrate_with_ghosts(chroms, analytes, moved, params, unresolved,
                                              ghost_scale)
    after = {p.label: p for p in peaks}

    standards = {label for label in landmarks if label in by_label} | set(internal)
    if not all(label in after and after[label].detected for label in found):
        return None
    if any(by_label[label].detected and not (label in after and after[label].detected)
           for label in standards):
        return None
    if len(lost) < 2 and not any(
            not by_label[label].detected and label in after and after[label].detected
            for label in internal):
        return None
    span = [a.label for a in analytes if a.role == "native" and a.label in moved_labels]
    before_n = sum(1 for label in span if label in by_label and by_label[label].detected)
    after_n = sum(1 for label in span if label in after and after[label].detected)
    if after_n <= before_n:
        return None
    for peak in peaks:
        if peak.label in moved_labels:
            peak.flags.append(REALIGNED_FLAG)
    record.peaks = peaks
    record.ghosts = ghosts
    record.is_areas = internal_standard_areas(peaks, analytes)
    record.is_window_max = _standard_window_max(chroms, analytes, moved, params)

    lo, hi = min(found.values(), default=0.0), max(found.values(), default=0.0)
    where = sorted(anchors[label] for label in moved_labels)
    biggest = max((drift[label] for label in moved_labels), key=abs)
    standards_note = (f"labelled standards {', '.join(sorted(found))} eluted "
                      + (f"{lo:+.2f} min" if hi - lo < 0.005 else f"{lo:+.2f} to {hi:+.2f} min")
                      + " from where the calibration puts them" if found else
                      "the labelled standards that were found show an offset from where the "
                      "calibration puts them")
    return (f"{record.entry.name}: {standards_note}, so the {len(moved_labels)} analytes from "
            f"{where[0]:.2f} to {where[-1]:.2f} min were searched where the drift through the "
            f"file's {len(points)} standards puts them, by up to {biggest:+.2f} min. "
            f"Native peaks found there: {before_n} before, {after_n} after")


#: How many blanks (or reference injections) the run's ghost-peak scale is the median of.
GHOST_SCALE_FILES = 4


def _measure_ghost_scale(
    quantified: Sequence[FileEntry], converted_dir: Path, analytes: list[Analyte],
    anchors_for_name: Callable[[str], dict[str, float]], params: IntegrationParams,
) -> GhostScale:
    """The run's internal-standard 204 ghost fraction (:func:`core.ghosts.run_scale`).

    Read in the run's blanks, where 204 stands nearly alone, each at the anchors of the
    calibration it is integrated against, until :data:`GHOST_SCALE_FILES` of them have measured
    both steps: the fraction is the instrument's, steady across a run, so more blanks only cost
    reads. The reference injections are read only for a step the blanks do not measure. Each
    file is read and released in turn, as in the integration loop.
    """
    measured: dict[str, list[dict[int, float]]] = {"blank": [], "ref": []}
    for kind in ("blank", "ref"):
        if kind == "ref" and all(any(step in f for f in measured["blank"]) for step in (1, 2)):
            break
        for entry in quantified:
            if entry.file_type != kind:
                continue
            if sum(1 for f in measured[kind] if 1 in f and 2 in f) >= GHOST_SCALE_FILES:
                break
            path = _mzml_path(entry, converted_dir)
            if not path.exists():
                continue
            chroms = read_srm_mzml(path)
            try:
                if not any(c.times.size > 1 for c in chroms.values()):
                    continue
                fractions = is_fraction(chroms, analytes, anchors_for_name(entry.name), params,
                                        min_r=SCALE_MIN_R)
            finally:
                del chroms
            if fractions:
                measured[kind].append(fractions)
    return run_scale(measured["blank"], measured["ref"])


def select_calibrations(
    fits: Mapping[str, "RunFitReport"], order: Sequence[str]
) -> tuple[list[str], dict[str, str]]:
    """Which calibration injections may anchor files and supply response factors.

    Every injection is asked the run-fit question. One whose retention-time table does not
    fit is a different chromatographic frame, or an injection that failed, and its anchors
    send every file integrated at them to the wrong retention times. Measured case: a run's
    first calibration sat 4.3 min late against its other three, and the four files nearest it
    in time, a sample, a reference injection and two blanks, were all integrated there and
    lost both internal standards. Such an injection is excluded, provided another fits. If
    none fits, all are kept and the run fails the fit check as it always did.

    Returns the usable names in ``order``, and ``{excluded name: reason}``.
    """
    usable = [name for name in order if name in fits and fits[name].fit_ok]
    if not usable:
        return [name for name in order if name in fits], {}
    excluded = {name: (fits[name].reason or "retention-time table does not fit")
                for name in order if name in fits and name not in usable}
    return usable, excluded


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
    fit_by_cal = {first_cal.name: fit}
    alignment_by_cal = {first_cal.name: alignment}
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
            extra_anchors = build_anchors(
                extra_chroms, extra_seeded, params, alignment=extra_alignment
            )
            anchors_by_cal[extra.name] = extra_anchors
            fit_by_cal[extra.name] = rt_fit_quality(extra_chroms, extra_seeded, extra_anchors, params)
            alignment_by_cal[extra.name] = extra_alignment
        finally:
            del extra_chroms

    # A calibration whose retention-time table does not fit anchors nothing and calibrates
    # nothing, as long as another one fits (see :func:`select_calibrations`). It is still
    # integrated, at its own anchors, so the Calibration view shows what it is.
    usable_cals, excluded_cals = select_calibrations(
        fit_by_cal, [e.name for e in cal_entries if e.name in anchors_by_cal]
    )
    own_anchors = {name: anchors_by_cal.pop(name) for name in excluded_cals}
    primary_cal = usable_cals[0] if usable_cals else first_cal.name
    alignment = alignment_by_cal.get(primary_cal, alignment)
    fit = fit_by_cal.get(primary_cal, fit)

    # Each calibration was anchored on its own, so each could independently slip a window
    # by one peak. Comparing them with each other catches that. See
    # :func:`core.calibrate.reconcile_anchors`. Every repair is recorded rather than applied
    # silently.
    anchors_by_cal, anchor_repairs = reconcile_anchors(
        anchors_by_cal, [e.name for e in cal_entries if e.name in anchors_by_cal]
    )
    anchors = anchors_by_cal.get(primary_cal, anchors)

    notes: list[str] = [
        f"{name}: {reason}. This calibration injection anchors no other file and contributes "
        "no response factor"
        for name, reason in excluded_cals.items()
    ]
    notes.extend(anchor_repairs)

    # Congeners this run's chromatography does not separate. Each usable calibration is
    # integrated once at its own anchors: the standard holds every congener, so one there with
    # no apex of its own, inside a neighbour's peak, eluted with that neighbour (see
    # :func:`core.integrate.unresolved_pairs`). Its anchor was only a prediction, so the
    # envelope is fitted as two components with free positions
    # (:func:`core.deconvolve.fit_unresolved_pair`), and the two fitted retention times become
    # both congeners' anchors for that calibration. Every file, the calibration included, then
    # splits the pair's shared envelope at that measured spacing, so the response factor and
    # the masses it scales are measured alike. A pair whose envelope the fit cannot describe
    # is left as before.
    by_name = {e.name: e for e in cal_entries}
    found: set[tuple[str, str]] = set()
    transition_of = {a.label: (a.q1, a.q3) for a in analytes}
    for name in list(anchors_by_cal):
        cal_anchors = anchors_by_cal[name]
        cal_path = _mzml_path(by_name[name], converted_dir) if name in by_name else None
        if cal_path is None or not cal_path.exists():
            continue
        chroms = read_srm_mzml(cal_path)
        try:
            if not any(c.times.size > 1 for c in chroms.values()):
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                first_pass = {p.label: p for p in integrate_run(chroms, analytes, cal_anchors, params)}
            measured = dict(cal_anchors)
            for one, two in sorted(unresolved_pairs(list(first_pass.values()), analytes, cal_anchors, params)):
                kept = first_pass[one] if first_pass[one].detected else first_pass[two]
                chrom = get_transition(chroms, *transition_of[one])
                if chrom is None:
                    continue
                pair_fit = fit_unresolved_pair(
                    chrom.times, chrom.intensities, cal_anchors[one], cal_anchors[two],
                    kept.i_start, kept.i_end,
                )
                if pair_fit is None:
                    continue
                measured[one], measured[two] = pair_fit.first.centre, pair_fit.second.centre
                found.add((one, two))
            anchors_by_cal[name] = measured
        finally:
            del chroms
    anchors = anchors_by_cal.get(primary_cal, anchors)
    unresolved = frozenset(found)
    notes.extend(
        f"congeners {one} and {two} are not resolved in this run's calibration: their retention "
        "times were measured by a two-peak fit there, and wherever they share one envelope it is "
        "split between them at that spacing (flag 'unresolved')"
        for one, two in sorted(unresolved)
    )

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
    for name in excluded_cals:
        calibration_for[name] = name

    def anchors_for_name(name: str) -> dict[str, float]:
        paired = calibration_for.get(name, "")
        return anchors_by_cal.get(paired) or own_anchors.get(paired) or anchors

    # Ghost peaks (see :mod:`core.ghosts`). Every congener's in-source fragments leave a copy of
    # its peak on lower homologs' channels, at a size the instrument sets. The scale is the
    # internal standard 204's copy, read where 204 stands nearly alone, in the run's blanks, and
    # it has to be known before the first file is integrated: one extra read per blank.
    ghost_scale: GhostScale | None = None
    if params.ghost_correction:
        ghost_scale = _measure_ghost_scale(quantified, converted_dir, analytes, anchors_for_name,
                                           params)
        notes.append(ghost_scale.note())

    anchors_for_file: dict[str, dict[str, float]] = {}
    for i, entry in enumerate(quantified):
        report("integrate", i, len(quantified), entry.sample_id)
        file_id = _unique_id(entry, set(records))
        record = FileRecord(entry=entry, mzml=_mzml_path(entry, converted_dir))
        file_anchors = anchors_for_name(entry.name)
        anchors_for_file[file_id] = file_anchors
        _integrate_one(record, analytes, file_anchors, params, unresolved, ghost_scale)
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

    # A reported file whose labelled standards show its chromatography drifted is searched where
    # the drift puts each analyte (see :func:`realign_shifted_file`). Every other file is left
    # exactly as integrated.
    landmarks = landmark_labels(analytes)
    typical = typical_heights(records, sample_ids + blank_ids, landmarks)
    for fid in sample_ids + blank_ids:
        if records[fid].usable:
            note = realign_shifted_file(
                records[fid], analytes, anchors_for_file[fid], params, typical, landmarks,
                unresolved, ghost_scale,
            )
            if note:
                notes.append(note)

    # A standard that was not detected never enters ``is_areas``, so nothing it normalises gets
    # a mass. Say so, file by file, and say whether the standard is missing from the injection
    # or was looked for in the wrong place: the two have different remedies.
    is_labels = [a.label for a in analytes if a.role in IS_ROLES]
    undetected = classify_missing_standards(records, sample_ids + blank_ids, is_labels)
    for fid, missing in undetected.items():
        record = records[fid]
        paired = calibration_for.get(record.entry.name, "")
        note = _missing_standard_note(
            record.entry.name, missing, is_labels,
            anchors_by_cal.get(paired) or own_anchors.get(paired) or anchors,
        )
        if note:
            notes.append(note)

    # An internal standard is spiked identically into every reported file, so one that
    # disagrees with its own run is a mis-integration, and it would scale every congener
    # that shares it. Rejected here, before any mass is computed.
    notes.extend(reject_inconsistent_standards(records, sample_ids + blank_ids))
    missing_standards = classify_missing_standards(records, sample_ids + blank_ids, is_labels)

    # ---- 3) Response factors, averaged over every usable calibration injection -----------
    report("quantify", 0, len(_QUANTIFY_STEPS), _QUANTIFY_STEPS[0])
    models = [
        build_calibration(records[cid].peaks, analytes, name=records[cid].entry.name)
        for cid in cal_ids
        if records[cid].usable and records[cid].entry.name not in excluded_cals
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
        + [f"{MISSING_IS_FLAG}:{label}={why}"
           for label, why in missing_standards.get(file_id, {}).items()]
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
        unresolved_pairs=tuple(sorted(unresolved)),
        excluded_calibrations=excluded_cals,
        ghost_scale=ghost_scale,
    )
