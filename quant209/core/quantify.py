"""Turn integrated peak areas into congener masses: RRFs, the mass formula, recovery.

The lab quantifies by internal standard. Every formula below was decoded from the
analyst's own quantification workbook and checked cell by cell against it (the ``cal``
sheet supplies the first, a sample sheet the second and third):

    RRF(analyte)   = (cal_mass / cal_area) / (cal_mass_IS / cal_area_IS)
    Mass(analyte)  = RRF(analyte) * Area_sample(analyte) * (IS_spike / IS_area_sample)
    recovery(homolog) = Mass_surrogate(sample) / Mass_surrogate(reference)
    Corrected      = Mass / recovery   if recovery < 1 else Mass

**The internal standard is chosen per analyte**, by chlorination: d-PCB30 ("30D") normalises
the low-chlorination homologs and PCB204 ("204") the high ones. The choice is recorded in
each analyte's ``is_norm`` field and resolved to a table row by
:func:`core.standards.internal_standard_label`.

**The IS area in the mass formula is *this sample's* own**, not the calibration injection's.
The internal standards are spiked just before injection, so the ratio
``IS_spike / IS_area_sample`` carries whatever the instrument's response did between the
calibration and this sample. The congener mass is corrected for it.
An IS peak that is missing, zero, or that the presence test rejected therefore makes the
sample unquantifiable rather than merely imprecise. See :func:`internal_standard_areas`.

**Recovery correction is one-sided.** EPA-1668 practice is to correct for extraction losses
and never to deflate an apparent gain, so a recovery at or above 1 leaves the mass alone.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

from .config import Analyte
from .peaks import PeakResult
from .standards import HOMOLOG_SURROGATE, IS_SPIKE_NG, internal_standard_label

#: Roles in the canonical table that are internal standards (as opposed to 13C surrogates,
#: which are spiked pre-extraction and measure recovery rather than normalising response).
IS_ROLES = frozenset({"IS_d30", "IS_204"})

#: Absolute floor (counts·s) below which a *native* calibration peak is treated as a failed
#: integration rather than a measurement. See :func:`build_calibration`.
CAL_AREA_MIN = 500.0

#: ...and the same floor expressed relative to the run: this fraction of the median positive
#: native calibration area. The larger of the two applies.
CAL_AREA_FRAC = 0.02


@dataclass(frozen=True, eq=False)
class CalibrationModel:
    """The response factors one run quantifies against, and the evidence behind them.

    ``rrf`` is the reported model: analyte label -> relative response factor. ``cal_area`` and
    ``is_response`` are what it was computed from, kept so any RRF in a report can be
    recomputed by hand from two numbers. ``invalid`` records every analyte that was *refused*
    an RRF and why. An empty entry there and an absent entry in ``rrf`` mean different
    things, and an analyst reading a blank column is entitled to the difference.

    ``eq=False``: identity equality, matching the other result containers in this package.
    """

    rrf: dict[str, float] = field(default_factory=dict)
    cal_area: dict[str, float] = field(default_factory=dict)
    #: Internal-standard label -> ``cal_mass / cal_area`` on the calibration injection.
    is_response: dict[str, float] = field(default_factory=dict)
    #: Analyte label -> human-readable reason it has no RRF.
    invalid: dict[str, str] = field(default_factory=dict)
    #: Names of the calibration injections that contributed, in the order they were averaged.
    injections: tuple[str, ...] = ()
    #: The native-area floor this model applied (see :func:`build_calibration`).
    area_floor: float = CAL_AREA_MIN

    def rrf_for(self, label: str) -> float | None:
        return self.rrf.get(label)


def _is_label(is_norm: str, analytes: Sequence[Analyte]) -> str | None:
    """The internal-standard analyte label for an ``is_norm`` tag, or None if absent.

    A canonical table that carries no row for one of the two internal standards is a broken
    table, but it is broken in a way the caller can report on (every analyte normalised by the
    missing standard simply gets no RRF) rather than one that should abort a whole run.
    """
    try:
        return internal_standard_label(is_norm, analytes)
    except ValueError:
        return None


def build_calibration(
    peaks: Iterable[PeakResult],
    analytes: Sequence[Analyte],
    name: str = "",
) -> CalibrationModel:
    """Derive one calibration injection's RRF table from its integrated peaks.

    ``peaks`` is the output of ``core.integrate.integrate_run`` on the calibration file.
    Every sample gets the same call, so a response factor cannot be measured under one
    integration rule and consumed under another.

    **The mis-integrated-cal-peak guard.** A congener's RRF is only trustworthy if its
    calibration peak was actually integrated. The calibration standard contains all 209
    congeners, so a near-zero *native* cal area is a failed integration, not a real absence.
    ``RRF = (cal_mass / cal_area) / is_response`` explodes as the area goes to zero. Any
    sample signal at all would then turn into an absurd mass. Two tests refuse such
    analytes:

    * An area below ``max(CAL_AREA_MIN, CAL_AREA_FRAC * median(positive native cal areas))``.
    * A peak the presence test rejected. A rejected peak's area is a bounded estimate of
      whatever occupies the analyte's RT territory, so it is not a response measurement at
      any magnitude. It is judged on :attr:`core.peaks.PeakResult.detected` and never on
      ``sn``, which absent peaks report large and meaningless (documented at
      ``core.peaks.PeakResult.sn``).

    Both tests apply to natives only. A labelled standard sitting low on a quiet channel is
    ordinary, and there is no "the standard contains it" argument to make about a channel
    that carries exactly one analyte.

    A refused analyte is reported nowhere rather than reported wrongly: it gets no RRF, so
    :func:`quantify_masses` produces no mass for it in any file of the run, and ``invalid``
    carries the reason for the report.
    """
    by_label = {p.label: p for p in peaks}
    cal_area = {label: float(p.area) for label, p in by_label.items()}

    # IS response: cal_mass / cal_area for each internal standard the table declares. Gated
    # on detection for the same reason natives are. An unresolved IS peak would rescale every
    # RRF in the run.
    is_response: dict[str, float] = {}
    for analyte in analytes:
        if analyte.role not in IS_ROLES:
            continue
        result = by_label.get(analyte.label)
        if result is None or not result.detected or result.area <= 0 or np.isnan(analyte.cal_mass):
            continue
        is_response[analyte.label] = float(analyte.cal_mass) / float(result.area)

    native_areas = [
        cal_area[a.label]
        for a in analytes
        if a.role == "native" and cal_area.get(a.label, 0.0) > 0
    ]
    floor = (
        max(CAL_AREA_MIN, CAL_AREA_FRAC * float(np.median(native_areas)))
        if native_areas
        else CAL_AREA_MIN
    )

    rrf: dict[str, float] = {}
    invalid: dict[str, str] = {}
    for analyte in analytes:
        result = by_label.get(analyte.label)
        is_label = _is_label(analyte.is_norm, analytes)
        response = is_response.get(is_label) if is_label else None
        if result is None:
            invalid[analyte.label] = "no calibration peak (transition or RT anchor missing)"
            continue
        if np.isnan(analyte.cal_mass):
            invalid[analyte.label] = "no certified mass in the calibration standard"
            continue
        if response is None:
            invalid[analyte.label] = f"no usable {analyte.is_norm} internal-standard cal response"
            continue
        if result.area <= 0:
            invalid[analyte.label] = "calibration peak integrated to zero area"
            continue
        if analyte.role == "native" and not result.detected:
            invalid[analyte.label] = "calibration peak not detected (presence test rejected it)"
            continue
        if analyte.role == "native" and result.area < floor:
            invalid[analyte.label] = (
                f"calibration area {result.area:.0f} below the {floor:.0f} counts*s floor"
            )
            continue
        rrf[analyte.label] = (float(analyte.cal_mass) / float(result.area)) / response

    return CalibrationModel(
        rrf=rrf,
        cal_area=cal_area,
        is_response=is_response,
        invalid=invalid,
        injections=(name,) if name else (),
        area_floor=floor,
    )


def average_calibrations(models: Sequence[CalibrationModel]) -> CalibrationModel:
    """Average replicate calibration injections' RRFs, per analyte.

    A run normally carries several calibration injections. Averaging their response factors
    limits the effect of any one injection's error. An analyte contributes wherever it
    produced a usable response. A congener mis-integrated in one injection is still
    quantified from the others, and is refused only when *no* injection produced a usable
    response.

    ``cal_area``, ``is_response`` and ``area_floor`` are the *first* injection's. They are
    the audit trail for a single file, and averaging them would produce numbers that
    reconstruct none of the injections actually acquired.
    """
    if not models:
        raise ValueError("average_calibrations needs at least one calibration model")
    if len(models) == 1:
        return models[0]

    labels: set[str] = set()
    for model in models:
        labels.update(model.rrf)
    rrf = {
        label: float(np.mean([m.rrf[label] for m in models if label in m.rrf]))
        for label in labels
    }
    invalid: dict[str, str] = {}
    for model in models:
        for label, reason in model.invalid.items():
            if label not in rrf:
                invalid.setdefault(label, reason)

    base = models[0]
    return CalibrationModel(
        rrf=rrf,
        cal_area=base.cal_area,
        is_response=base.is_response,
        invalid=invalid,
        injections=tuple(name for m in models for name in m.injections),
        area_floor=base.area_floor,
    )


def internal_standard_areas(
    peaks: Iterable[PeakResult], analytes: Sequence[Analyte]
) -> dict[str, float]:
    """This file's usable internal-standard areas: ``{IS label: area}``.

    "Usable" is a detection question, not a magnitude one. An IS peak the presence test
    rejected is not an IS measurement. Because the whole sample is scaled by
    ``IS_spike / IS_area``, a bad one does not degrade the result gracefully. It multiplies
    every congener in the sample by a number derived from noise. Undetected and zero-area
    standards are therefore omitted entirely, which propagates as "this sample could not be
    quantified" rather than as a plausible wrong answer.

    Never gated on ``sn``: an absent peak reports S/N far above ``min_sn`` (see
    ``core.peaks.PeakResult.sn``), so an S/N test here would accept exactly the peaks this
    one rejects.
    """
    is_labels = {a.label for a in analytes if a.role in IS_ROLES}
    return {
        p.label: float(p.area)
        for p in peaks
        if p.label in is_labels and p.detected and p.area > 0
    }


def is_factor(
    is_areas: Mapping[str, float], is_norm: str, analytes: Sequence[Analyte]
) -> float | None:
    """``IS_spike / IS_area`` for the internal-standard group an analyte belongs to.

    Returns None when that standard was not usable in this file, which is the signal that
    nothing normalised by it can be quantified here.
    """
    label = _is_label(is_norm, analytes)
    area = is_areas.get(label) if label else None
    if not area:
        return None
    return IS_SPIKE_NG / float(area)


def quantify_masses(
    peaks: Iterable[PeakResult],
    analytes: Sequence[Analyte],
    calibration: CalibrationModel,
    is_areas: Mapping[str, float],
) -> dict[str, float]:
    """``Mass = RRF * Area * (IS_spike / IS_area)`` for every analyte this file can support.

    An analyte is skipped, reported nowhere rather than reported as zero, when it has no
    response factor (refused at the calibration step) or when its internal standard was not
    usable in this file. Zero means "measured and below the limit". Absent means "not
    measurable". The two must not be conflated in a report.

    Detection is deliberately *not* applied here. The analyst integrates and quantifies a
    shoulder they chose to split off a dominant neighbour, and so does this chain. Each
    peak's :attr:`core.peaks.PeakResult.detected` travels alongside the mass
    (``RunResult.detected``), so the report, the GUI and the QA/QC flags can gate on it
    where that is the right call.
    """
    by_label = {a.label: a for a in analytes}
    factors: dict[str, float | None] = {}
    masses: dict[str, float] = {}
    for result in peaks:
        analyte = by_label.get(result.label)
        if analyte is None:
            continue
        rrf = calibration.rrf_for(result.label)
        if rrf is None:
            continue
        if analyte.is_norm not in factors:
            factors[analyte.is_norm] = is_factor(is_areas, analyte.is_norm, analytes)
        factor = factors[analyte.is_norm]
        if factor is None:
            continue
        masses[result.label] = rrf * float(result.area) * factor
    return masses


def compute_recoveries(
    sample_masses: Mapping[str, float],
    reference_masses: Mapping[str, float],
    reference_detected: Mapping[str, bool] | None = None,
) -> dict[str, float]:
    """Per-homolog surrogate recovery: this sample's surrogate mass over the reference's.

    The reference injection *defines* 100% recovery. It carries the same surrogate spike but
    skipped the extraction and cleanup this sample went through, so the ratio is exactly the
    fraction of each homolog that survived them. Keyed by homolog (via
    :data:`core.standards.HOMOLOG_SURROGATE`) because that is what the correction consumes and
    what the analyst's Summary sheet reports.

    A homolog whose surrogate is missing from either file, or zero in the reference, yields no
    entry at all rather than a zero or an infinity: there is no evidence about it.

    ``reference_detected`` extends that rule to the reference peak's own detection. The
    denominator is a *measurement*. If the reference injection's own surrogate peak was not
    detected, the integrator flagged it ``absent``, so the number is a bounded estimate of
    whatever sat at that retention time rather than a peak. Dividing by it then does not
    produce a poor recovery. It produces a meaningless one. Measured: one reference injection
    in a 2024 run had its 13C 180 anchor land 0.9 min off its peak, and every hepta recovery
    measured against it came out near 80 instead of near 1. Refusing the division leaves the
    homolog uncorrected and visibly missing, which an analyst can act on. Reporting 8000%
    recovery silently deflates every hepta mass in the batch by the same factor.
    """
    recoveries: dict[str, float] = {}
    for homolog, surrogate in HOMOLOG_SURROGATE.items():
        sample = sample_masses.get(surrogate)
        reference = reference_masses.get(surrogate)
        if sample is None or not reference:
            continue
        if reference_detected is not None and not reference_detected.get(surrogate, True):
            continue
        recoveries[homolog] = float(sample) / float(reference)
    return recoveries


def correct_for_recovery(
    masses: Mapping[str, float],
    recoveries: Mapping[str, float],
    analytes: Sequence[Analyte],
) -> dict[str, float]:
    """Divide each native mass by its homolog's surrogate recovery, when that recovery < 1.

    One-sided by design (EPA-1668): correcting for extraction losses is standard practice,
    deflating a mass because a surrogate over-recovered is not. An apparent gain is far more
    often an integration or spiking artefact than a real one.

    A recovery of exactly 1, above 1, at or below 0, or missing altogether leaves the mass
    untouched. The zero/negative case cannot arise from a real measurement (mass and
    reference mass are both non-negative), but dividing by it would emit an infinity or flip
    a sign into a report, so the guard stays. Labelled standards carry no homolog and are
    never corrected, which is why the internal standard still reads exactly its own spike
    amount after this step.
    """
    by_label = {a.label: a for a in analytes}
    corrected: dict[str, float] = {}
    for label, mass in masses.items():
        analyte = by_label.get(label)
        recovery = recoveries.get(analyte.homolog) if analyte is not None else None
        if recovery is not None and 0.0 < recovery < 1.0:
            corrected[label] = float(mass) / recovery
        else:
            corrected[label] = float(mass)
    return corrected
