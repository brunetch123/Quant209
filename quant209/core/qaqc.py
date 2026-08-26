"""QA/QC: the limit of quantification, below-LOQ censoring, and per-sample flags.

Three separable judgements, kept separate because a lab may change any one of them without
touching the others:

    LOQ(analyte)  = mean(blank masses) + k * SD(blank masses)      (k parameterised)
    reported mass = 0 where mass ≤ LOQ                             (censor, never subtract)
    flags         = surrogate recovery outside 50-150%, plus the count of analytes whose
                    integration the presence test rejected

**LOQ is defined by this run's own blanks**, not by a stored detection limit. Every analyte
carries some blank level (laboratory background, carryover, a trace of the solvent), and a
mass that does not stand above the spread of the run's own blanks is not evidence of the
congener being in the sample. Published practice differs on ``k``. A 99% one-sided
interval (mean + 3SD) is common for air matrices and a 95% interval for food matrices. That
is why it is a parameter and not a constant.

**Censoring rather than blank subtraction** mirrors the lab's practice: a below-LOQ analyte is
reported as 0, and an above-LOQ one is reported at its full measured mass. Subtracting the
blank mean instead would push some reported masses negative and would silently apply a
correction to every congener in the run on the strength of a handful of blank injections.

**Recoveries are reported and flagged, never used to reject a sample.** The analyst decides
what a 40% hexa recovery means for their dataset. This module's job is to make sure they see
it.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .peaks import PeakResult


@dataclass(frozen=True)
class QaqcParams:
    """The three QA/QC judgements' tunables, with the lab's air-work defaults."""

    #: LOQ = mean + ``k_sd`` * SD of the blank masses. 3 is roughly a one-sided 99% interval.
    k_sd: float = 3.0
    #: Whether at-or-below-LOQ masses are reported as 0. Off reports the raw mass and leaves
    #: the censoring decision to whoever reads the report.
    censor_below_loq: bool = True
    #: Surrogate recovery bounds outside which a sample is flagged (fractions, not percent).
    recovery_low: float = 0.50
    recovery_high: float = 1.50
    #: Fewest blanks before an LOQ is claimed at all. Two is the minimum at which "the spread
    #: of the blanks" means anything. See :func:`compute_loq`.
    min_blanks_for_loq: int = 2


def compute_loq(blanks: pd.DataFrame, params: QaqcParams | None = None) -> pd.Series:
    """Per-analyte ``mean + k*SD`` over the run's blank masses, analyte-indexed.

    ``blanks`` is a label-by-sample frame of the blanks' *recovery-corrected* masses, falling
    back to the uncorrected mass when no reference injection exists to correct against. The
    LOQ has to live in the same units as the numbers it censors.

    A blank that reported no mass for an analyte contributes **0**, not a gap, the same as
    ``dict.get(label, 0.0)`` would return. The analyte was looked for in that blank and
    nothing was quantifiable, which is evidence of a low blank level rather than an absence
    of evidence. Dropping the cell instead would raise the LOQ of every analyte that is clean
    in most blanks and dirty in one.

    Fewer than ``min_blanks_for_loq`` blanks yields an empty Series, and
    :func:`censor_below_loq` then censors nothing: with one blank there is no spread to
    measure, so "mean + 3SD" is that blank's own mass and censoring against it would delete
    every sample result at or below a single injection's background. This function's
    ``ddof=1 if size > 1 else 0`` branch leaves a reachable one-blank path. Lowering
    ``min_blanks_for_loq`` to 1 still reaches it here, and still gives a zero-width limit.
    Sample standard deviation (``ddof=1``) otherwise, since the blanks are a sample of the
    run's background and not its population.
    """
    params = params or QaqcParams()
    if blanks.shape[1] < params.min_blanks_for_loq or blanks.empty:
        return pd.Series(dtype=float)
    values = blanks.astype(float).fillna(0.0)
    ddof = 1 if values.shape[1] > 1 else 0
    return values.mean(axis=1) + params.k_sd * values.std(axis=1, ddof=ddof)


def censor_below_loq(
    masses: pd.DataFrame, loq: pd.Series, params: QaqcParams | None = None
) -> pd.DataFrame:
    """Return ``masses`` with every at-or-below-LOQ cell set to 0. Never mutates its input.

    The comparison is ≤: a mass exactly at the limit has not cleared it. An analyte with
    no LOQ (one absent from the blanks entirely) is left alone rather than censored against
    an assumed limit of zero, and a NaN (the analyte was not quantifiable in that file at
    all) stays NaN, since "not measured" and "measured below the limit" are different report
    entries.
    """
    params = params or QaqcParams()
    if not params.censor_below_loq or loq.empty:
        return masses.copy()
    # ``axis=0`` broadcasts the analyte-indexed limits down the rows. Without it, pandas
    # aligns a bare ``DataFrame op Series`` on the *columns* (sample ids here), which would
    # silently match nothing and censor nothing. A NaN on either side compares False, so an
    # unmeasured mass stays NaN and an analyte with no limit is left alone. These are exactly
    # the two cases that must not be censored.
    below = masses.le(loq.reindex(masses.index), axis=0)
    return masses.mask(below, 0.0)


def recovery_flags(
    recoveries: Mapping[str, float], params: QaqcParams | None = None
) -> list[str]:
    """Flags for surrogate recoveries outside the acceptable band, one per homolog."""
    params = params or QaqcParams()
    flags: list[str] = []
    for homolog, recovery in recoveries.items():
        if recovery is None or not np.isfinite(recovery):
            continue
        if recovery < params.recovery_low:
            flags.append(f"low_recovery:{homolog}={recovery:.0%}")
        elif recovery > params.recovery_high:
            flags.append(f"high_recovery:{homolog}={recovery:.0%}")
    return flags


def integration_flags(peaks: Iterable[PeakResult]) -> list[str]:
    """One flag counting the analytes in this file whose integration is low-confidence.

    Counted on :attr:`core.peaks.PeakResult.detected`, the presence test's verdict, and
    **never** on ``sn``. An S/N-based count would double-count the same defect through a
    measure that does not survive it. An absent analyte's apex may sit outside its own
    integrated region, so its S/N is computed off an extrapolated baseline and comes out
    large (observed up to ~4.9e4). Absent peaks therefore almost never carry ``low_sn``, and
    an S/N-based count silently misses exactly the peaks it was meant to find. See
    ``core.peaks.PeakResult.sn``.

    ``low_sn`` on a *detected* peak is real information, a genuine peak near the detection
    limit, but it is a grading of a peak that exists, not a doubt about whether it does, so
    it belongs in the per-analyte flags the report already carries rather than in this count.
    """
    undetected = sum(1 for peak in peaks if not peak.detected)
    return [f"{undetected}_undetected_analytes"] if undetected else []


def sample_flags(
    recoveries: Mapping[str, float],
    peaks: Iterable[PeakResult],
    params: QaqcParams | None = None,
) -> list[str]:
    """Every QA/QC flag for one file: recovery band breaches, then the integration count."""
    return recovery_flags(recoveries, params) + integration_flags(peaks)
