"""Benchmark the automated quantification chain against the analyst's own numbers.

The area benchmark establishes that the engine's peak *areas* agree with analyst
manual integration. This measures what happens to that agreement once the areas are
turned into congener masses and corrected for recovery, the two links that carry every
calibration, internal-standard and surrogate decision in the method:

    area  --RRF, IS drift-->  mass (ng)  --surrogate recovery-->  corrected mass

Both links are scored against the analyst independently, because they can fail in opposite
directions and cancel in the final number:

* **Masses.** ``RunResult.masses`` (uncorrected, i.e. ``RRF * area * (IS_spike / IS_area)``),
  joined to the workbook's "Mass (ng)" column (``report.truth`` column J), which is the
  analyst's own pre-recovery-correction mass. Scored where the *analyst's* mass clears
  ``--min-mass``. Below that the analyst is quantifying an integration of noise, and a
  relative difference between two such numbers says nothing about the chain.
* **Recoveries.** ``RunResult.recoveries`` joined to the workbook Summary sheet's Recoveries
  block, per homolog per sample. Scored as absolute error in recovery *fraction* rather than
  as a relative difference: a recovery is already a ratio, and 0.02 on a recovery of 0.8 and
  0.02 on a recovery of 0.2 are the same measurement error, which a relative statistic would
  report as 2.5% and 10%.

Two further sections score the files those two links are *made of*, which the workbook holds
no sample sheet for and which nothing above therefore reaches:

* **Response factors**, against the analyst's own ``cal`` sheets. Every mass in the run is
  ``RRF * area * (IS_spike / IS_area)`` and the RRF comes from the calibration injection (the
  most crowded chromatogram in the sequence). An RRF error multiplies every sample, and it is
  invisible in the mass tables whenever it cancels against an area error of the opposite sign,
  which is exactly what it did on a crowded calibration window before the cluster baseline
  (see ``docs/METHODS.md`` section 2.4).
* **Reference injections**, against the analyst's ``ref`` sheets: the surrogate masses that are
  the denominator of every recovery in the batch.

Neither is gated. They are diagnostics that say *where* a mass disagreement came from, and both
are reported next to the analyst's own workbook-to-workbook reproducibility, below which a
disagreement is no longer a statement about the engine's integration.

The whole run goes through ``core.pipeline.run_analysis`` (the production entry point, with
production defaults) rather than through a bench-local reassembly of the same steps. A
benchmark that builds its own pipeline measures the benchmark.

Retention times come from the **aligned** path (``core.calibrate.align_rt_map``), which is
``run_analysis``'s only behaviour and ``bench_areas.py``'s default. There is no ``@canonical``
or ``@workbook`` switch here, because the point of this benchmark is to measure the pipeline as
it ships. See ``tools/bench_areas.py`` for those diagnostics on the area layer.

Gates (native analytes, detectable subset)
-------------------------------------------------------------------
    median |RPD| of uncorrected mass vs analyst  ≤  10%
    mean absolute error of recovery vs analyst   ≤  0.05

A run acquired under a different GC method is reported as a second table on the same
aligned path, with no tuning of
any kind: it is a *different GC method*, and it is here to show what the pipeline does with a
run the canonical retention-time table does not describe, not to be optimised against.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from dataclasses import dataclass, fields, replace
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

# Allow ``python tools/bench_quant.py`` from a source checkout without installing.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quant209.core.calibrate import align_rt_map, build_anchors  # noqa: E402
from quant209.core.config import load_canonical_peaks  # noqa: E402
from quant209.core.integrate import integrate_run  # noqa: E402
from quant209.core.mzml import read_srm_mzml  # noqa: E402
from quant209.core.peaks import IntegrationParams  # noqa: E402
from quant209.core.pipeline import RunResult, run_analysis  # noqa: E402
from quant209.core.qaqc import QaqcParams  # noqa: E402
from quant209.core.quantify import build_calibration  # noqa: E402
from quant209.core.sequence import find_data_files, scan_inputs  # noqa: E402
from quant209.core.standards import HOMOLOG_SURROGATE  # noqa: E402
from quant209.report.truth import read_run_truth  # noqa: E402
from quant209.report.validation import (  # noqa: E402
    BAND,
    GATE_MEDIAN_RPD,
    GATE_RECOVERY_MAE,
    HOMOLOG_ORDER,
    NEAR_LOQ_MULTIPLE,
    InjectionMatch,
    Metrics as MassMetrics,
    RecoveryMetrics,
    RecoveryRow,
    RefRow,
    RrfMetrics,
    RrfRow,
    MassRow,
    analyst_rrf_spread as _analyst_rrf_spread,
    classify_mass_miss as classify_miss,
    join_masses as _join_masses,
    join_recoveries as _join_recoveries,
    join_refs as _join_refs,
    join_rrfs as _join_rrfs,
    match_areas as _match_areas,
    match_injections as _match_injections,
    median_abs_rpd as _median_abs_rpd,
    metrics as mass_metrics,
    per_injection_rrf as _per_injection_rrf,
    recovery_metrics,
    rrf_metrics,
)

DEFAULT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = Path("docs/validation/quant.md")

def discover_runs(root: Path) -> list[str]:
    """Run-folder names under ``root``: any directory holding ``Converted`` or ``.D`` files.

    Used when no run is named on the command line. The archive is not distributed with the
    source and its folder names are not published, so there is no useful hardcoded default.
    """
    found: list[str] = []
    for d in sorted(root.iterdir()):
        if not d.is_dir() or d.name.startswith((".", "_")):
            continue
        if (d / "Converted").is_dir() or any(
            p.is_dir() and p.suffix.upper() == ".D" for p in d.iterdir()
        ):
            found.append(d.name)
    return found


#: Analyst mass (ng) below which a row is treated as noise rather than a measurement.
DEFAULT_MIN_MASS = 1.0

# BAND, GATE_MEDIAN_RPD, GATE_RECOVERY_MAE, NEAR_LOQ_MULTIPLE and HOMOLOG_ORDER live in
# ``quant209.report.validation`` (imported above), identical values shared with every
# other benchmark so a threshold changed in one place changes everywhere it is used.


# --------------------------------------------------------------------------------------
# Joined rows
# --------------------------------------------------------------------------------------
#
# MassRow, RecoveryRow, RrfRow, RefRow and InjectionMatch live in
# ``quant209.report.validation`` (imported above), unchanged, and shared with the
# per-run report generator. Only ``RunReport`` (this tool's own rendering aggregate)
# stays here.


@dataclass
class RunReport:
    """Everything measured for one analysis run.

    Deliberately holds only the joined rows and a handful of scalars, never the
    ``RunResult``: two runs are benchmarked in one process on a host with ~2 GB, and a
    ``RunResult`` carries every ``PeakResult`` of every file in the run.
    """

    run: str
    cal_files: list[str]
    ref_files: list[str]
    reference_policy: str
    #: Reference-injection file id -> the quantified files bracketed by it. The single
    #: orchestration decision that most moves a reported recovery, so it is printed.
    reference_batches: dict[str, list[str]]
    status: str
    status_reason: str
    notes: list[str]
    fit_ok: bool
    fit_summary: str
    rt_identity: bool
    rt_provenance: dict[str, int]
    aliases: dict[str, str]
    mass_rows: list[MassRow]
    recovery_rows: list[RecoveryRow]
    #: Response-factor agreement against the analyst's own ``cal`` sheets, natives only.
    rrf_rows: list[RrfRow]
    #: Which calibration injection each analyst workbook's ``cal`` sheet holds.
    injection_matches: list[InjectionMatch]
    #: The engine's RRF from one injection vs the analyst's from the workbook that used it:
    #: (workbook, cal id, n, median \|RPD\|, count outside ±25%).
    per_injection_rrf: list[tuple[str, str, int, float, int]]
    #: Median \|RPD\| between the analyst's *own* workbooks' RRFs for the same analyte, the
    #: reproducibility floor any agreement claim has to be read against.
    analyst_rrf_spread: float
    #: Reference-injection surrogate-mass agreement against the analyst's ``ref`` sheets.
    ref_rows: list[RefRow]
    n_files: int
    n_samples: int
    n_blanks: int
    n_natives: int
    #: Calibration injections that actually *contributed* a response factor, which is not the
    #: same as the number of cal files: one that decoded to nothing is listed but averaged over
    #: nothing (``CalibrationModel.injections``).
    n_cal_injections: int
    n_invalid_rrf: int
    loq_analytes: int
    censored_cells: int
    truth_only: list[str]
    file_only: list[str]
    recovery_unmatched: list[str]
    skipped_files: list[str]


# --------------------------------------------------------------------------------------
# Measurement
# --------------------------------------------------------------------------------------
#
# _join_masses, _join_recoveries, _median_abs_rpd, _match_injections, _match_areas,
# _analyst_rrf_spread, _join_rrfs, _per_injection_rrf and _join_refs live in
# ``quant209.report.validation`` (imported above, under these same names), unchanged.
# This tool's own contribution starts at ``sweep_cal_areas`` below, which uses
# ``_match_areas`` to identify a workbook/injection pairing without building a ``RunResult``.


def sweep_cal_areas(
    folder: Path, field: str, values: list[float], params: IntegrationParams
) -> list[str]:
    """Re-integrate this run's calibration injections at each value of an ``IntegrationParams``
    field and score their areas against the analyst's ``cal`` sheets.

    The threshold-selection evidence for ``overlap_frac``, generated rather than hand-run,
    so the result can be reproduced from source instead of taken on trust.

    Only the calibration files are re-integrated: they are what a cal-area comparison needs,
    and re-running the whole pipeline once per swept value would read every file of the run
    that many times for numbers no sample contributes to. Retention-time anchors come from the
    first calibration injection at the *baseline* parameters and are held fixed across the
    sweep, so the only thing varying between rows is the field being swept.
    """
    entries = scan_inputs(find_data_files(folder))
    cal_entries = [e for e in entries if e.file_type == "cal"]
    if not cal_entries:
        raise SystemExit(f"{folder.name}: no calibration injection to sweep")
    converted = folder / "Converted"
    analytes = load_canonical_peaks()
    natives = {a.label for a in analytes if a.role == "native"}

    first = converted / f"{cal_entries[0].path.stem}.mzML"
    if not first.exists():
        raise SystemExit(f"{folder.name}: {first.name} not found")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        chroms = read_srm_mzml(first)
        try:
            anchors = build_anchors(
                chroms, analytes, params, alignment=align_rt_map(chroms, analytes, params)
            )
        finally:
            del chroms

    # {cal id: {swept value: {label: area}}}, one file read per calibration injection.
    areas: dict[str, dict[float, pd.Series]] = {}
    for entry in cal_entries:
        path = converted / f"{entry.path.stem}.mzML"
        if not path.exists():
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            chroms = read_srm_mzml(path)
            try:
                per_value: dict[float, pd.Series] = {}
                for value in values:
                    swept = replace(params, **{field: value})
                    peaks = integrate_run(chroms, analytes, anchors, swept)
                    per_value[value] = pd.Series(
                        {p.label: float(p.area) for p in peaks}, dtype=float
                    )
                areas[entry.sample_id] = per_value
            finally:
                del chroms
    if not areas:
        raise SystemExit(f"{folder.name}: no calibration injection decoded")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        truth = read_run_truth(folder)

    # Pair each workbook with an injection at the *baseline* value, so every row of the table
    # scores the same pairing and the rows are comparable down the column.
    baseline = values[0]
    known = set(next(iter(areas.values()))[baseline].index)
    matches = _match_areas(
        truth.cal_areas, {cid: v[baseline] for cid, v in areas.items()}, natives, known
    )
    if not matches:
        raise SystemExit(f"{folder.name}: no workbook cal sheet could be matched to an injection")

    lines = [f"### {folder.name}", "", f"Pairing (at `{field}={baseline:g}`): "
             + ", ".join(f"`{m.workbook}` -> `{m.cal_id}`" for m in matches), ""]
    header = [f"`{field}`"] + [f"`{m.cal_id}` (vs `{m.workbook}`)" for m in matches]
    body = []
    for value in values:
        cells = [f"{value:g}"]
        for match in matches:
            column = truth.cal_areas[match.workbook]
            ours = areas[match.cal_id][value]
            labels = [
                label
                for label in column.index
                if label in natives
                and label in ours.index
                and np.isfinite(column.get(label, np.nan))
                and float(column[label]) > 0
            ]
            a = ours.loc[labels].to_numpy(dtype=float)
            b = column.loc[labels].to_numpy(dtype=float)
            mean = (a + b) / 2.0
            rpds = np.abs(np.divide(a - b, mean, out=np.zeros_like(mean), where=mean > 0))
            cells.append(f"{100 * float(np.median(rpds)):.2f}% / {int((rpds > BAND).sum())}")
        body.append(cells)
    lines += _table(header, body)
    lines += [
        "",
        "Each cell is *median \\|RPD\\| of native calibration area vs the analyst's `cal` "
        "sheet* / *count outside +/-25%*. A value above 1 disables clustering entirely (a "
        "valley cannot stand higher than the lower of the two apexes it lies between), so that "
        "row is the pre-cluster engine.",
        "",
    ]
    return lines


def bench_run(
    folder: Path,
    params: IntegrationParams,
    qaqc: QaqcParams,
    aliases: dict[str, str] | None = None,
) -> RunReport:
    """Run the production pipeline over ``folder`` and join its output to analyst truth.

    ``aliases`` maps a data file's id to the truth column it should be scored against.
    DIAGNOSTIC ONLY, for the case where the *ground truth* is provably attached to the wrong
    file: it lets the
    corrected pairing be measured and reproduced without editing anyone's workbook. It must
    never be used to shop for a better-agreeing column.
    """
    aliases = aliases or {}
    entries = scan_inputs(find_data_files(folder))
    if not entries:
        raise SystemExit(f"{folder.name}: no data files found")

    result = run_analysis(entries, folder / "Converted", params=params, qaqc=qaqc)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # unreadable-workbook warnings are surfaced separately
        truth = read_run_truth(folder)

    mass_rows, truth_only, file_only = _join_masses(result, truth.masses, aliases)
    recovery_rows, recovery_unmatched = _join_recoveries(result, truth.recoveries, aliases)

    # The calibration and reference injections themselves, against the analyst's own ``cal``
    # and ``ref`` sheets. These are not scored by the mass tables above (the workbook holds
    # no *sample* sheet for them), yet every mass in the run is computed from them, so an
    # error here reaches every number in the report while showing up nowhere else in it.
    # ``aliases`` is deliberately not applied: it re-pairs *sample* columns, and these sheets
    # are joined by injection, not by sample.
    natives = {a.label for a in result.analytes if a.role == "native"}
    injection_matches = _match_injections(result, truth.cal_areas)
    rrf_rows = _join_rrfs(result, truth.cal_rrfs, truth.cal_areas, injection_matches)
    per_injection = _per_injection_rrf(result, truth.cal_rrfs, injection_matches)
    analyst_spread = _analyst_rrf_spread(truth.cal_rrfs, natives)
    ref_rows = _join_refs(result, truth.ref_surrogate_masses, truth.workbook_samples)

    provenance: dict[str, int] = {}
    for analyte in result.analytes:
        key = result.alignment.source_of(analyte)
        provenance[key] = provenance.get(key, 0) + 1

    natives = [a for a in result.analytes if a.role == "native"]
    censored = int(
        (
            (result.masses_censored == 0.0)
            & (result.masses_corrected.fillna(0.0) != 0.0)
        )
        .to_numpy()
        .sum()
    )
    batches: dict[str, list[str]] = {}
    for file_id, ref_id in sorted(result.reference_for.items()):
        batches.setdefault(ref_id, []).append(file_id)

    return RunReport(
        run=folder.name,
        cal_files=list(result.cal_files),
        ref_files=list(result.ref_files),
        reference_policy=result.reference_policy,
        reference_batches=batches,
        status=result.status,
        status_reason=result.status_reason,
        notes=list(result.notes),
        fit_ok=result.fit.fit_ok,
        fit_summary=result.fit.summary,
        rt_identity=result.alignment.is_identity,
        rt_provenance=provenance,
        aliases=dict(aliases),
        mass_rows=mass_rows,
        recovery_rows=recovery_rows,
        rrf_rows=rrf_rows,
        injection_matches=injection_matches,
        per_injection_rrf=per_injection,
        analyst_rrf_spread=analyst_spread,
        ref_rows=ref_rows,
        n_files=len(result.quantified_entries),
        n_samples=len(result.sample_ids),
        n_blanks=len(result.blank_ids),
        n_natives=len(natives),
        n_cal_injections=len(result.calibration.injections),
        n_invalid_rrf=sum(1 for a in natives if a.label not in result.calibration.rrf),
        loq_analytes=int(result.loq.notna().sum()),
        censored_cells=censored,
        truth_only=truth_only,
        file_only=file_only,
        recovery_unmatched=recovery_unmatched,
        skipped_files=[
            fid for fid, rec in result.files.items() if not rec.usable
        ],
    )


# --------------------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------------------
#
# MassMetrics, RecoveryMetrics, mass_metrics, RrfMetrics, rrf_metrics, recovery_metrics and
# classify_miss live in ``quant209.report.validation`` (imported above as ``Metrics as
# MassMetrics``, ``metrics as mass_metrics`` and ``classify_mass_miss as classify_miss``),
# shared unchanged with the area benchmark and the per-run report.


def _split(rows: list[MassRow]) -> tuple[list[MassRow], list[MassRow]]:
    """(native rows, labelled-standard rows)."""
    return (
        [r for r in rows if r.role == "native"],
        [r for r in rows if r.role != "native"],
    )


def gate(report: RunReport, min_mass: float) -> tuple[bool, MassMetrics, RecoveryMetrics]:
    """The mass-agreement verdict for one run, and the two metrics it is taken on."""
    m = mass_metrics(_split(report.mass_rows)[0], min_mass)
    r = recovery_metrics(report.recovery_rows)
    passed = (
        np.isfinite(m.median_abs_rpd)
        and m.median_abs_rpd <= GATE_MEDIAN_RPD
        and np.isfinite(r.mae)
        and r.mae <= GATE_RECOVERY_MAE
    )
    return bool(passed), m, r


# --------------------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------------------


def _pct(x: float) -> str:
    return "n/a" if not np.isfinite(x) else f"{100 * x:.1f}%"


def _num(x: float, places: int = 4) -> str:
    return "n/a" if not np.isfinite(x) else f"{x:.{places}f}"


def _table(header: list[str], rows: list[list[str]]) -> list[str]:
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return out


def _mass_metrics_table(m: MassMetrics) -> list[str]:
    return _table(
        ["metric", "value"],
        [
            ["joined rows", f"{m.n_joined}"],
            ["detectable rows (analyst mass >= threshold)", f"{m.n_detectable}"],
            ["median \\|RPD\\|", _pct(m.median_abs_rpd)],
            ["within +/-25%", _pct(m.within_band)],
            ["mass-weighted total ratio (sum auto / sum manual)", _num(m.total_ratio)],
            ["median signed RPD (bias)", _pct(m.median_rpd)],
            ["detectable rows with no computed mass", f"{m.n_unquantified}"],
        ],
    )


def _recovery_metrics_table(r: RecoveryMetrics) -> list[str]:
    return _table(
        ["metric", "value"],
        [
            ["joined (sample, homolog) pairs", f"{r.n}"],
            ["mean absolute error (recovery fraction)", _num(r.mae)],
            ["median signed error (bias)", _num(r.median_error, 4)],
            ["max absolute error", _num(r.max_abs_error)],
            ["within +/-0.05", _pct(r.within_005)],
        ],
    )


def render_rrf(report: RunReport) -> list[str]:
    """The response-factor section: the calibration injections, scored against the analyst.

    Every mass in this run is ``RRF * area * (IS_spike / IS_area)``, and the RRF comes from the
    calibration injection (the most crowded chromatogram in the sequence, and the one file
    the analyst's workbook holds no sample sheet for). An error here multiplies every sample,
    and it is invisible in the mass tables above whenever it happens to cancel against an area
    error of the opposite sign.
    """
    lines: list[str] = ["", "### Response factors vs the analyst's `cal` sheets", ""]
    if not report.rrf_rows:
        return lines + ["No workbook in this run holds a populated `cal` sheet."]

    m = rrf_metrics(report.rrf_rows)
    lines += [
        "The analyst side is the **median across this run's workbooks**: ours is averaged over "
        "every calibration injection, and each workbook's `cal` sheet is one injection. The "
        "per-injection table below compares single injections like for like.",
        "",
    ]
    lines += _table(
        ["metric", "value"],
        [
            ["native analytes compared", f"{m.n}"],
            ["median \\|RPD\\| of RRF", _pct(m.median_abs_rpd)],
            ["within +/-25%", _pct(m.within_band)],
            [
                "outside +/-25%",
                f"**{m.n_off}** ({m.n_high} read high, {m.n_low} read low)",
            ],
            [
                "median \\|RPD\\| of the calibration areas behind them",
                _pct(m.median_abs_area_rpd),
            ],
            [
                "analyst's own workbook-to-workbook RRF spread (reproducibility floor)",
                _pct(report.analyst_rrf_spread),
            ],
        ],
    )
    lines += [
        "",
        "An RRF reading *high* is a calibration peak that integrated **short**: "
        "`RRF = (cal_mass / cal_area) / is_response` is inversely proportional to the area. "
        "The area column is therefore the diagnostic, and the analyst's own spread is the "
        "floor below which a disagreement stops being about our integration.",
    ]

    if report.injection_matches:
        lines += ["", "#### Which calibration injection each workbook used", ""]
        lines += [
            "The `cal` sheet does not name its injection (its \"Source:\" cell names the "
            "workbook it was copied from), so the injection is identified by scoring this "
            "run's own calibration files against the sheet's native areas.",
            "",
        ]
        lines += _table(
            ["workbook", "matched injection", "median \\|RPD\\| of areas", "runner-up"],
            [
                [
                    f"`{match.workbook}`", f"`{match.cal_id}`",
                    _pct(match.score), _pct(match.runner_up),
                ]
                for match in report.injection_matches
            ],
        )
        # A match whose runner-up ties it at printed precision has not actually identified the
        # injection. Said here rather than left for a reader to notice two equal numbers.
        ties = [
            m
            for m in report.injection_matches
            if np.isfinite(m.runner_up) and _pct(m.score) == _pct(m.runner_up)
        ]
        if ties:
            lines += [
                "",
                "**Not discriminating** for "
                + ", ".join(f"`{m.workbook}`" for m in ties)
                + ". The runner-up injection scores the same at this precision, so which "
                "injection that sheet holds is not resolved. No conclusion below depends on it. "
                "Both candidates are calibration injections of this run, and the "
                "per-injection RRF comparison is like-for-like either way, but the match "
                "should not be quoted as if it were established.",
            ]

    if report.per_injection_rrf:
        lines += ["", "#### Per-injection RRF agreement (no averaging on either side)", ""]
        lines += _table(
            ["workbook", "injection", "n", "median \\|RPD\\|", "outside +/-25%"],
            [
                [f"`{book}`", f"`{cal}`", f"{n}", _pct(med), f"{off}"]
                for book, cal, n, med, off in report.per_injection_rrf
            ],
        )

    lines += ["", "#### Worst 20 response factors (largest \\|RPD\\|)", ""]
    worst = sorted(report.rrf_rows, key=lambda r: abs(r.rpd), reverse=True)[:20]
    lines += _table(
        ["#", "analyte", "homolog", "our RRF", "analyst RRF", "ratio", "our cal area",
         "analyst cal area", "area ratio", "flags", "method"],
        [
            [
                f"{i}", r.label, r.homolog, f"{r.ours:.4f}", f"{r.analyst:.4f}",
                f"{r.ours / r.analyst:.2f}" if r.analyst else "n/a",
                f"{r.our_cal_area:.0f}", f"{r.analyst_cal_area:.0f}",
                f"{r.our_cal_area / r.analyst_cal_area:.2f}" if r.analyst_cal_area else "n/a",
                r.flags, r.method,
            ]
            for i, r in enumerate(worst, start=1)
        ],
    )

    grouped: dict[str, list[RrfRow]] = {}
    for r in report.rrf_rows:
        grouped.setdefault(r.homolog, []).append(r)
    lines += ["", "#### RRF agreement by homolog", ""]
    lines += _table(
        ["homolog", "n", "median \\|RPD\\|", "within +/-25%", "outside +/-25%"],
        [
            [
                homolog, f"{len(rs)}",
                _pct(float(np.median(np.abs([r.rpd for r in rs])))),
                _pct(float(np.mean(np.abs([r.rpd for r in rs]) <= BAND))),
                f"{int((np.abs([r.rpd for r in rs]) > BAND).sum())}",
            ]
            for homolog, rs in (
                (h, grouped[h]) for h in HOMOLOG_ORDER if h in grouped
            )
        ],
    )
    return lines


def render_refs(report: RunReport) -> list[str]:
    """The reference-injection section: the denominator of every recovery in the run."""
    lines: list[str] = ["", "### Reference injections vs the analyst's `ref` sheets", ""]
    if not report.ref_rows:
        return lines + [
            "No workbook in this run could be paired with a reference injection "
            "(see `_join_refs`: pairing requires every sample in a workbook to bracket to the "
            "same reference)."
        ]
    lines += [
        "Each workbook's `ref` sheet holds the 13C-surrogate masses that define 100% recovery "
        "for its samples. The sheet does not name its injection, so the pairing is derived "
        "from this run's own reference bracketing rather than assumed. These masses are the "
        "denominator of every recovery in the run: an error here moves every corrected mass "
        "in the batch by the same factor.",
        "",
    ]
    rpds = np.abs([r.rpd for r in report.ref_rows])
    lines += _table(
        ["metric", "value"],
        [
            ["(reference, homolog) pairs compared", f"{len(report.ref_rows)}"],
            ["median \\|RPD\\| of surrogate mass", _pct(float(np.median(rpds)))],
            ["within +/-25%", _pct(float(np.mean(rpds <= BAND)))],
            ["max \\|RPD\\|", _pct(float(np.max(rpds)))],
        ],
    )

    refs = sorted({r.ref_id for r in report.ref_rows})
    lines += ["", "#### Surrogate mass RPD by homolog and reference injection", ""]
    body = []
    for homolog in HOMOLOG_ORDER:
        cells = []
        for ref in refs:
            match = [r for r in report.ref_rows if r.ref_id == ref and r.homolog == homolog]
            cells.append(f"{100 * match[0].rpd:+.1f}%" if match else "-")
        if any(c != "-" for c in cells):
            body.append([homolog, *cells])
    lines += _table(["homolog", *[f"`{r}`" for r in refs]], body)
    pairs = ", ".join(
        f"`{r}` <- `{w}`"
        for r, w in sorted({(row.ref_id, row.workbook) for row in report.ref_rows})
    )
    lines += ["", f"Pairing used: {pairs}."]
    return lines


def render_run(report: RunReport, min_mass: float) -> list[str]:
    native, labeled = _split(report.mass_rows)
    m_native = mass_metrics(native, min_mass)
    m_labeled = mass_metrics(labeled, min_mass)
    rec = recovery_metrics(report.recovery_rows)
    passed, _, _ = gate(report, min_mass)

    title = report.run + (" (diagnostic: truth columns re-paired)" if report.aliases else "")
    lines: list[str] = [f"## {title}", ""]
    lines += [
        f"- Files integrated: **{report.n_files}** "
        f"({report.n_samples} sample, {report.n_blanks} blank, "
        f"{len(report.cal_files)} cal, {len(report.ref_files)} ref; hexane blanks excluded)",
        f"- RT anchors + alignment from: `{report.cal_files[0]}` "
        f"(response factors averaged over {report.n_cal_injections} contributing cal "
        f"injection(s) of {len(report.cal_files)})",
        f"- Recovery denominator policy: **{report.reference_policy}**"
        + (
            " (each file measured against the reference injection that most recently "
            "preceded it in acquisition order)"
            if report.reference_policy == "bracketed"
            else " (one reference injection for the whole run)"
        ),
        f"- Run status: **{report.status}**"
        + (f" - {report.status_reason}" if report.status_reason else ""),
        f"- Run-fit guard (`core.calibrate.rt_fit_quality`): "
        f"**{'fit_ok' if report.fit_ok else 'FIT NOT OK'}** - {report.fit_summary}",
    ]
    if report.rt_identity:
        lines.append(
            "- RT alignment: **identity** - every fitted channel sits inside `IDENTITY_TOL` "
            "of no-change, so the canonical retention times are used unaltered"
        )
    else:
        order = ("matched", "interpolated", "line")
        parts = [
            f"{report.rt_provenance.get(k, 0)} {k}"
            for k in order
            if report.rt_provenance.get(k, 0)
        ]
        lines.append(
            "- RT alignment provenance: "
            + ", ".join(parts)
            + " (`matched` = identified on a peak of this run's cal injection, "
            "`interpolated` = placed between two identified peaks, `line` = the channel's "
            "fitted line alone)"
        )
    lines += [
        f"- Native analytes with no response factor (cal-peak guard): "
        f"**{report.n_invalid_rrf}** of {report.n_natives}",
        f"- LOQ derived for **{report.loq_analytes}** analytes; "
        f"**{report.censored_cells}** non-zero cells censored to 0",
        f"- Detectability threshold: analyst mass >= {min_mass:g} ng",
    ]
    if report.aliases:
        pairs = ", ".join(f"`{k}` -> truth `{v}`" for k, v in sorted(report.aliases.items()))
        lines.append(f"- Truth columns re-paired (diagnostic): {pairs}")
    for note in report.notes:
        lines.append(f"- Note: {note}")

    if report.reference_batches:
        lines += ["", "#### Reference bracketing (recovery denominator per file)", ""]
        lines += _table(
            ["reference injection", "n files", "files bracketed by it"],
            [
                [f"`{ref_id}`", f"{len(members)}", ", ".join(f"`{m}`" for m in members)]
                for ref_id, members in sorted(report.reference_batches.items())
            ],
        )

    lines += ["", "### Uncorrected mass - native analytes", ""]
    lines += _mass_metrics_table(m_native)
    lines += ["", "### Uncorrected mass - labeled standards (13C surrogates / IS)", ""]
    lines += _mass_metrics_table(m_labeled)
    lines += ["", "### Surrogate recovery vs analyst", ""]
    lines += _recovery_metrics_table(rec)

    # Per-homolog recovery
    lines += ["", "#### Recovery error by homolog", ""]
    grouped_rec: dict[str, list[RecoveryRow]] = {}
    for r in report.recovery_rows:
        grouped_rec.setdefault(r.homolog, []).append(r)
    body = []
    for homolog in [h for h in HOMOLOG_ORDER if h in grouped_rec]:
        rs = grouped_rec[homolog]
        errs = np.array([r.error for r in rs])
        body.append([
            homolog, f"{len(rs)}", _num(float(np.mean(np.abs(errs)))),
            f"{np.mean([r.manual for r in rs]):.3f}", f"{np.mean([r.auto for r in rs]):.3f}",
            _num(float(np.max(np.abs(errs)))),
        ])
    lines += _table(
        ["homolog", "n", "MAE", "mean analyst", "mean auto", "max \\|error\\|"], body
    )

    # Per-homolog mass
    lines += ["", "### Uncorrected mass by homolog (native, detectable subset)", ""]
    grouped: dict[str, list[MassRow]] = {}
    for r in native:
        if r.manual >= min_mass:
            grouped.setdefault(r.homolog, []).append(r)
    order = [h for h in HOMOLOG_ORDER if h in grouped] + sorted(set(grouped) - set(HOMOLOG_ORDER))
    body = []
    for homolog in order:
        rs = grouped[homolog]
        rpds = np.abs([r.rpd for r in rs])
        ratio = sum(r.auto for r in rs) / sum(r.manual for r in rs)
        body.append([
            homolog, f"{len(rs)}", _pct(float(np.median(rpds))),
            _pct(float(np.mean(rpds <= BAND))), f"{ratio:.4f}",
        ])
    lines += _table(["homolog", "n", "median \\|RPD\\|", "within +/-25%", "total ratio"], body)

    # Worst 20
    lines += ["", "### Worst 20 (native, detectable subset, largest \\|RPD\\|)", ""]
    worst = sorted(
        (r for r in native if r.manual >= min_mass), key=lambda r: abs(r.rpd), reverse=True
    )[:20]
    lines += _table(
        ["#", "sample", "analyte", "homolog", "RT", "manual ng", "auto ng", "RPD", "flags",
         "method"],
        [
            [
                f"{i}", r.sample, r.label, r.homolog, f"{r.rt_apex:.3f}",
                f"{r.manual:.3f}", f"{r.auto:.3f}", f"{100 * r.rpd:+.0f}%", r.flags, r.method,
            ]
            for i, r in enumerate(worst, start=1)
        ],
    )

    # Miss breakdown
    detectable = [r for r in native if r.manual >= min_mass]
    misses = [r for r in detectable if abs(r.rpd) > BAND]
    lines += ["", "### Miss breakdown (native, detectable, outside +/-25%)", ""]
    if misses:
        counts: dict[str, list[MassRow]] = {}
        for r in misses:
            counts.setdefault(classify_miss(r, min_mass), []).append(r)
        manual_total = sum(r.manual for r in detectable)
        lines += _table(
            ["failure class", "n", "% of misses", "% of detectable rows",
             "% of detectable analyst mass"],
            [
                [
                    cls, f"{len(rs)}", f"{100 * len(rs) / len(misses):.0f}%",
                    f"{100 * len(rs) / len(detectable):.1f}%",
                    f"{100 * sum(r.manual for r in rs) / manual_total:.1f}%",
                ]
                for cls, rs in sorted(counts.items(), key=lambda kv: -len(kv[1]))
            ],
        )
        low = sum(1 for r in misses if r.rpd < 0)
        inside = [r for r in detectable if abs(r.rpd) <= BAND]
        ratio_inside = (
            sum(r.auto for r in inside) / sum(r.manual for r in inside) if inside else float("nan")
        )
        lines += [
            "",
            f"Direction: {low} of {len(misses)} misses read *low* against the analyst "
            f"({len(misses) - low} read high). "
            f"{sum(1 for r in misses if r.auto == 0)} produced exactly zero mass.",
            "",
            "Mass-weighted total ratio over the rows *inside* the band: "
            f"**{_num(ratio_inside)}** (vs {_num(m_native.total_ratio)} over all detectable "
            "rows). The gap between the two is the total-burden cost of the misses, as "
            "distinct from any systematic bias on the congeners that do agree.",
        ]
    else:
        lines.append("No detectable native row fell outside the band.")

    lines += render_rrf(report)
    lines += render_refs(report)

    # Coverage
    lines += ["", "### Coverage", ""]

    def _listing(title: str, items: list[str]) -> None:
        lines.append(f"- {title}: " + (", ".join(f"`{i}`" for i in items) if items else "none"))

    _listing("Data files that decoded to nothing", report.skipped_files)
    _listing("Truth columns with data but no scored data file", report.truth_only)
    _listing("Quantified data files with no truth column", report.file_only)
    _listing(
        "Analyst recoveries with no automated counterpart", report.recovery_unmatched
    )

    lines += [
        "",
        f"**Gate** (median \\|RPD\\| <= {_pct(GATE_MEDIAN_RPD)} AND recovery MAE <= "
        f"{GATE_RECOVERY_MAE:g}): **{'PASS' if passed else 'FAIL'}** - "
        f"median \\|RPD\\| {_pct(m_native.median_abs_rpd)}, recovery MAE {_num(rec.mae)}.",
        "",
    ]
    return lines


def render_params(params: IntegrationParams, qaqc: QaqcParams) -> list[str]:
    lines = ["## Parameters", "", "### Integration", ""]
    default = IntegrationParams()
    lines += _table(
        ["parameter", "value", "default"],
        [
            [f.name, f"{getattr(params, f.name)}", f"{getattr(default, f.name)}"]
            for f in fields(params)
        ],
    )
    changed = [f.name for f in fields(params) if getattr(params, f.name) != getattr(default, f.name)]
    lines += [
        "",
        "Changes from the defaults calibrated against analyst areas: "
        + (", ".join(f"`{c}`" for c in changed) if changed else "**none**")
        + ".",
        "",
        "### QA/QC",
        "",
    ]
    qdefault = QaqcParams()
    lines += _table(
        ["parameter", "value", "default"],
        [
            [f.name, f"{getattr(qaqc, f.name)}", f"{getattr(qdefault, f.name)}"]
            for f in fields(qaqc)
        ],
    )
    qchanged = [f.name for f in fields(qaqc) if getattr(qaqc, f.name) != getattr(qdefault, f.name)]
    lines += [
        "",
        "Changes from the defaults: "
        + (", ".join(f"`{c}`" for c in qchanged) if qchanged else "**none**")
        + ".",
        "",
    ]
    return lines


def render_summary(reports: list[RunReport], min_mass: float) -> list[str]:
    body = []
    rrf_body = []
    for report in reports:
        passed, m, r = gate(report, min_mass)
        label = report.run + (" (re-paired)" if report.aliases else "")
        body.append([
            label,
            report.status,
            f"{m.n_detectable}",
            _pct(m.median_abs_rpd),
            _pct(m.within_band),
            _num(m.total_ratio),
            f"{r.n}",
            _num(r.mae),
            "**PASS**" if passed else "**FAIL**",
        ])
        # The cal/ref sheets are joined by injection, so a run scored twice (as labelled and
        # re-paired) would contribute the same row twice: sample aliasing cannot change it.
        if not report.aliases:
            rm = rrf_metrics(report.rrf_rows)
            ref_rpds = np.abs([x.rpd for x in report.ref_rows]) if report.ref_rows else None
            rrf_body.append([
                report.run,
                f"{rm.n}",
                _pct(rm.median_abs_rpd),
                f"{rm.n_off}",
                _pct(rm.median_abs_area_rpd),
                _pct(report.analyst_rrf_spread),
                f"{len(report.ref_rows)}",
                _pct(float(np.median(ref_rpds))) if ref_rpds is not None else "n/a",
            ])
    lines = (
        ["## Verdict", ""]
        + _table(
            ["run", "status", "detectable rows", "mass median \\|RPD\\|", "within +/-25%",
             "total ratio", "recovery pairs", "recovery MAE", "gate"],
            body,
        )
        + [
            "",
            f"Gate: median \\|RPD\\| of uncorrected native mass <= {_pct(GATE_MEDIAN_RPD)} "
            f"**and** surrogate-recovery MAE <= {GATE_RECOVERY_MAE:g}, on "
            "the run named on the command line.",
            "",
        ]
    )
    if rrf_body:
        lines += [
            "### Calibration and reference injections",
            "",
            "The two files every mass in a run is computed from, and the two the analyst's "
            "workbook holds no sample sheet for, so nothing in the table above scores them "
            "directly. Not gated: they are diagnostics of where a mass disagreement comes "
            "from, and the analyst's own reproducibility floor is printed beside them.",
            "",
        ]
        lines += _table(
            ["run", "RRF n", "RRF median \\|RPD\\|", "RRF outside +/-25%",
             "cal-area median \\|RPD\\|", "analyst's own spread", "ref pairs",
             "ref median \\|RPD\\|"],
            rrf_body,
        )
        lines.append("")
    return lines


def render(
    reports: list[RunReport], params: IntegrationParams, qaqc: QaqcParams,
    min_mass: float, command: str = "",
) -> str:
    lines = [
        "# Masses and recoveries: automated against analyst",
        "",
        f"Generated {date.today().isoformat()} by `tools/bench_quant.py`"
        + (f":\n\n    {command}\n" if command else "."),
        "",
        "> Every table below is generated. Prose under an **Interpretation** heading is "
        "hand-written from those tables and is *not* regenerated. Re-running the command "
        "overwrites this file, so re-attach it.",
        "",
        "The whole run goes through `core.pipeline.run_analysis`, the production entry point, "
        "with production defaults. Its **uncorrected** masses "
        "(`RRF * area * (IS_spike / IS_area)`) are joined to the analyst workbook's "
        "\"Mass (ng)\" column and its per-homolog surrogate recoveries to the workbook "
        "Summary sheet's Recoveries block (`report.truth.read_run_truth`).",
        "",
        "`RPD = (auto - manual) / ((auto + manual) / 2)` on masses; recoveries are scored as "
        "**absolute error in recovery fraction**, since a recovery is already a ratio. "
        "Headline metrics cover native analytes whose *analyst* mass clears the detectability "
        "threshold; labeled standards and sub-threshold rows are reported separately rather "
        "than folded in.",
        "",
        "The analyst workbook holds no worksheet for the calibration or reference injections, "
        "so those files are quantified by the pipeline (they are what the response factors "
        "and the recovery denominator are made of) but are not scored here. Hexane carryover "
        "blanks are not quantified at all.",
        "",
    ]
    lines += render_summary(reports, min_mass)
    lines += render_params(params, qaqc)
    for report in reports:
        lines += render_run(report, min_mass)
    return "\n".join(lines).rstrip() + "\n"


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------


def _apply_overrides(params, overrides: list[str]):
    """Apply ``name=value`` overrides to a params dataclass, coerced to the field's type.

    A list of labels is given comma-separated, and empty to switch its rule off:
    ``--param forced_riders=83,117`` or ``--param forced_present=``.
    """
    known = {f.name for f in fields(params)}
    for item in overrides:
        name, _, raw = item.partition("=")
        if name not in known:
            raise SystemExit(f"unknown {type(params).__name__} field: {name}")
        current = getattr(params, name)
        if isinstance(current, bool):
            value = raw.strip().lower() in {"1", "true", "yes"}
        elif isinstance(current, int):
            value = int(raw)
        elif isinstance(current, tuple):
            value = tuple(part.strip() for part in raw.split(",") if part.strip())
        else:
            value = float(raw)
        object.__setattr__(params, name, value)
    return params


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", nargs="*",
                    help="run folder name(s) under --root (default: every run folder found). "
                         "Suffix a run with '@repaired' to score it against the re-paired "
                         "truth columns given by --alias. The same folder may appear twice, "
                         "e.g. RUN RUN@repaired.")
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT,
                    help="directory holding the analysis-run folders")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT,
                    help="markdown report to write (relative to --root); '-' to skip")
    ap.add_argument("--min-mass", type=float, default=DEFAULT_MIN_MASS,
                    help="analyst mass (ng) below which a row is treated as noise")
    ap.add_argument("--param", action="append", default=[], metavar="NAME=VALUE",
                    help="override an IntegrationParams field (repeatable)")
    ap.add_argument("--qaqc", action="append", default=[], metavar="NAME=VALUE",
                    help="override a QaqcParams field (repeatable)")
    ap.add_argument("--sweep", metavar="FIELD=V1,V2,...",
                    help="instead of the full report, re-integrate each run's calibration "
                         "injections at each value of an IntegrationParams field and score "
                         "their areas against the analyst's `cal` sheets. This is the "
                         "threshold-selection evidence for `overlap_frac`; a value above 1 "
                         "disables clustering and gives the pre-cluster engine. Example: "
                         "--sweep overlap_frac=10,0.02,0.05,0.10,0.15")
    ap.add_argument("--alias", action="append", default=[], metavar="FILEID=TRUTHCOL",
                    help="score the data file with id FILEID against truth column TRUTHCOL "
                         "(diagnostic, for provably mis-paired ground truth; repeatable). "
                         "Applies only to runs suffixed '@repaired'.")
    args = ap.parse_args(argv)

    params = _apply_overrides(IntegrationParams(), args.param)
    qaqc = _apply_overrides(QaqcParams(), args.qaqc)
    aliases: dict[str, str] = {}
    for item in args.alias:
        file_id, sep, truth_col = item.partition("=")
        if not sep or not file_id or not truth_col:
            raise SystemExit(f"--alias expects FILEID=TRUTHCOL, got {item!r}")
        aliases[file_id] = truth_col

    def _resolve(run: str) -> Path:
        folder = Path(run)
        if not folder.is_dir():
            folder = args.root / run
        if not folder.is_dir():
            raise SystemExit(f"run folder not found: {run}")
        return folder

    if args.sweep:
        field, sep, raw = args.sweep.partition("=")
        if not sep:
            raise SystemExit(f"--sweep expects FIELD=V1,V2,..., got {args.sweep!r}")
        if field not in {f.name for f in fields(IntegrationParams)}:
            raise SystemExit(f"unknown IntegrationParams field: {field}")
        try:
            values = [float(v) for v in raw.split(",") if v.strip()]
        except ValueError:
            raise SystemExit(f"--sweep values must be numbers, got {raw!r}") from None
        if not values:
            raise SystemExit("--sweep needs at least one value")
        command = "python tools/bench_quant.py " + " ".join(
            f'"{a}"' if "@" in a or " " in a else a
            for a in (argv if argv is not None else sys.argv[1:])
        )
        lines = [
            f"## Calibration-area agreement across `{field}`",
            "",
            f"Generated by `tools/bench_quant.py`:\n\n    {command}\n",
            "",
        ]
        for spec in args.runs or discover_runs(Path(args.root)):
            folder = _resolve(spec.partition("@")[0])
            print(f"[bench] sweeping {folder.name} ...", file=sys.stderr, flush=True)
            lines += sweep_cal_areas(folder, field, values, params)
        text = "\n".join(lines).rstrip() + "\n"
        print(text)
        return 0

    reports: list[RunReport] = []
    for spec in args.runs or discover_runs(Path(args.root)):
        run, sep, suffix = spec.partition("@")
        tokens = [t for t in suffix.split("+") if t] if sep else []
        run_aliases: dict[str, str] = {}
        for token in tokens:
            if token == "repaired":
                if not aliases:
                    raise SystemExit(f"{spec!r} asks for '@repaired' but no --alias was given")
                run_aliases = aliases
            else:
                raise SystemExit(f"unknown run modifier {token!r} in {spec!r}")
        folder = _resolve(run)
        note = "re-paired truth" if run_aliases else "as labelled"
        print(f"[bench] {folder.name} ({note}) ...", file=sys.stderr, flush=True)
        reports.append(bench_run(folder, params, qaqc, run_aliases))

    command = "python tools/bench_quant.py " + " ".join(
        f'"{a}"' if "@" in a or " " in a else a
        for a in (argv if argv is not None else sys.argv[1:])
    )
    text = render(reports, params, qaqc, args.min_mass, command)
    print(text)
    if str(args.out) != "-":
        out = args.out if args.out.is_absolute() else args.root / args.out
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text, encoding="utf-8")
        print(f"[bench] wrote {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
