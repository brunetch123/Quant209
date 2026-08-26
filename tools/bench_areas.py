"""Benchmark automated peak areas against the analyst's manual integration.

The engine's integrated areas must agree with the areas the analyst typed into the
quantification workbook. Areas are the first link of the calculation chain (area ->
RRF-scaled mass -> recovery correction), so every later number inherits whatever bias
lives here. This benchmark is the gate that bias has to clear.

What is compared
----------------
RT anchors are built from the *first* calibration injection (``Cal01``) and reused for
every other file in the run. Replicate cal injections are averaged for *RRFs* but never
for anchors, so a single injection is used.

That same calibration injection also supplies the run's RT *map*
(``core.calibrate.align_rt_map``) before anchors are measured, and answers the run-fit
guard afterwards. That is the default here (``--rt-source aligned``, reported on every
run section). This tool and ``tests/test_calibrate.py`` are presently the only things
exercising it. ``@canonical`` reproduces the unaligned behaviour and ``@workbook`` the
analyst-table control (see :data:`RT_SOURCES`).

Every file whose detected type is ``sample`` or ``blank`` is then integrated and joined to
the workbook areas on (sample id, analyte label). Deliberately excluded:

    - ``hex``: hexane carryover checks. The analyst never quantifies them.
    - ``cal``: the calibration injections themselves. They define the anchors.
    - ``ref``: ``read_run_truth`` returns *no* ``Ref04/05/06`` columns for either
      benchmark run (the reference injection's own worksheet is workbook scaffolding,
      excluded by ``truth._NON_SAMPLE``), so there is nothing to join a reference file
      against.

Metrics
-------
``RPD = (auto - manual) / ((auto + manual) / 2)``, the symmetric relative percent
difference.

Headline metrics are computed on *native* analytes only, over the **detectable
subset**: rows whose manual area is at least ``--min-area`` counts·s. Below that
threshold the analyst is integrating noise, and an RPD between two noise measurements
says nothing about the engine. Labelled standards (13C surrogates, injection standards)
are reported as a separate table. The workbook records areas for them too, and they
are a useful independent check, but they are single dominant peaks on quiet channels
and would flatter the headline number.

The area-weighted total ratio (sum(auto)/sum(manual)) is the number that matters for
total PCB burden. Median |RPD| can look fine while a systematic bias on the biggest
peaks walks the totals off, and vice versa.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from dataclasses import dataclass, fields, replace
from datetime import date
from pathlib import Path

import numpy as np
import openpyxl

# Allow ``python tools/bench_areas.py`` from a source checkout without installing.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quant209.core.calibrate import (  # noqa: E402
    RunFitReport,
    align_rt_map,
    build_anchors,
    rt_fit_quality,
)
from quant209.core.config import Analyte, load_canonical_peaks  # noqa: E402
from quant209.core.integrate import integrate_run  # noqa: E402
from quant209.core.mzml import read_srm_mzml  # noqa: E402
from quant209.core.peaks import IntegrationParams, PeakResult  # noqa: E402
from quant209.core.sequence import FileEntry, find_data_files, scan_inputs  # noqa: E402
from quant209.report.truth import read_run_truth  # noqa: E402
from quant209.report.validation import (  # noqa: E402
    BAND,
    GATE_MEDIAN_RPD,
    GATE_WITHIN_FRAC,
    HOMOLOG_ORDER,
    NEAR_LOQ_MULTIPLE,
    AreaRow as Row,
    Metrics,
    classify_area_miss as classify_miss,
    join_areas,
    metrics,
)

DEFAULT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = Path("docs/validation/areas.md")

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


#: File types the analyst quantifies, and which this benchmark therefore scores.
QUANTIFIED_TYPES = ("sample", "blank")

#: Manual area (counts·s) below which a row is treated as noise, not a measurement.
DEFAULT_MIN_AREA = 1000.0

# BAND, GATE_MEDIAN_RPD, GATE_WITHIN_FRAC, NEAR_LOQ_MULTIPLE and HOMOLOG_ORDER live in
# ``quant209.report.validation`` (imported above), identical values shared with every
# other benchmark so a threshold changed in one place changes everywhere it is used.


@dataclass
class Skipped:
    """A data file that could not be scored, and why."""

    name: str
    reason: str


@dataclass
class RunReport:
    """Everything measured for one analysis run."""

    run: str
    cal_file: str
    rt_source: str
    fit: RunFitReport
    # How each analyte's expected RT was arrived at when ``rt_source`` is ``aligned``:
    # ``RtAlignment.source_of`` counted over the canonical table. Empty on the other sources,
    # which do not align. This is the auditor's line: "identified on a peak" is direct
    # evidence, "interpolated" is inference from neighbours, "line" is the channel fit alone.
    rt_provenance: dict[str, int]
    rt_identity: bool  # the alignment came out as the identity: canonical RTs left untouched
    aliases: dict[str, str]
    rows: list[Row]
    skipped: list[Skipped]
    truth_only: list[str]          # truth columns with data but no data file
    truth_empty: list[str]         # truth columns that are entirely blank
    file_only: list[str]           # quantified data files with no truth column
    unmatched_labels: list[str]    # truth analyte rows with no canonical analyte
    n_files: int


# --------------------------------------------------------------------------------------
# Retention-time table
# --------------------------------------------------------------------------------------

#: Where a run's pre-anchor retention times come from.
#:
#: ``aligned``: the shipped canonical table, mapped onto this run by
#:     ``core.calibrate.align_rt_map`` from its own calibration injection. **This tool's
#:     default.** No production caller wires it in yet, so an unqualified benchmark
#:     number measured through it is the closest thing to production behaviour there
#:     currently is.
#: ``canonical``: the canonical table applied raw, with no alignment, kept so the two
#:     can be compared on one command line (``RUN RUN@canonical``).
#: ``workbook``: the run's own analyst "Cal RT" column (see :func:`workbook_cal_rts`).
#:     The diagnostic control: the retention times the analyst actually used, which
#:     bound what any alignment could achieve.
RT_SOURCES = ("aligned", "canonical", "workbook")

#: Per-sample sheet coordinates of the analyst's calibration retention time. Column C is the
#: analyte label and column D the "Cal RT" the analyst calibrated that run against (row 15
#: is the header. See ``report.truth`` for the rest of the layout).
_WB_COL_LABEL, _WB_COL_CAL_RT, _WB_DATA_START = 3, 4, 16


def workbook_cal_rts(folder: Path) -> dict[str, float]:
    """Read one run's own analyst "Cal RT" column: ``{analyte label: RT minutes}``.

    DIAGNOSTIC ONLY: a control, not a workaround. ``core.config``'s canonical ``cal_rt``
    values are method- and column-specific. They were extracted from a single run's
    calibration sheet, and ``build_anchors`` refines them by well under a minute, enough
    for column-build drift but not for a different GC oven program. Reading a run's own
    table separates "the integration engine disagrees with the analyst" from "the engine
    was pointed at the wrong retention times", which are otherwise indistinguishable in
    the headline metrics.

    ``core.calibrate.align_rt_map`` derives the run's retention times from its own
    calibration injection (``--rt-source aligned``, the default). This function stays as
    the ceiling that alignment is measured against (the retention times the analyst
    themself worked to), and must never be a production path, since it reads them out of
    the very workbook the benchmark scores against.
    """
    rts: dict[str, float] = {}
    for xlsx in sorted(folder.glob("*.xlsx")):
        if xlsx.name.startswith("~$") or "sequence" in xlsx.name.lower():
            continue
        wb = openpyxl.load_workbook(xlsx, data_only=True, read_only=True)
        try:
            for ws in wb.worksheets:
                header = next(
                    ws.iter_rows(min_row=15, max_row=15, max_col=11, values_only=True), ()
                )
                if not any(isinstance(v, str) and "Cal RT" in v for v in header if v):
                    continue
                for row in ws.iter_rows(
                    min_row=_WB_DATA_START, max_col=_WB_COL_CAL_RT, values_only=True
                ):
                    label, rt = row[_WB_COL_LABEL - 1], row[_WB_COL_CAL_RT - 1]
                    if label is None:
                        continue
                    try:
                        rts[str(label)] = float(rt)  # type: ignore[arg-type]
                    except (TypeError, ValueError):
                        continue
                if rts:
                    return rts  # one sheet is enough: Cal RT is per-run, not per-sample
        finally:
            wb.close()
    return rts


def apply_rt_table(analytes: list[Analyte], rts: dict[str, float]) -> list[Analyte]:
    """Return ``analytes`` with ``cal_rt`` replaced wherever ``rts`` has a value."""
    return [replace(a, cal_rt=rts[a.label]) if a.label in rts else a for a in analytes]


# --------------------------------------------------------------------------------------
# Measurement
# --------------------------------------------------------------------------------------

def _usable(chroms: dict) -> bool:
    """Whether a parsed mzML yielded anything integrable.

    Several conversions in this dataset were interrupted, leaving a file that is valid
    XML down to a point and then simply stops. ``read_srm_mzml`` tolerates that and
    returns the prefix it managed to decode, which can be nothing at all. Emptiness is
    the check, not file size.
    """
    return any(c.times.size > 1 for c in chroms.values())


def _integrate_file(
    path: Path, analytes: list[Analyte], anchors: dict[str, float], params: IntegrationParams
) -> tuple[list[PeakResult], int]:
    """Integrate one file. Returns its peaks and how many analytes went unplaced.

    ``integrate_run`` warns (rather than raises) when analytes fall out: no anchor, or the
    transition missing from a partially converted file. Here that count is data: a file that
    silently lost half the canonical table should show up in the report, not in a log nobody
    reads.
    """
    chroms = read_srm_mzml(path)
    if not _usable(chroms):
        return [], len(analytes)
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        results = integrate_run(chroms, analytes, anchors, params)
    del chroms
    return results, len(analytes) - len({r.label for r in results})


def bench_run(
    folder: Path,
    params: IntegrationParams,
    rt_source: str = "aligned",
    aliases: dict[str, str] | None = None,
) -> RunReport:
    """Integrate every quantified file in ``folder`` and join the areas to analyst truth.

    Every joined row is kept, sub-threshold ones included: the detectability threshold is a
    *reporting* choice applied in :func:`metrics`, so it can be varied (``--min-area``)
    without re-integrating 20 files.

    ``rt_source`` selects where the pre-anchor retention times come from. See
    :data:`RT_SOURCES`. The default ``aligned`` is the production path.

    ``aliases`` maps a data file's sample id to the truth column it should be scored
    against. DIAGNOSTIC ONLY, for the case where the *ground truth* is provably attached to
    the wrong file (the archive holds one sample/blank pair that looks swapped): it lets
    the corrected pairing be measured and reproduced without editing anyone's workbook. It
    must never be used to shop for a better-agreeing column.
    """
    aliases = aliases or {}
    analytes = load_canonical_peaks()
    if rt_source == "workbook":
        rts = workbook_cal_rts(folder)
        if not rts:
            raise SystemExit(f"{folder.name}: no analyst Cal RT column found in any workbook")
        analytes = apply_rt_table(analytes, rts)
    by_label = {a.label: a for a in analytes}

    entries: list[FileEntry] = scan_inputs(find_data_files(folder))
    skipped: list[Skipped] = []

    cal_entries = [e for e in entries if e.file_type == "cal"]
    if not cal_entries:
        raise SystemExit(f"{folder.name}: no calibration injection found")
    cal = cal_entries[0]  # Cal01 (see module docstring)
    cal_chroms = read_srm_mzml(cal.path)
    if not _usable(cal_chroms):
        raise SystemExit(f"{folder.name}: calibration file {cal.name} parsed to no chromatograms")
    # ``aligned`` is the production path: map the canonical retention times onto this run
    # before anchoring. The other two sources deliberately skip it: ``canonical`` to
    # reproduce the pre-alignment behaviour, ``workbook`` because the analyst's own table is
    # already in this run's retention-time frame.
    alignment = align_rt_map(cal_chroms, analytes, params) if rt_source == "aligned" else None
    anchors = build_anchors(cal_chroms, analytes, params, alignment=alignment)
    provenance: dict[str, int] = {}
    if alignment is not None:
        for analyte in analytes:
            key = alignment.source_of(analyte)
            provenance[key] = provenance.get(key, 0) + 1
    # The run-fit guard, reported rather than enforced: this tool exists to measure what the
    # engine does, including when it is pointed at the wrong retention times.
    fit = rt_fit_quality(cal_chroms, analytes, anchors, params)
    del cal_chroms

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # unreadable-workbook warnings are re-reported below
        truth = read_run_truth(folder)

    targets = [e for e in entries if e.file_type in QUANTIFIED_TYPES]
    truth_cols = set(truth.areas.columns)
    populated = {c for c in truth.areas.columns if truth.areas[c].notna().any()}

    # Files whose truth column has no data are skipped *before* integration: doing so would
    # be expensive I/O for a row nothing would join against. Everything that does get
    # integrated is handed to ``join_areas`` (``quant209.report.validation``) afterwards, in
    # the same order, to build the actual Row objects: that function is shared with the
    # per-run report generator and is not itself an I/O decision, so it takes peaks already
    # in hand rather than a folder.
    peaks_by_file: dict[str, list[PeakResult]] = {}
    file_types: dict[str, str] = {}
    scored: set[str] = set()
    file_only: list[str] = []

    for entry in sorted(targets, key=lambda e: e.sample_id):
        truth_id = aliases.get(entry.sample_id, entry.sample_id)
        if truth_id not in populated:
            file_only.append(entry.sample_id)
            continue
        peaks, missing = _integrate_file(entry.path, analytes, anchors, params)
        if not peaks:
            skipped.append(
                Skipped(entry.name, "no chromatograms decoded (truncated/header-only conversion)")
            )
            continue
        if missing:
            skipped.append(
                Skipped(
                    entry.name,
                    f"partially decoded: {missing} of {len(analytes)} analytes had no trace "
                    f"(scored on the {len(analytes) - missing} that did)",
                )
            )
        scored.add(entry.sample_id)
        peaks_by_file[entry.sample_id] = peaks
        file_types[entry.sample_id] = entry.file_type

    rows = join_areas(peaks_by_file, file_types, analytes, truth.areas, aliases)

    return RunReport(
        run=folder.name,
        cal_file=cal.name,
        rt_source=rt_source,
        fit=fit,
        rt_provenance=provenance,
        rt_identity=alignment.is_identity if alignment is not None else False,
        aliases=dict(aliases),
        rows=rows,
        skipped=skipped,
        truth_only=sorted(
            populated - {aliases.get(e.sample_id, e.sample_id) for e in targets}
        ),
        truth_empty=sorted(truth_cols - populated),
        file_only=sorted(file_only),
        unmatched_labels=sorted(set(truth.areas.index) - set(by_label)),
        n_files=len(scored),
    )


# --------------------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------------------

# Metrics, metrics() and classify_miss() live in ``quant209.report.validation`` (imported
# above as ``Metrics``, ``metrics`` and ``classify_area_miss as classify_miss``), shared
# unchanged with the other benchmarks.


def _split(rows: list[Row]) -> tuple[list[Row], list[Row]]:
    """(native rows, labelled-standard rows)."""
    native = [r for r in rows if r.role == "native"]
    return native, [r for r in rows if r.role != "native"]


# --------------------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------------------

def _pct(x: float) -> str:
    return "n/a" if not np.isfinite(x) else f"{100 * x:.1f}%"


def _table(header: list[str], rows: list[list[str]]) -> list[str]:
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return out


def _metrics_table(name: str, m: Metrics) -> list[str]:
    return _table(
        ["metric", name],
        [
            ["joined rows", f"{m.n_joined}"],
            ["detectable rows (manual area >= threshold)", f"{m.n_detectable}"],
            ["median \\|RPD\\|", _pct(m.median_abs_rpd)],
            ["within +/-25%", _pct(m.within_band)],
            ["area-weighted total ratio (sum auto / sum manual)",
             "n/a" if not np.isfinite(m.total_ratio) else f"{m.total_ratio:.4f}"],
            ["median signed RPD (bias)", _pct(m.median_rpd)],
        ],
    )


def render_run(report: RunReport, min_area: float) -> list[str]:
    native, labeled = _split(report.rows)
    m_native = metrics(native, min_area)
    m_labeled = metrics(labeled, min_area)

    diagnostics = []
    if report.rt_source != "aligned":
        diagnostics.append(f"RT table from {report.rt_source}")
    if report.aliases:
        diagnostics.append("truth columns re-paired")
    title = report.run + (f" (diagnostic: {', '.join(diagnostics)})" if diagnostics else "")
    lines: list[str] = [f"## {title}", ""]
    lines += [
        f"- Files scored: **{report.n_files}** ({', '.join(QUANTIFIED_TYPES)} types)",
        f"- RT anchors refined from: `{report.cal_file}`",
        f"- Pre-anchor retention times: **{report.rt_source}**",
    ]
    if report.rt_source == "aligned":
        if report.rt_identity:
            lines.append(
                "- RT alignment: **identity** - every fitted channel's line sits inside "
                "`IDENTITY_TOL` of no-change, so the canonical retention times are used "
                "unaltered and the anchors are those of `canonical`"
            )
        else:
            # ``core.calibrate.RtAlignment.source_of``, counted over the canonical table.
            order = ("matched", "interpolated", "line")
            parts = [
                f"{report.rt_provenance.get(k, 0)} {k}"
                for k in order
                if report.rt_provenance.get(k, 0)
            ]
            lines.append(
                "- RT alignment provenance (how each analyte's expected RT was arrived at): "
                + ", ".join(parts)
                + " (`matched` = identified on a peak of this run's cal injection, "
                "`interpolated` = placed between two identified peaks, `line` = the "
                "channel's fitted line alone)"
            )
    lines += [
        f"- Run-fit guard (`core.calibrate.rt_fit_quality`, on the cal injection): "
        f"**{'fit_ok' if report.fit.fit_ok else 'FIT NOT OK'}** - {report.fit.summary}",
        f"- Detectability threshold: manual area >= {min_area:g} counts*s",
    ]
    if report.aliases:
        pairs = ", ".join(f"`{k}` -> truth `{v}`" for k, v in sorted(report.aliases.items()))
        lines.append(f"- Truth columns re-paired (diagnostic): {pairs}")
    lines += [
        "",
        "### Headline - native analytes",
        "",
    ]
    lines += _metrics_table("value", m_native)
    lines += ["", "### Labeled standards (13C surrogates / injection standards)", ""]
    lines += _metrics_table("value", m_labeled)

    # Per-homolog
    lines += ["", "### Per-homolog (native, detectable subset)", ""]
    grouped: dict[str, list[Row]] = {}
    for r in native:
        if r.manual >= min_area:
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
        (r for r in native if r.manual >= min_area), key=lambda r: abs(r.rpd), reverse=True
    )[:20]
    # ``RT`` is the apex the engine actually integrated at. It is the single most diagnostic
    # column here: an RT far from the analyte's expected elution means the run was scored
    # against the wrong retention-time table, which looks identical to bad integration in
    # every other column.
    lines += _table(
        ["#", "sample", "analyte", "homolog", "RT", "manual", "auto", "RPD", "flags",
         "method", "S/N"],
        [
            [
                f"{i}", r.sample, r.label, r.homolog, f"{r.rt_apex:.3f}",
                f"{r.manual:.0f}", f"{r.auto:.0f}", f"{100 * r.rpd:+.0f}%", r.flags, r.method,
                f"{r.sn:.1f}",
            ]
            for i, r in enumerate(worst, start=1)
        ],
    )

    # Where the misses live
    detectable = [r for r in native if r.manual >= min_area]
    misses = [r for r in detectable if abs(r.rpd) > BAND]
    lines += ["", "### Miss breakdown (native, detectable, outside +/-25%)", ""]
    if misses:
        counts: dict[str, list[Row]] = {}
        for r in misses:
            counts.setdefault(classify_miss(r, min_area), []).append(r)
        manual_total = sum(r.manual for r in detectable)
        lines += _table(
            ["failure class", "n", "% of misses", "% of detectable rows",
             "% of detectable manual area"],
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
            f"{sum(1 for r in misses if r.auto == 0)} integrated to exactly zero area.",
            "",
            "Area-weighted total ratio over the rows *inside* the band: "
            f"**{ratio_inside:.4f}** (vs {m_native.total_ratio:.4f} over all detectable rows) "
            "The gap between the two is the total-burden cost of the misses, as distinct "
            "from any systematic bias on the peaks that do agree.",
        ]
    else:
        lines.append("No detectable native row fell outside the band.")

    # Coverage / skips
    lines += ["", "### Coverage", ""]
    def _listing(title: str, items: list[str]) -> None:
        lines.append(f"- {title}: " + (", ".join(f"`{i}`" for i in items) if items else "none"))

    if report.skipped:
        lines.append("- Skipped or partial data files:")
        for s in report.skipped:
            lines.append(f"  - `{s.name}` - {s.reason}")
    else:
        lines.append("- Skipped or partial data files: none")
    _listing("Truth columns with data but no data file", report.truth_only)
    _listing("Truth columns present but entirely blank", report.truth_empty)
    _listing("Quantified data files with no truth column", report.file_only)
    _listing("Truth analyte rows outside the canonical table", report.unmatched_labels)

    # Gate
    passed = (
        np.isfinite(m_native.median_abs_rpd)
        and m_native.median_abs_rpd <= GATE_MEDIAN_RPD
        and m_native.within_band >= GATE_WITHIN_FRAC
    )
    lines += [
        "",
        f"**Gate** (median \\|RPD\\| <= {_pct(GATE_MEDIAN_RPD)} AND >= {_pct(GATE_WITHIN_FRAC)} "
        f"within +/-25%): **{'PASS' if passed else 'FAIL'}** - "
        f"median \\|RPD\\| {_pct(m_native.median_abs_rpd)}, "
        f"within band {_pct(m_native.within_band)}.",
        "",
    ]
    return lines


def render_params(params: IntegrationParams) -> list[str]:
    default = IntegrationParams()
    changed = [
        f.name for f in fields(params)
        if getattr(params, f.name) != getattr(default, f.name)
    ]
    lines = ["## Integration parameters", ""]
    lines += _table(
        ["parameter", "value", "default"],
        [
            [f.name, f"{getattr(params, f.name)}", f"{getattr(default, f.name)}"]
            for f in fields(params)
        ],
    )
    lines += [
        "",
        "Changes from the defaults calibrated against analyst areas: "
        + (", ".join(f"`{c}`" for c in changed) if changed else "**none**")
        + ".",
        "",
    ]
    return lines


def render_summary(reports: list[RunReport], min_area: float) -> list[str]:
    """One row per benchmarked run: the whole verdict at a glance."""
    body = []
    for report in reports:
        m = metrics(_split(report.rows)[0], min_area)
        passed = (
            np.isfinite(m.median_abs_rpd)
            and m.median_abs_rpd <= GATE_MEDIAN_RPD
            and m.within_band >= GATE_WITHIN_FRAC
        )
        label = report.run + (" (re-paired)" if report.aliases else "")
        body.append([
            label, report.rt_source, "yes" if report.fit.fit_ok else "**no**",
            f"{m.n_detectable}",
            _pct(m.median_abs_rpd), _pct(m.within_band),
            "n/a" if not np.isfinite(m.total_ratio) else f"{m.total_ratio:.4f}",
            "**PASS**" if passed else "**FAIL**",
        ])
    return (
        ["## Verdict", ""]
        + _table(
            ["run", "RT source", "fit_ok", "detectable rows", "median \\|RPD\\|",
             "within +/-25%", "total ratio", "gate"],
            body,
        )
        + [
            "",
            f"Gate: median \\|RPD\\| <= {_pct(GATE_MEDIAN_RPD)} **and** "
            f">= {_pct(GATE_WITHIN_FRAC)} of detectable native rows within +/-25%.",
            "",
        ]
    )


def render(
    reports: list[RunReport], params: IntegrationParams, min_area: float, command: str = ""
) -> str:
    lines = [
        "# Peak areas: automated against analyst",
        "",
        f"Generated {date.today().isoformat()} by `tools/bench_areas.py`"
        + (f":\n\n    {command}\n" if command else "."),
        "",
        "> Every table below is generated. Prose under an **Interpretation** heading is "
        "hand-written from those tables and is *not* regenerated. Re-running the command "
        "overwrites this file, so re-attach it.",
        "",
        "Automated integration (`core.integrate.integrate_run`, RT anchors from the run's "
        "first calibration injection) is compared against the areas the analyst manually "
        "integrated and typed into the quantification workbook "
        "(`report.truth.read_run_truth`).",
        "",
        "`RPD = (auto - manual) / ((auto + manual) / 2)`. Headline metrics cover native "
        "analytes whose *manual* area clears the detectability threshold; labeled standards "
        "and sub-threshold rows are reported separately rather than folded in.",
        "",
        "The analyst workbook holds no worksheet for the reference injections "
        "(`Ref04/05/06`), so reference files are not scored; hexane carryover blanks and the "
        "calibration injections are likewise out of scope.",
        "",
    ]
    lines += render_summary(reports, min_area)
    lines += render_params(params)
    for report in reports:
        lines += render_run(report, min_area)
    return "\n".join(lines).rstrip() + "\n"


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------

def _apply_overrides(params: IntegrationParams, overrides: list[str]) -> IntegrationParams:
    """Apply ``name=value`` parameter overrides, coerced to the field's declared type."""
    known = {f.name for f in fields(params)}
    for item in overrides:
        name, _, raw = item.partition("=")
        if name not in known:
            raise SystemExit(f"unknown IntegrationParams field: {name}")
        current = getattr(params, name)
        if isinstance(current, bool):
            value = raw.strip().lower() in {"1", "true", "yes"}
        elif isinstance(current, int):
            value = int(raw)
        else:
            value = float(raw)
        setattr(params, name, value)
    return params


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", nargs="*",
                    help="run folder name(s) under --root (default: every run folder found). "
                         "Suffix a run with '@' and '+'-joined modifiers to control that one "
                         "section: 'aligned'/'canonical'/'workbook' set its RT source "
                         "(overriding --rt-source) and 'repaired' applies --alias. The same "
                         "folder may appear more than once, e.g. RUN RUN@canonical.")
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT,
                    help="directory holding the analysis-run folders")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT,
                    help="markdown report to write (relative to --root); '-' to skip")
    ap.add_argument("--min-area", type=float, default=DEFAULT_MIN_AREA,
                    help="manual area (counts*s) below which a row is treated as noise")
    ap.add_argument("--param", action="append", default=[], metavar="NAME=VALUE",
                    help="override an IntegrationParams field (repeatable)")
    ap.add_argument("--alias", action="append", default=[], metavar="FILEID=TRUTHCOL",
                    help="score the data file with sample id FILEID against truth column "
                         "TRUTHCOL (diagnostic, for provably mis-paired ground truth; "
                         "repeatable). Applies only to runs suffixed '@repaired'.")
    ap.add_argument("--rt-source", choices=RT_SOURCES, default=RT_SOURCES[0],
                    help="pre-anchor retention times: 'aligned' (default, the production "
                         "path, canonical table mapped onto the run by "
                         "core.calibrate.align_rt_map), 'canonical' (raw canonical table, "
                         "the pre-alignment behaviour) or 'workbook' (the run's own analyst "
                         "Cal RT column, the diagnostic control, see workbook_cal_rts)")
    args = ap.parse_args(argv)

    params = _apply_overrides(IntegrationParams(), args.param)
    aliases: dict[str, str] = {}
    for item in args.alias:
        file_id, sep, truth_col = item.partition("=")
        if not sep or not file_id or not truth_col:
            raise SystemExit(f"--alias expects FILEID=TRUTHCOL, got {item!r}")
        aliases[file_id] = truth_col
    reports = []
    for spec in args.runs or discover_runs(Path(args.root)):
        # "RUN", "RUN@canonical", "RUN@repaired", "RUN@workbook+repaired": the same folder
        # may appear more than once under different diagnostics, which is how one command
        # produces the as-measured result and its diagnostics side by side.
        run, sep, suffix = spec.partition("@")
        tokens = [t for t in suffix.split("+") if t] if sep else []
        rt_source = args.rt_source
        run_aliases: dict[str, str] = {}
        for token in tokens:
            if token in RT_SOURCES:
                rt_source = token
            elif token == "repaired":
                if not aliases:
                    raise SystemExit(f"{spec!r} asks for '@repaired' but no --alias was given")
                run_aliases = aliases
            else:
                raise SystemExit(f"unknown run modifier {token!r} in {spec!r}")
        folder = Path(run)
        if not folder.is_dir():
            folder = args.root / run
        if not folder.is_dir():
            raise SystemExit(f"run folder not found: {run}")
        note = rt_source + " RTs" + (", re-paired truth" if run_aliases else "")
        print(f"[bench] {folder.name} ({note}) ...", file=sys.stderr, flush=True)
        reports.append(bench_run(folder, params, rt_source, run_aliases))

    command = "python tools/bench_areas.py " + " ".join(
        f'"{a}"' if "@" in a or " " in a else a for a in (argv if argv is not None else sys.argv[1:])
    )
    text = render(reports, params, args.min_area, command)
    print(text)
    if str(args.out) != "-":
        out = args.out if args.out.is_absolute() else args.root / args.out
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text, encoding="utf-8")
        print(f"[bench] wrote {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
