"""Run the pipeline over every analysis run in the archive and pool the comparison.

Produces one row per (run, file, analyte) for areas and for masses, plus per-run and
per-homolog rollups, so that failure modes can be looked for across the whole archive
rather than one run at a time. Runs are processed one at a time and written out as they
finish, so a crash or an out-of-memory kill costs one run, not the batch.

    python tools/master_bench.py                       # all runs
    python tools/master_bench.py --out docs/validation/master
    python tools/master_bench.py --runs <run-folder> <run-folder>
"""

from __future__ import annotations

import argparse
import csv
import time
import traceback
import warnings
from dataclasses import asdict
from pathlib import Path

import quant209  # noqa: F401  (BLAS guards before numpy)
import pandas as pd  # noqa: E402  (after the guards)
from quant209.core.config import load_canonical_peaks
from quant209.core.pipeline import run_analysis
from quant209.core.sequence import find_data_files, scan_inputs
from quant209.report.truth import read_run_truth
from quant209.report.validation import (
    join_areas, join_masses, join_recoveries, rebase_group_masses, verify_pairing,
)

NOT_A_RUN = {
    ".git", ".venv", ".pytest_cache", "build", "docs", "quant209", "tests",
    "tools", "quant209.egg-info", "Background", "MSconvert", "__pycache__",
}

AREA_FIELDS = ["run", "sample", "file_type", "label", "homolog", "role",
               "manual", "auto", "flags", "method", "sn", "rt_apex"]
MASS_FIELDS = ["run", "sample", "file_type", "label", "homolog", "role",
               "manual", "auto", "quantified", "flags", "method", "rt_apex", "detected"]
REC_FIELDS = ["run", "sample", "homolog", "manual", "auto"]
PAIRING_FIELDS = ["run", "sheet", "file", "own_distance", "best_other", "best_distance", "verdict"]
RRF_FIELDS = ["run", "label", "homolog", "is_norm", "ours", "analyst"]
RUN_FIELDS = ["run", "status", "files", "samples", "truth_cols", "joined_areas",
              "joined_masses", "fit_ok", "fit_reason", "seconds", "note"]


def run_folders(root: Path) -> list[Path]:
    out = []
    for d in sorted(root.iterdir()):
        if not d.is_dir() or d.name in NOT_A_RUN or d.name.startswith("."):
            continue
        if (d / "Converted").exists() or any(
            p.is_dir() and p.suffix.upper() == ".D" for p in d.iterdir()
        ):
            out.append(d)
    return out


def bench_one(folder: Path, analytes) -> tuple[list, list, list, dict, list, list]:
    """Integrate, quantify and join one run. Never raises: failures are reported as rows.

    Every analyst sheet is first checked against the data file it is about to be scored
    against (``report.validation.verify_pairing``), and a proven crossing is scored against
    the right file. Masses are not scored for a file the program correctly refuses: one whose
    internal standard is not in the injection, so that it reports an error and no mass. Both
    decisions are written to the run's note and to ``pairing.csv``.
    """
    started = time.time()
    summary = {k: "" for k in RUN_FIELDS}
    summary["run"] = folder.name
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            truth = read_run_truth(folder)
        populated = [c for c in truth.areas.columns if truth.areas[c].notna().any()]
        summary["truth_cols"] = len(populated)
        if not populated:
            summary["status"] = "no-truth"
            summary["note"] = "no analyst quantification in this folder"
            summary["seconds"] = round(time.time() - started, 1)
            return [], [], [], summary, [], []

        entries = scan_inputs(find_data_files(folder))
        summary["files"] = len(entries)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = run_analysis(entries, converted_dir=folder / "Converted")
        summary["samples"] = len(result.sample_ids)
        summary["fit_ok"] = result.fit.fit_ok
        summary["fit_reason"] = "" if result.fit.fit_ok else result.fit.reason

        peaks = {fid: rec.peaks for fid, rec in result.files.items()}
        types = {fid: rec.entry.file_type for fid, rec in result.files.items()}
        aliases, pairing = verify_pairing(result, truth.areas)
        areas = join_areas(peaks, types, analytes, truth.areas, aliases)
        truth_masses, rebased = rebase_group_masses(truth, analytes)
        masses, _truth_only, _file_only = join_masses(result, truth_masses, aliases)
        recs, _unmatched = join_recoveries(result, truth.recoveries, aliases)
        masses, unscored = _drop_refused(masses, result, analytes)
        notes = list(rebased)
        if aliases:
            notes.append("re-paired by labelled-standard fingerprint: "
                         + ", ".join(f"{f}->{s}" for f, s in sorted(aliases.items())))
        if unscored:
            notes.append("masses not scored, internal standard not in the injection: "
                         + ", ".join(sorted(unscored)))
        summary["note"] = "; ".join(notes)
        summary["joined_areas"] = len(areas)
        summary["joined_masses"] = len(masses)
        summary["status"] = "ok"
        summary["seconds"] = round(time.time() - started, 1)
        analyst_rrf = (truth.cal_rrfs.median(axis=1, skipna=True)
                       if not truth.cal_rrfs.empty else pd.Series(dtype=float))
        rrfs = [{"label": a.label, "homolog": a.homolog, "is_norm": a.is_norm,
                 "ours": result.calibration.rrf.get(a.label, ""),
                 "analyst": analyst_rrf.get(a.label, "")}
                for a in analytes if a.role == "native"]
        return areas, masses, recs, summary, pairing, rrfs
    except Exception as exc:  # a bad run must not stop the archive sweep
        summary["status"] = "error"
        summary["note"] = f"{type(exc).__name__}: {exc}"[:200]
        summary["seconds"] = round(time.time() - started, 1)
        traceback.print_exc()
        return [], [], [], summary, [], []


def _drop_refused(masses: list, result, analytes) -> tuple[list, set[str]]:
    """Mass rows minus those the program refused because a standard is not in the injection.

    Such a file gets an error and no mass by design, so there is nothing of the program's to
    score, and the analyst's number for it rests on standard areas that were not measured. A
    standard the program merely failed to find (``not_found``) stays scored: that is the
    program's failure, and the benchmark must keep counting it.
    """
    norm = {a.label: a.is_norm for a in analytes}
    absent = {
        (fid, flag.split(":", 1)[1].split("=", 1)[0])
        for fid, flags in (result.flags or {}).items()
        for flag in flags
        if flag.startswith("missing_internal_standard:") and flag.endswith("=absent")
    }
    if not absent:
        return masses, set()
    kept, dropped = [], set()
    for row in masses:
        fid = row.sample.split("->", 1)[0]
        if (fid, norm.get(row.label)) in absent:
            dropped.add(fid)
        else:
            kept.append(row)
    return kept, dropped


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--out", type=Path, default=Path("docs/validation/master"))
    ap.add_argument("--runs", nargs="*", default=None)
    args = ap.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)
    analytes = load_canonical_peaks()
    folders = run_folders(args.root)
    if args.runs:
        wanted = set(args.runs)
        folders = [f for f in folders if f.name in wanted]

    area_f = (args.out / "areas.csv").open("w", newline="", encoding="utf-8")
    mass_f = (args.out / "masses.csv").open("w", newline="", encoding="utf-8")
    rec_f = (args.out / "recoveries.csv").open("w", newline="", encoding="utf-8")
    run_f = (args.out / "runs.csv").open("w", newline="", encoding="utf-8")
    pair_f = (args.out / "pairing.csv").open("w", newline="", encoding="utf-8")
    rrf_f = (args.out / "rrfs.csv").open("w", newline="", encoding="utf-8")
    aw = csv.DictWriter(area_f, AREA_FIELDS); aw.writeheader()
    mw = csv.DictWriter(mass_f, MASS_FIELDS); mw.writeheader()
    rw = csv.DictWriter(rec_f, REC_FIELDS); rw.writeheader()
    sw = csv.DictWriter(run_f, RUN_FIELDS); sw.writeheader()
    pw = csv.DictWriter(pair_f, PAIRING_FIELDS); pw.writeheader()
    fw = csv.DictWriter(rrf_f, RRF_FIELDS); fw.writeheader()

    started = time.time()
    try:
        for i, folder in enumerate(folders, start=1):
            areas, masses, recs, summary, pairing, rrfs = bench_one(folder, analytes)
            for row in rrfs:
                fw.writerow({"run": folder.name, **row})
            for row in pairing:
                pw.writerow({"run": folder.name, **asdict(row)})
            for row in areas:
                aw.writerow({"run": folder.name, **asdict(row)})
            for row in masses:
                mw.writerow({"run": folder.name, **asdict(row)})
            for row in recs:
                d = asdict(row)
                rw.writerow({"run": folder.name, "sample": d.get("sample", ""),
                             "homolog": d.get("homolog", ""),
                             "manual": d.get("manual", ""), "auto": d.get("auto", "")})
            sw.writerow(summary)
            for fh in (area_f, mass_f, rec_f, run_f, pair_f, rrf_f):
                fh.flush()
            print(f"[{i}/{len(folders)}] {folder.name}: {summary['status']} "
                  f"areas={summary['joined_areas']} masses={summary['joined_masses']} "
                  f"{summary['seconds']}s", flush=True)
    finally:
        for fh in (area_f, mass_f, rec_f, run_f, pair_f, rrf_f):
            fh.close()

    print(f"\nswept {len(folders)} runs in {(time.time() - started) / 60:.1f} min -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
