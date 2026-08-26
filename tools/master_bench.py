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
from quant209.core.config import load_canonical_peaks
from quant209.core.pipeline import run_analysis
from quant209.core.sequence import find_data_files, scan_inputs
from quant209.report.truth import read_run_truth
from quant209.report.validation import join_areas, join_masses, join_recoveries

NOT_A_RUN = {
    ".git", ".venv", ".pytest_cache", "build", "docs", "quant209", "tests",
    "tools", "quant209.egg-info", "Background", "MSconvert", "__pycache__",
}

AREA_FIELDS = ["run", "sample", "file_type", "label", "homolog", "role",
               "manual", "auto", "flags", "method", "sn", "rt_apex"]
MASS_FIELDS = ["run", "sample", "file_type", "label", "homolog", "role",
               "manual", "auto", "quantified", "flags", "method", "rt_apex", "detected"]
REC_FIELDS = ["run", "sample", "homolog", "manual", "auto"]
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


def bench_one(folder: Path, analytes) -> tuple[list, list, list, dict]:
    """Integrate, quantify and join one run. Never raises: failures are reported as rows."""
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
            return [], [], [], summary

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
        areas = join_areas(peaks, types, analytes, truth.areas)
        masses, _truth_only, _file_only = join_masses(result, truth.masses)
        recs, _unmatched = join_recoveries(result, truth.recoveries)
        summary["joined_areas"] = len(areas)
        summary["joined_masses"] = len(masses)
        summary["status"] = "ok"
        summary["seconds"] = round(time.time() - started, 1)
        return areas, masses, recs, summary
    except Exception as exc:  # a bad run must not stop the archive sweep
        summary["status"] = "error"
        summary["note"] = f"{type(exc).__name__}: {exc}"[:200]
        summary["seconds"] = round(time.time() - started, 1)
        traceback.print_exc()
        return [], [], [], summary


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
    aw = csv.DictWriter(area_f, AREA_FIELDS); aw.writeheader()
    mw = csv.DictWriter(mass_f, MASS_FIELDS); mw.writeheader()
    rw = csv.DictWriter(rec_f, REC_FIELDS); rw.writeheader()
    sw = csv.DictWriter(run_f, RUN_FIELDS); sw.writeheader()

    started = time.time()
    try:
        for i, folder in enumerate(folders, start=1):
            areas, masses, recs, summary = bench_one(folder, analytes)
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
            for fh in (area_f, mass_f, rec_f, run_f):
                fh.flush()
            print(f"[{i}/{len(folders)}] {folder.name}: {summary['status']} "
                  f"areas={summary['joined_areas']} masses={summary['joined_masses']} "
                  f"{summary['seconds']}s", flush=True)
    finally:
        for fh in (area_f, mass_f, rec_f, run_f):
            fh.close()

    print(f"\nswept {len(folders)} runs in {(time.time() - started) / 60:.1f} min -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
