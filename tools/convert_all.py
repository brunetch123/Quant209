"""Convert every quantifiable .D acquisition in a directory of run folders to mzML.

Hex-blank injections are skipped: they are inter-sample rinses that no calculation
consumes (``core.pipeline`` never opens one), and they are 44% of the acquisitions
in this archive, so converting them would nearly double the wall time for nothing.

Existing outputs are reused unless they are truncated (an interrupted MSConvert
leaves a header-only file, which ``convert.msconvert.is_truncated`` detects).

    python tools/convert_all.py                # convert everything under the repo root
    python tools/convert_all.py --workers 3    # msconvert runs are short. A few in flight is fine
    python tools/convert_all.py --dry-run      # just report the plan
"""

from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import quant209  # noqa: F401  (sets the BLAS guards before numpy is imported)
from quant209.convert.msconvert import convert_file, is_truncated
from quant209.core.sequence import detect_type

#: Repository infrastructure that is never an analysis run.
NOT_A_RUN = {
    ".git", ".venv", ".pytest_cache", "build", "docs", "quant209", "tests",
    "tools", "quant209.egg-info", "Background", "MSconvert", "__pycache__",
}

#: File types the quantification chain never opens (see ``core.pipeline.QUANTIFIED_TYPES``).
SKIP_TYPES = {"hex"}


def run_folders(root: Path) -> list[Path]:
    """Every immediate subdirectory that holds at least one ``.D`` acquisition."""
    out = []
    for d in sorted(root.iterdir()):
        if not d.is_dir() or d.name in NOT_A_RUN or d.name.startswith("."):
            continue
        if any(p.is_dir() and p.suffix.upper() == ".D" for p in d.iterdir()):
            out.append(d)
    return out


def pending(folder: Path) -> list[Path]:
    """The ``.D`` files in ``folder`` that still need converting."""
    converted = folder / "Converted"
    todo = []
    for src in sorted(folder.iterdir()):
        if not (src.is_dir() and src.suffix.upper() == ".D"):
            continue
        if detect_type(src.stem) in SKIP_TYPES:
            continue
        out = converted / f"{src.stem}.mzML"
        if out.exists() and not is_truncated(out):
            continue
        todo.append(src)
    return todo


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--msconvert", type=Path, default=Path("MSconvert/pwiz-bin/msconvert.exe"))
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    if not args.msconvert.exists():
        print(f"msconvert not found: {args.msconvert}", file=sys.stderr)
        return 1

    jobs: list[tuple[Path, Path]] = []
    for folder in run_folders(args.root):
        todo = pending(folder)
        if todo:
            print(f"{folder.name}: {len(todo)} to convert", flush=True)
        jobs.extend((src, folder / "Converted") for src in todo)

    print(f"\n{len(jobs)} acquisitions to convert", flush=True)
    if args.dry_run or not jobs:
        return 0

    started = time.time()
    done = failures = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(convert_file, src, args.msconvert, out): src for src, out in jobs
        }
        for fut in as_completed(futures):
            src = futures[fut]
            done += 1
            try:
                result = fut.result()
            except Exception as exc:  # a crashed msconvert must not stop the batch
                failures += 1
                print(f"  !! {src.name}: {exc!r}", flush=True)
                continue
            if not result.ok:
                failures += 1
                print(f"  !! {src.name}: {result.message}", flush=True)
            if done % 25 == 0 or done == len(jobs):
                rate = done / max(time.time() - started, 1e-9)
                left = (len(jobs) - done) / rate if rate else 0
                print(
                    f"  {done}/{len(jobs)} converted, {failures} failed, "
                    f"{rate * 60:.0f}/min, ~{left / 60:.0f} min left",
                    flush=True,
                )

    print(f"\ndone: {done - failures} converted, {failures} failed, "
          f"{(time.time() - started) / 60:.1f} min", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
