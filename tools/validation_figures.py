"""Draw the archive validation figures from a pooled comparison.

Reads the CSVs that ``tools/master_bench.py`` writes and draws every data figure of the
validation deck as SVG: agreement by method and flag, the area and mass scatters, Aroclor
fingerprint agreement, and the breakdowns of where and why the two integrations differ.

    python tools/validation_figures.py --data docs/validation/master --out docs/images/validation

The figures name run folders and plot per-sample results, so the output is laboratory data.
Keep it out of anything published.

The Aroclor figures need the congener-level reference table (``--aroclor``). The run table
fills its file counts from the run folders under ``--root``; a run whose folder is absent
gets blank counts.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import quant209  # noqa: F401  (BLAS guards before numpy does any work)

sys.path.insert(0, str(Path(__file__).resolve().parent))

from figlib import agreement, detection, drivers, errors, fingerprint, performance  # noqa: E402
from figlib.data import load  # noqa: E402

MODULES = [performance, agreement, fingerprint, errors, detection, drivers]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, default=Path("docs/validation/master"),
                    help="folder holding areas.csv, masses.csv, recoveries.csv and runs.csv")
    ap.add_argument("--out", type=Path, default=Path("docs/images/validation"))
    ap.add_argument("--aroclor", type=Path, default=Path("rushneckAroclors_173list.csv"),
                    help="congener-level Aroclor reference table")
    ap.add_argument("--root", type=Path, default=Path("."), help="archive of run folders")
    ap.add_argument("--only", nargs="*", default=None, metavar="NUMBER",
                    help="draw only these figures, by number (for example 34 41 50a)")
    ap.add_argument("--as-first-drawn", action="store_true",
                    help="leave out the legends, axis titles and notes added since the figures were first drawn")
    args = ap.parse_args(argv)

    started = time.time()
    ds = load(args.data, aroclor=args.aroclor, root=args.root)
    ds.extras = not args.as_first_drawn
    if not ds.aroclors:
        print(f"no Aroclor reference at {args.aroclor}: the fingerprint figures are skipped")

    args.out.mkdir(parents=True, exist_ok=True)
    wanted = set(args.only) if args.only else None
    written = 0
    for module in MODULES:
        for figure in module.FIGURES:
            number = figure.__name__.removeprefix("fig_")
            if wanted is not None and number not in wanted:
                continue
            if getattr(figure, "needs_aroclors", False) and not ds.aroclors:
                continue
            name, svg = figure(ds)
            (args.out / name).write_text(svg, encoding="utf-8", newline="\n")
            written += 1
            print(f"  {name}")
    print(f"wrote {written} figures to {args.out} in {time.time() - started:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
