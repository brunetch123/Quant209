"""Re-take the screenshots in the README, against the synthetic demo run.

This drives the real application (the same widgets a user sees) rather than mocking
anything, so a screenshot cannot drift from the interface it claims to show. It needs a run
folder. ``tools/make_demo_run.py`` writes a fictional one.

    python tools/make_demo_run.py --out .
    python tools/make_screenshots.py

Writes ``docs/images/review-overview.png`` and ``docs/images/review-samples.png``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

DEMO_RUN = "DEMO01-06_20240115"

#: Window size the images are taken at: a comfortable laptop window, not a maximised 4K one.
WIDTH, HEIGHT = 2000, 1250

#: Which Review tab each image shows.
SHOTS = ((0, "review-overview"), (1, "review-samples"))

#: The Samples tab opens on the first transition, which is the three-analyte mono channel.
#: The penta channel carries 37 analytes and is what the interface is actually for.
SAMPLES_TRANSITION = "penta"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", default=DEMO_RUN, help="run folder to open")
    ap.add_argument("--out", default="docs/images", help="directory to write the PNGs to")
    args = ap.parse_args(argv)

    run_dir = ROOT / args.run
    if not (run_dir / "Converted").is_dir():
        print(f"no converted run at {run_dir}. Run tools/make_demo_run.py first")
        return 1

    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QApplication

    from quant209.app.main import MainWindow, configure_app
    from quant209.app.state import STEP_REVIEW
    from quant209.core.pipeline import run_analysis
    from quant209.core.sequence import find_data_files, scan_inputs

    QApplication.setHighDpiScaleFactorRoundingPolicy(
        Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
    )
    app = QApplication.instance() or QApplication(sys.argv[:1])
    configure_app(app)

    print(f"analysing {args.run} ...")
    entries = scan_inputs(find_data_files(run_dir))
    result = run_analysis(entries, run_dir / "Converted")
    print(f"  {result.status}: {len(result.sample_ids)} samples, {len(result.blank_ids)} blanks")

    window = MainWindow()
    window.resize(WIDTH, HEIGHT)
    window.state.set_run_folder(run_dir)
    window.state.set_entries(entries)
    window.state.set_result(result)
    window.goto_step(STEP_REVIEW)
    window.show()
    for _ in range(10):
        app.processEvents()

    review = window.page(STEP_REVIEW)
    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    for index, name in SHOTS:
        tab = review.ensure_tab(index)
        review.tabs.setCurrentIndex(index)
        for _ in range(20):
            app.processEvents()
        combo = getattr(tab, "transition_combo", None)
        if combo is not None:
            for i in range(combo.count()):
                if combo.itemText(i).startswith(SAMPLES_TRANSITION):
                    combo.setCurrentIndex(i)
                    break
            for _ in range(20):
                app.processEvents()
        path = out_dir / f"{name}.png"
        window.grab().save(str(path))
        print(f"  wrote {path.relative_to(ROOT)}")

    window.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
