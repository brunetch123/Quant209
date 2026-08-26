"""Measure a run's processing time, and whether the window stays alive while it happens.

The two gates
--------------
**Wall time** from "Run analysis" to a result on the state, through the real GUI path: a
shown :class:`~quant209.app.main.MainWindow`, a real
:class:`~quant209.app.workers.PipelineWorker` on the global thread pool, real progress
signals driving the real Process page. Not ``run_analysis`` called directly: the question
is what the analyst waits for, and that includes the marshalling and the repaints.

**UI-thread responsiveness**, which is the gate that matters more. An application that
stops repainting for a second while it works is one the operating system offers to close
for you. The budget is **no stall over 100 ms**: below about that, a window still feels
attached to the mouse.

How a stall is measured
------------------------
A :class:`QTimer` at :data:`HEARTBEAT_MS` with a :class:`QElapsedTimer` restarted on every
tick. The gap it reports is how long the UI thread went without servicing its event loop.
That *is* what "the window froze" means, and it is not the same as how long the worker took.
A pool thread doing arithmetic does not stall the UI. A pool thread holding the GIL through
a long call does, and only the heartbeat can tell the two apart.

``Qt.PreciseTimer`` is used because Windows' default coarse timer has ~15 ms granularity,
which is most of the budget and would report quantisation as jitter.

Reading the output
-------------------
``median`` should sit at the nominal interval, and its distance from ``max`` is the whole
story. A handful of gaps in the 20-50 ms range is ordinary scheduling on a loaded desktop.
Anything over 100 ms is a real freeze and is counted separately.

``--runs N`` processes the run N times **in one process**, which is a different and harsher
measurement than N fresh launches: this host has ~2 GB of RAM, and consecutive runs leave it
under enough memory pressure to page. One run per process is what an analyst does. Several
in one process is where an accumulation bug would show up. Both matter, so both are
runnable.

The navigation that follows is measured too, separately and out of budget: opening Review
builds a page, and that cost belongs to a click rather than to processing.
"""

from __future__ import annotations

import statistics
import sys
import time
from pathlib import Path

from PySide6.QtCore import QElapsedTimer, Qt, QTimer
from PySide6.QtWidgets import QApplication

from quant209.app.main import MainWindow, configure_app
from quant209.app.state import STEP_PROCESS, STEP_REVIEW
from quant209.core.sequence import find_data_files, scan_inputs

#: Heartbeat interval. Ten milliseconds samples a 100 ms budget ten times over.
HEARTBEAT_MS = 10

#: A gap longer than this is a visible freeze, and the gate.
STALL_MS = 100.0

def default_run() -> Path | None:
    """The first run folder found beside the repository, when none is given on the command line."""
    root = Path(__file__).resolve().parents[1]
    for d in sorted(root.iterdir()):
        if d.is_dir() and (d / "Converted").is_dir():
            return d
    return None


def summarise_gaps(label: str, gaps: list[float], wall: float | None = None) -> list[float]:
    """Print one line of heartbeat statistics. Return the gaps that broke the budget."""
    over = [g for g in gaps if g > STALL_MS]
    line = f"{label:24s}"
    if wall is not None:
        line += f" wall {wall:6.2f} s "
    if gaps:
        ordered = sorted(gaps)
        line += (
            f" ticks {len(gaps):4d}  median {statistics.median(gaps):5.1f}"
            f"  p99 {ordered[int(0.99 * (len(ordered) - 1))]:5.1f}"
            f"  max {max(gaps):6.1f} ms  >{STALL_MS:.0f}ms: {len(over)}"
        )
    print(line)
    return over


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    runs = 1
    if "--runs" in args:
        runs = int(args[args.index("--runs") + 1])
    folder = Path(args[0]) if args and not args[0].startswith("--") else default_run()

    app = QApplication([])
    configure_app(app)
    window = MainWindow()
    window.resize(1600, 1000)
    window.show()
    app.processEvents()

    state = window.state
    state.set_run_folder(folder)
    files = find_data_files(folder)
    if not files:
        print(f"no data files in {folder}")
        return 1
    state.set_entries(scan_inputs(files))
    if not window.goto_step(STEP_PROCESS):
        print("the Process step is locked. Are the mzML converted?")
        return 1
    page = window.page(STEP_PROCESS)
    app.processEvents()

    print(f"{folder.name}: {len(files)} files, {runs} run(s) in this process")
    print(f"heartbeat {HEARTBEAT_MS} ms   stall gate {STALL_MS:.0f} ms\n")

    gaps: list[float] = []
    beat = QElapsedTimer()
    heartbeat = QTimer()
    heartbeat.setInterval(HEARTBEAT_MS)
    heartbeat.setTimerType(Qt.TimerType.PreciseTimer)
    heartbeat.timeout.connect(lambda: gaps.append(float(beat.restart())))

    walls: list[float] = []
    stalls: list[float] = []
    progress = {"i": 0, "t0": 0.0}

    def start_run() -> None:
        gaps.clear()
        beat.start()
        heartbeat.start()
        progress["t0"] = time.perf_counter()
        page.run_analysis()

    def on_result(result: object) -> None:
        if result is None:
            return
        wall = time.perf_counter() - progress["t0"]
        heartbeat.stop()
        walls.append(wall)
        stalls.extend(summarise_gaps(f"  process run {progress['i'] + 1}", list(gaps), wall))
        progress["i"] += 1
        if progress["i"] < runs:
            state.set_result(None)
            QTimer.singleShot(200, start_run)
        else:
            QTimer.singleShot(200, navigate)

    state.result_changed.connect(on_result)

    def navigate() -> None:
        """Opening Review: a page build, out of budget but worth knowing."""
        gaps.clear()
        beat.start()
        heartbeat.start()
        started = time.perf_counter()
        window.goto_step(STEP_REVIEW)
        app.processEvents()
        wall = time.perf_counter() - started
        heartbeat.stop()
        print()
        summarise_gaps("  open Review", list(gaps), wall)
        QTimer.singleShot(0, app.quit)

    QTimer.singleShot(400, start_run)
    app.exec()

    if not walls:
        print("\nno result. The run did not finish")
        return 1
    print(
        f"\nprocessing wall: median {statistics.median(walls):.2f} s "
        f"(min {min(walls):.2f}, max {max(walls):.2f}) over {len(walls)} run(s)"
    )
    passed = not stalls
    print(
        f"UI-thread stalls > {STALL_MS:.0f} ms during processing: {len(stalls)}  "
        f"({'PASS' if passed else 'FAIL'})"
    )
    if stalls:
        print(f"  worst: {['%.0f ms' % s for s in sorted(stalls, reverse=True)[:5]]}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
