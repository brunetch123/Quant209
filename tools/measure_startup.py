"""Measure Quant 209's cold start against the product's 3-second budget.

The gate
---------
**Median wall-clock time from process launch to a painted window must be under 3000 ms.**

Three seconds is roughly the point at which a desktop application stops feeling like it
opened and starts feeling like it is loading.

What is measured
-----------------
Each run launches ``quant209-gui --measure-startup`` as a subprocess. That flag makes the
application quit as soon as its first frame has been painted (see
``quant209.app.main.MainWindow.paintEvent``), so the process's own lifetime *is* the
cold-start time. Two numbers come out of each run:

``wall``
    Measured here, by this script, around the whole subprocess. **This is the gate.** It
    includes the interpreter start, ``site`` processing, and the Qt plugin load: every
    millisecond between the user's double-click and the window, including the ones the
    application itself cannot see.
``in-process``
    Reported by the application, from the moment ``quant209.app.main`` began importing to
    the first painted frame. Always smaller than ``wall``. The difference is interpreter
    startup. Useful for attributing a regression: if ``wall`` grew and ``in-process`` did
    not, something was added to the import path outside this package.

"Cold-ish", honestly
---------------------
A genuinely cold start (nothing of Python, Qt or this package in the OS file cache) cannot
be reproduced from a script without clearing the system cache, which needs administrator
rights. The first run below is the coldest of the three and is usually a few hundred
milliseconds slower than the rest. The median is the steady-state number and is what the
gate is applied to. If the *first* run is over budget while the median is not, the machine's
disk is the story, not this application's imports.

Telling a slow application from a slow machine (``--calibrate``)
----------------------------------------------------------------
Most of this number is not the application's own cost. An empty ``QMainWindow``
(interpreter, ``site``, PySide6, the Qt platform plugin, one window) is the floor no
application can get under, and on a loaded or memory-pressured host that floor moves a
*lot*: measured on the development host at 591 ms early in a session and 1822 ms two hours
later, with no code change between them. An absolute number read on its own will therefore
accuse the application of a regression it did not have.

``--calibrate`` interleaves runs of an empty Qt window with runs of the real one and reports
both medians and the difference. **The difference is the application's own cost and is the
number to watch across commits.** Roughly 1 s of it is building and first-painting the
window, which is normal for a Qt application of this size. If the gate fails while the
difference is unchanged, the machine is the story: re-measure on an idle host before
touching the code.

When it fails
--------------
Almost always an eager import. Check the import graph first::

    .venv\\Scripts\\python.exe -X importtime -c "import quant209.app.main" 2>&1 | sort -k2 -n -r | head

``quant209.app.main`` must reach a painted window without numpy, pandas, pyqtgraph,
openpyxl or ``core.pipeline``. ``tests/test_app_shell.py`` asserts exactly that, and it is
the assertion that usually catches the regression before this script does.
"""

from __future__ import annotations

import re
import statistics
import subprocess
import sys
import time
from pathlib import Path

#: The product budget, in milliseconds. Applied to the median wall-clock time.
GATE_MS = 3000.0

RUNS = 3

_STARTUP_RE = re.compile(r"startup_ms=([0-9.]+)")


def gui_command() -> list[str]:
    """The command that launches the GUI, preferring the installed console script.

    The ``quant209-gui`` executable is what a user actually runs, so it is what should be
    timed. Its launcher stub is a real (small) part of the cost. ``python -m
    quant209.app.main`` is the fallback for a source checkout that was never installed.
    """
    scripts = Path(sys.executable).parent
    for name in ("quant209-gui.exe", "quant209-gui"):
        candidate = scripts / name
        if candidate.exists():
            return [str(candidate)]
    return [sys.executable, "-m", "quant209.app.main"]


def measure_once(command: list[str]) -> tuple[float, float | None]:
    """Launch the GUI with ``--measure-startup``. Return ``(wall_ms, in_process_ms)``."""
    start = time.perf_counter()
    proc = subprocess.run(
        [*command, "--measure-startup"], capture_output=True, text=True, timeout=120
    )
    wall_ms = (time.perf_counter() - start) * 1000.0
    if proc.returncode != 0:
        raise RuntimeError(
            f"{command[0]} exited {proc.returncode}\n{proc.stderr.strip() or proc.stdout.strip()}"
        )
    match = _STARTUP_RE.search(proc.stdout or "")
    return wall_ms, float(match.group(1)) if match else None


#: An empty Qt window: the floor no application can get under on this machine, right now.
_BARE_QT = (
    "from PySide6.QtWidgets import QApplication, QMainWindow;"
    "app=QApplication([]);w=QMainWindow();w.show();app.processEvents()"
)


def measure_bare_qt() -> float:
    """Wall-clock milliseconds for an empty ``QMainWindow``: interpreter and Qt only."""
    start = time.perf_counter()
    subprocess.run([sys.executable, "-c", _BARE_QT], capture_output=True, timeout=120)
    return (time.perf_counter() - start) * 1000.0


def main() -> int:
    calibrate = "--calibrate" in sys.argv[1:]
    command = gui_command()
    print(f"launching: {' '.join(command)} --measure-startup")
    print(f"runs: {RUNS}   gate: median wall < {GATE_MS:.0f} ms\n")

    walls: list[float] = []
    floors: list[float] = []
    for run in range(1, RUNS + 1):
        # Interleaved, not measured in a block: the host drifts over the length of a run,
        # and a floor measured before the application runs would be a different machine's.
        if calibrate:
            floors.append(measure_bare_qt())
        wall_ms, in_process_ms = measure_once(command)
        walls.append(wall_ms)
        reported = f"{in_process_ms:8.1f} ms" if in_process_ms is not None else "       n/a"
        line = f"  run {run}   wall {wall_ms:8.1f} ms   in-process {reported}"
        if calibrate:
            line += f"   bare-Qt floor {floors[-1]:8.1f} ms"
        print(line)

    median = statistics.median(walls)
    if calibrate:
        floor = statistics.median(floors)
        print(
            f"\nbare-Qt floor: {floor:.1f} ms"
            f"   Quant 209's own cost: {median - floor:.1f} ms"
            "   <- this is the number to watch across commits"
        )
    passed = median < GATE_MS
    print(f"\nmedian wall: {median:.1f} ms   ({'PASS' if passed else 'FAIL'}, gate {GATE_MS:.0f} ms)")
    if not passed:
        print(
            "\nOver budget. Check the import graph first. It is almost always an eager "
            'import:\n  python -X importtime -c "import quant209.app.main"'
        )
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
