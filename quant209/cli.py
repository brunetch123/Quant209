"""quant209 command-line interface.

Headless entry point for the pipeline: scan a run folder, optionally convert raw Agilent
``.D`` acquisitions to ``.mzML`` first, run the quantification pipeline, write the
analyst-format Excel report, and (optionally) the validation report against the analyst's
ground truth. Also usable stand-alone as a conversion-only tool (``quant209 convert``).

Import discipline
------------------
This module is kept import-cheap at module scope on purpose. ``quant209.core.pipeline``
(numpy/pandas), ``quant209.core.qaqc`` and ``quant209.report.*`` (openpyxl/matplotlib) are
imported **inside** each subcommand handler, not here, so ``quant209 --help`` and argument
parsing stay instant. The PySide6 app's own cold-start budget depends on the same discipline.
``quant209.core.sequence`` and ``quant209.convert.msconvert`` are cheap (no numpy/pandas) and
imported at module scope.

Two different monkeypatch targets
---------------------------------
Because the heavy dependencies above are imported *inside* the functions that use them (a
fresh ``from module import name`` on every call), a test that wants to stub them must patch
the **source module's** attribute (e.g. ``quant209.core.pipeline.run_analysis``). Patching
``quant209.cli.run_analysis`` would have no effect, since no such module-level name exists
here. The lightweight names imported at module scope (``scan_inputs``, ``plan_conversion``,
...) are the opposite: a test patches ``quant209.cli.<name>`` directly. See
``tests/test_cli.py`` for both patterns in use.

Error handling and exit codes
-------------------------------
``main()`` never lets an exception escape as a traceback (unless ``--debug`` is passed). Any
``Exception`` raised by a subcommand handler is caught, printed as one line
(``error: ...``), and mapped to exit code 1. This covers a bad run folder, a run with no
calibration injection (``core.pipeline.run_analysis``'s own ``ValueError``), and a failed
conversion alike. ``UsageError`` is the one exception mapped to exit code 2 instead (a bad
combination of flags this module catches itself, e.g. ``--from-raw`` without
``--msconvert``). An unrecognized subcommand or missing required argument is caught the same
way via argparse's own ``SystemExit``, which ``main()`` translates into a return value so
callers (and tests) never have to catch ``SystemExit`` themselves.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import quant209  # noqa: F401 (sets BLAS thread-limit env vars before anything numpy-adjacent)

from quant209.convert.msconvert import ConversionPlan, convert_many, plan_conversion
from quant209.core.sequence import find_data_files, scan_inputs

EXIT_OK = 0
EXIT_FAILURE = 1
EXIT_USAGE = 2

#: Where ``quant209 run``/``quant209 validate`` write their output when ``--out`` is not
#: given, relative to the current working directory.
DEFAULT_OUT_DIR = "output"


class UsageError(Exception):
    """A bad combination of CLI arguments that argparse itself can't express as a
    per-option constraint (e.g. ``--from-raw`` without ``--msconvert``). Mapped to exit
    code 2 by :func:`main`, same family as argparse's own usage errors."""


# --------------------------------------------------------------------------------------
# Small formatting/progress helpers
# --------------------------------------------------------------------------------------


def _pct(x: float) -> str:
    """``0.031 -> "3.1%"``. NaN (compared against itself, the branchless NaN test) -> "n/a"."""
    return "n/a" if x != x else f"{100 * x:.1f}%"


def _pipeline_progress():
    """``progress(phase, i, n, label)`` -> one plain, log-friendly line per call.

    Plain lines rather than a carriage-return-overwritten one: this callback's output is as
    likely to land in a redirected log file as a terminal, and a plain line is unambiguous
    in both. ``i`` is 0-based per ``core.pipeline``'s own contract. It is printed 1-based.
    """

    def progress(phase: str, i: int, n: int, label: str) -> None:
        print(f"[{phase}] {i + 1}/{max(n, 1)} {label}")

    return progress


def _convert_progress():
    """``progress(i, n, name)`` (``convert.msconvert.convert_many``'s own contract)."""

    def progress(i: int, n: int, name: str) -> None:
        print(f"[convert] {i + 1}/{max(n, 1)} {name}")

    return progress


def _require_folder(folder: str | Path) -> Path:
    folder = Path(folder)
    if not folder.is_dir():
        raise ValueError(f"not a run folder (not a directory): {folder}")
    return folder


# --------------------------------------------------------------------------------------
# Conversion, shared between `run --from-raw` and `convert`
# --------------------------------------------------------------------------------------


def _convert_folder(folder: Path, msconvert_exe: Path, out_dir: Path) -> ConversionPlan:
    """Plan and run MSConvert over ``folder``'s raw ``.D`` acquisitions into ``out_dir``.

    Prints the plan, per-file progress, and a results summary. Raises ``ValueError`` (a
    conversion error, exit code 1, see the module docstring) if any pending conversion
    failed. The per-file failure lines are already on stdout by the time it raises, so
    callers should not repeat them.
    """
    d_files = find_data_files(folder, prefer_mzml=False)
    entries = scan_inputs(d_files)
    plan = plan_conversion(entries, out_dir)
    print(
        f"[convert] plan: {len(plan.to_convert)} to convert, "
        f"{len(plan.reconvert)} to reconvert (truncated), {len(plan.skip)} already converted"
    )
    pending = plan.pending
    if not pending:
        print("[convert] nothing to do")
        return plan

    results = convert_many(
        [e.path for e in pending], msconvert_exe, out_dir, progress=_convert_progress()
    )
    failures = [r for r in results if not r.ok]
    print(f"[convert] {len(results) - len(failures)}/{len(results)} converted successfully")
    for r in failures:
        print(f"[convert]   FAILED {r.source.name}: {r.message}")
    if failures:
        raise ValueError(
            f"{len(failures)} of {len(results)} conversion(s) failed: "
            + ", ".join(r.source.name for r in failures)
        )
    return plan


# --------------------------------------------------------------------------------------
# Validation, shared between `run --validate` and `validate`
# --------------------------------------------------------------------------------------


def _print_validation_headline(report, out_dir: Path) -> None:
    print(
        f"[validate] area    median |RPD| {_pct(report.area.median_abs_rpd)}"
        f"  within +/-25% {_pct(report.area.within_band)}"
    )
    print(
        f"[validate] mass    median |RPD| {_pct(report.mass.median_abs_rpd)}"
        f"  within +/-25% {_pct(report.mass.within_band)}"
    )
    mae = report.recovery.mae
    print(f"[validate] recovery MAE {'n/a' if mae != mae else f'{mae:.4f}'}")
    print(f"[validate] verdict: {'PASS' if report.passed else 'FAIL'}")
    print(f"[validate] report written to: {out_dir}")


def _run_validation(folder: Path, result, out_dir: Path, *, required: bool) -> None:
    """Read ground truth for ``folder`` and, if present, write the validation report.

    ``required=True`` (the ``validate`` subcommand) raises ``ValueError`` (exit code 1) when
    no analyst workbook is found. ``required=False`` (``run --validate``) treats it as a note
    and skips the report, since ``--validate`` on ``run`` is opportunistic rather than the
    whole point of the command.
    """
    from quant209.report.truth import read_run_truth
    from quant209.report.validation import generate_validation_report

    truth = read_run_truth(folder)
    if truth.areas.empty:
        message = f"no analyst quantification workbook found in {folder}"
        if required:
            raise ValueError(message)
        print(f"[validate] {message}; skipping validation report")
        return

    val_out_dir = out_dir / "validation" / folder.name
    report = generate_validation_report(result, truth, val_out_dir, folder.name)
    _print_validation_headline(report, val_out_dir)


# --------------------------------------------------------------------------------------
# `run`
# --------------------------------------------------------------------------------------


def _print_run_headline(result, report_path: Path) -> None:
    print(f"[run] samples: {len(result.sample_ids)}  blanks: {len(result.blank_ids)}")
    status = result.status
    if result.status_reason:
        status = f"{status} ({result.status_reason})"
    print(f"[run] fit: {status}")

    rec = result.recoveries
    if rec.size and bool(rec.notna().any().any()):
        lo, hi = float(rec.min().min()), float(rec.max().max())
        print(f"[run] recovery range: {lo:.0%}-{hi:.0%}")
    else:
        print("[run] recovery range: n/a")

    n_flags = sum(len(v) for v in result.flags.values())
    print(f"[run] flags: {n_flags}")
    print(f"[run] report written to: {report_path}")


def cmd_run(args: argparse.Namespace) -> int:
    folder = _require_folder(args.folder)

    if args.from_raw:
        if args.msconvert is None:
            raise UsageError("--from-raw requires --msconvert EXE")
        try:
            _convert_folder(folder, args.msconvert, folder / "Converted")
        except ValueError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return EXIT_FAILURE

    entries = scan_inputs(find_data_files(folder))
    if not entries:
        raise ValueError(f"no data files found in {folder}")

    from quant209.core.pipeline import run_analysis
    from quant209.core.qaqc import QaqcParams

    qaqc = QaqcParams(k_sd=args.k) if args.k is not None else QaqcParams()
    result = run_analysis(
        entries, folder / "Converted", progress=_pipeline_progress(), qaqc=qaqc
    )

    from quant209.report.excel import write_report

    out_dir = Path(args.out) if args.out else Path(DEFAULT_OUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / f"{folder.name}_report.xlsx"
    write_report(result, report_path, run_name=folder.name)

    _print_run_headline(result, report_path)

    if args.validate:
        _run_validation(folder, result, out_dir, required=False)

    return EXIT_OK


# --------------------------------------------------------------------------------------
# `convert`
# --------------------------------------------------------------------------------------


def cmd_convert(args: argparse.Namespace) -> int:
    folder = _require_folder(args.folder)
    out_dir = Path(args.out) if args.out else folder / "Converted"
    try:
        _convert_folder(folder, args.msconvert, out_dir)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_FAILURE
    return EXIT_OK


# --------------------------------------------------------------------------------------
# `validate`
# --------------------------------------------------------------------------------------


def cmd_validate(args: argparse.Namespace) -> int:
    folder = _require_folder(args.folder)
    entries = scan_inputs(find_data_files(folder))
    if not entries:
        raise ValueError(f"no data files found in {folder}")

    from quant209.core.pipeline import run_analysis

    result = run_analysis(entries, folder / "Converted", progress=_pipeline_progress())

    out_dir = Path(args.out) if args.out else Path(DEFAULT_OUT_DIR)
    _run_validation(folder, result, out_dir, required=True)
    return EXIT_OK


# --------------------------------------------------------------------------------------
# argparse surface
# --------------------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="quant209",
        description="Headless PCB GC-MS/MS quantification pipeline.",
    )
    parser.add_argument(
        "--debug", action="store_true",
        help="on error, raise the full traceback instead of printing one line",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser(
        "run", help="scan (and optionally convert) a run folder, quantify it, write the Excel report"
    )
    p_run.add_argument("folder", type=Path, help="analysis-run folder")
    p_run.add_argument(
        "--out", type=Path, default=None,
        help=f"output directory for the report (default: {DEFAULT_OUT_DIR}/)",
    )
    p_run.add_argument(
        "--from-raw", action="store_true",
        help="convert raw .D acquisitions to mzML (into <folder>/Converted) before analyzing",
    )
    p_run.add_argument("--msconvert", type=Path, default=None, help="path to msconvert.exe (with --from-raw)")
    p_run.add_argument(
        "--validate", action="store_true",
        help="also write the validation report if an analyst quantification workbook is present",
    )
    p_run.add_argument("--k", type=float, default=None, help="LOQ = mean + k*SD of blanks (default 3.0)")
    p_run.set_defaults(func=cmd_run)

    p_convert = sub.add_parser("convert", help="convert a run folder's raw .D acquisitions to mzML")
    p_convert.add_argument("folder", type=Path, help="analysis-run folder")
    p_convert.add_argument("--msconvert", type=Path, required=True, help="path to msconvert.exe")
    p_convert.add_argument(
        "--out", type=Path, default=None,
        help="output directory for the converted mzML files (default: <folder>/Converted)",
    )
    p_convert.set_defaults(func=cmd_convert)

    p_validate = sub.add_parser(
        "validate", help="quantify a run folder and write the validation report against ground truth"
    )
    p_validate.add_argument("folder", type=Path, help="analysis-run folder")
    p_validate.add_argument(
        "--out", type=Path, default=None,
        help=f"output directory for the validation report (default: {DEFAULT_OUT_DIR}/)",
    )
    p_validate.set_defaults(func=cmd_validate)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        code = exc.code
        return code if isinstance(code, int) else EXIT_USAGE

    try:
        return args.func(args)
    except UsageError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE
    except Exception as exc:  # noqa: BLE001 (converted to a clean CLI error unless --debug)
        if args.debug:
            raise
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_FAILURE


if __name__ == "__main__":
    sys.exit(main())
