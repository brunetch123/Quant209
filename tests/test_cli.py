"""Tests for quant209.cli: the headless command-line interface.

Wiring tests, not integration tests: ``run_analysis``, ``write_report`` and
``generate_validation_report`` are monkeypatched to fast stubs for every test except the
"no calibration injection" case, which deliberately uses the real ``run_analysis`` (see that
test's docstring) since the failure it exercises happens before any file is read at all.

Two different monkeypatch targets are used, on purpose, matching how ``quant209.cli`` itself
imports things (see its module docstring):

* ``quant209.core.pipeline.run_analysis``, ``quant209.report.excel.write_report``,
  ``quant209.report.validation.generate_validation_report`` and
  ``quant209.report.truth.read_run_truth`` are imported *inside* the cli.py functions that
  use them (a fresh import on every call), so tests patch the **source module's** attribute.
* ``plan_conversion`` / ``convert_many`` are imported at cli.py's module top level, so tests
  patch ``quant209.cli.<name>`` directly.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import quant209.cli as cli
from quant209.convert.msconvert import ConversionPlan, ConversionResult
from quant209.core.sequence import FileEntry

# --------------------------------------------------------------------------------------
# A minimal, hand-built RunResult/ValidationReport pair the stubs hand back. It is enough
# for cli.py's headline-printing code (sample_ids/blank_ids/status/recoveries/flags) to run
# for real, so these tests also check the printed summary, not just that functions were
# called.
# --------------------------------------------------------------------------------------


def _stub_run_result():
    from quant209.core.calibrate import RtAlignment, RunFitReport
    from quant209.core.config import Analyte
    from quant209.core.peaks import IntegrationParams
    from quant209.core.pipeline import STATUS_OK, RunResult
    from quant209.core.qaqc import QaqcParams
    from quant209.core.quantify import CalibrationModel

    analytes = [
        Analyte(label="1", q1=188.0, q3=153.1, homolog="mono", order=1, role="native",
                is_norm="30D", cal_mass=29.06, cal_rt=19.563),
    ]
    return RunResult(
        entries=[],
        quantified_entries=[],
        analytes=analytes,
        params=IntegrationParams(),
        qaqc=QaqcParams(),
        cal_files=["Cal01"],
        ref_files=["Ref01"],
        alignment=RtAlignment(),
        anchors={"1": 19.563},
        fit=RunFitReport(fit_ok=True, native_found_frac=1.0, labeled_found_frac=1.0,
                          n_native=1, n_native_found=1, n_labeled=0, n_labeled_found=0),
        calibration=CalibrationModel(rrf={"1": 2.0}, cal_area={"1": 1e5}, is_response={},
                                      invalid={}, injections=("Cal01",)),
        status=STATUS_OK,
        sample_ids=["S1"],
        blank_ids=["LB1"],
        cal_ids=["Cal01"],
        ref_ids=["Ref01"],
        recoveries=pd.DataFrame({"S1": {"mono": 0.85}, "LB1": {"mono": 0.90}}),
        flags={"S1": [], "LB1": ["low_recovery:mono=90%"], "Cal01": [], "Ref01": []},
    )


def _stub_validation_report(out_dir: Path):
    from quant209.report.validation import Metrics, RecoveryMetrics, RrfMetrics, ValidationReport

    area = Metrics(n_joined=4, n_detectable=4, median_abs_rpd=0.02, within_band=1.0,
                    total_ratio=0.99, median_rpd=0.01)
    mass = Metrics(n_joined=4, n_detectable=4, median_abs_rpd=0.03, within_band=1.0,
                    total_ratio=0.98, median_rpd=0.01)
    recovery = RecoveryMetrics(n=2, mae=0.016, median_error=0.01, max_abs_error=0.02, within_005=1.0)
    rrf = RrfMetrics(n=1, median_abs_rpd=0.01, within_band=1.0, n_off=0, n_high=0, n_low=0,
                      median_abs_area_rpd=0.01)
    return ValidationReport(
        run="TESTRUN", out_dir=out_dir, area=area, mass=mass, recovery=recovery, rrf=rrf,
        n_bias_rows=0, area_pass=True, mass_pass=True, recovery_pass=True, passed=True,
    )


@pytest.fixture
def run_folder(tmp_path: Path) -> Path:
    """A run folder with one converted, non-cal mzML, enough for scan_inputs/find_data_files
    to find something, without needing a real MSConvert-produced file (the stubs never read
    file contents)."""
    folder = tmp_path / "TESTRUN"
    converted = folder / "Converted"
    converted.mkdir(parents=True)
    (converted / "AN_S1_20240101.mzML").write_bytes(b"stub")
    return folder


@pytest.fixture
def run_folder_with_cal(run_folder: Path) -> Path:
    (run_folder / "Converted" / "AN_Cal01_20240101.mzML").write_bytes(b"stub")
    return run_folder


# --------------------------------------------------------------------------------------
# argparse surface
# --------------------------------------------------------------------------------------


def test_unknown_command_exits_2(capsys):
    assert cli.main(["frobnicate"]) == cli.EXIT_USAGE


def test_no_command_exits_2():
    assert cli.main([]) == cli.EXIT_USAGE


def test_help_does_not_raise_and_exits_0():
    with pytest.raises(SystemExit) as excinfo:
        cli.build_parser().parse_args(["--help"])
    assert excinfo.value.code == 0


def test_missing_folder_argument_exits_2():
    assert cli.main(["run"]) == cli.EXIT_USAGE


def test_convert_without_required_msconvert_exits_2(tmp_path: Path):
    assert cli.main(["convert", str(tmp_path)]) == cli.EXIT_USAGE


# --------------------------------------------------------------------------------------
# `run`: folder validation
# --------------------------------------------------------------------------------------


def test_run_on_nonexistent_folder_exits_1_with_message(tmp_path: Path, capsys):
    missing = tmp_path / "does_not_exist"
    code = cli.main(["run", str(missing)])
    assert code == cli.EXIT_FAILURE
    err = capsys.readouterr().err
    assert "not a" in err.lower() or "not found" in err.lower()
    assert str(missing) in err


def test_run_on_folder_without_cal_exits_1_with_message(tmp_path: Path, capsys):
    """Uses the REAL run_analysis (not monkeypatched): core.pipeline.run_analysis raises
    ValueError for a missing calibration injection before it reads a single file, so a
    folder holding one non-cal, content-free stub mzML is enough to exercise this path
    end-to-end through the real pipeline wiring."""
    folder = tmp_path / "NOCAL"
    converted = folder / "Converted"
    converted.mkdir(parents=True)
    (converted / "AN_S1_20240101.mzML").write_bytes(b"not real mzml content")

    code = cli.main(["run", str(folder), "--out", str(tmp_path / "out")])

    assert code == cli.EXIT_FAILURE
    err = capsys.readouterr().err
    assert "calibration" in err.lower()


def test_run_with_no_data_files_exits_1(tmp_path: Path, capsys):
    folder = tmp_path / "EMPTY"
    folder.mkdir()
    code = cli.main(["run", str(folder)])
    assert code == cli.EXIT_FAILURE
    assert "no data files" in capsys.readouterr().err.lower()


# --------------------------------------------------------------------------------------
# `run`: wiring (stubbed pipeline)
# --------------------------------------------------------------------------------------


def test_run_writes_report_to_default_out_dir(run_folder_with_cal: Path, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    captured = {}

    def fake_run_analysis(entries, converted_dir, progress=None, qaqc=None, **kwargs):
        captured["entries_n"] = len(entries)
        captured["converted_dir"] = converted_dir
        captured["qaqc"] = qaqc
        if progress:
            progress("calibrate", 0, 1, "align")
        return _stub_run_result()

    def fake_write_report(result, out_path, run_name=None):
        captured["report_out_path"] = out_path
        captured["run_name"] = run_name
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_bytes(b"stub xlsx")
        return out_path

    monkeypatch.setattr("quant209.core.pipeline.run_analysis", fake_run_analysis)
    monkeypatch.setattr("quant209.report.excel.write_report", fake_write_report)

    code = cli.main(["run", str(run_folder_with_cal)])

    assert code == cli.EXIT_OK
    assert captured["converted_dir"] == run_folder_with_cal / "Converted"
    assert captured["run_name"] == run_folder_with_cal.name
    # cli.py never resolves DEFAULT_OUT_DIR to an absolute path. It stays relative to
    # whatever the process cwd is (here, chdir'd to tmp_path), matching how a real user's
    # shell would interpret a bare "output/" the same way.
    expected_relative = Path(cli.DEFAULT_OUT_DIR) / f"{run_folder_with_cal.name}_report.xlsx"
    assert captured["report_out_path"] == expected_relative
    assert (tmp_path / expected_relative).exists()

    out = capsys.readouterr().out
    assert "[calibrate]" in out  # progress callback was actually invoked
    assert "samples: 1" in out
    assert "blanks: 1" in out
    assert "flags: 1" in out
    assert str(expected_relative) in out


def test_run_honors_out_flag(run_folder_with_cal: Path, tmp_path, monkeypatch):
    monkeypatch.setattr(
        "quant209.core.pipeline.run_analysis",
        lambda entries, converted_dir, progress=None, qaqc=None, **kw: _stub_run_result(),
    )

    def fake_write_report(result, out_path, run_name=None):
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_bytes(b"x")
        return out_path

    monkeypatch.setattr("quant209.report.excel.write_report", fake_write_report)

    custom_out = tmp_path / "myout"
    code = cli.main(["run", str(run_folder_with_cal), "--out", str(custom_out)])

    assert code == cli.EXIT_OK
    assert (custom_out / f"{run_folder_with_cal.name}_report.xlsx").exists()


def test_run_passes_k_into_qaqc_params(run_folder_with_cal: Path, tmp_path, monkeypatch):
    from quant209.core.qaqc import QaqcParams

    captured = {}

    def fake_run_analysis(entries, converted_dir, progress=None, qaqc=None, **kwargs):
        captured["qaqc"] = qaqc
        return _stub_run_result()

    def fake_write_report(result, out_path, run_name=None):
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_bytes(b"x")
        return out_path

    monkeypatch.setattr("quant209.core.pipeline.run_analysis", fake_run_analysis)
    monkeypatch.setattr("quant209.report.excel.write_report", fake_write_report)

    code = cli.main(
        ["run", str(run_folder_with_cal), "--out", str(tmp_path / "out"), "--k", "4.5"]
    )

    assert code == cli.EXIT_OK
    assert isinstance(captured["qaqc"], QaqcParams)
    assert captured["qaqc"].k_sd == pytest.approx(4.5)


def test_run_default_k_when_not_given(run_folder_with_cal: Path, tmp_path, monkeypatch):
    from quant209.core.qaqc import QaqcParams

    captured = {}

    def fake_run_analysis(entries, converted_dir, progress=None, qaqc=None, **kwargs):
        captured["qaqc"] = qaqc
        return _stub_run_result()

    def fake_write_report(result, out_path, run_name=None):
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_bytes(b"x")
        return out_path

    monkeypatch.setattr("quant209.core.pipeline.run_analysis", fake_run_analysis)
    monkeypatch.setattr("quant209.report.excel.write_report", fake_write_report)

    cli.main(["run", str(run_folder_with_cal), "--out", str(tmp_path / "out")])

    assert captured["qaqc"].k_sd == QaqcParams().k_sd


# --------------------------------------------------------------------------------------
# `run --validate`
# --------------------------------------------------------------------------------------


def test_run_validate_writes_validation_report_when_truth_present(
    run_folder_with_cal: Path, tmp_path, monkeypatch, capsys
):
    monkeypatch.setattr(
        "quant209.core.pipeline.run_analysis",
        lambda entries, converted_dir, progress=None, qaqc=None, **kw: _stub_run_result(),
    )

    def fake_write_report(result, out_path, run_name=None):
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_bytes(b"x")
        return out_path

    monkeypatch.setattr("quant209.report.excel.write_report", fake_write_report)

    from quant209.report.truth import RunTruth

    non_empty_truth = RunTruth(
        areas=pd.DataFrame({"S1": {"1": 1000.0}}),
        masses=pd.DataFrame({"S1": {"1": 1.0}}),
        recoveries=pd.DataFrame({"S1": {"mono": 0.9}}),
    )
    monkeypatch.setattr("quant209.report.truth.read_run_truth", lambda folder: non_empty_truth)

    captured = {}

    def fake_generate_validation_report(result, truth, out_dir, run_name, **kwargs):
        captured["out_dir"] = out_dir
        captured["run_name"] = run_name
        return _stub_validation_report(Path(out_dir))

    monkeypatch.setattr(
        "quant209.report.validation.generate_validation_report", fake_generate_validation_report
    )

    out_dir = tmp_path / "out"
    code = cli.main(["run", str(run_folder_with_cal), "--out", str(out_dir), "--validate"])

    assert code == cli.EXIT_OK
    assert captured["out_dir"] == out_dir / "validation" / run_folder_with_cal.name
    assert captured["run_name"] == run_folder_with_cal.name

    out = capsys.readouterr().out
    assert "[validate]" in out
    assert "PASS" in out


def test_run_validate_skips_cleanly_when_no_truth_workbook(
    run_folder_with_cal: Path, tmp_path, monkeypatch, capsys
):
    monkeypatch.setattr(
        "quant209.core.pipeline.run_analysis",
        lambda entries, converted_dir, progress=None, qaqc=None, **kw: _stub_run_result(),
    )

    def fake_write_report(result, out_path, run_name=None):
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_bytes(b"x")
        return out_path

    monkeypatch.setattr("quant209.report.excel.write_report", fake_write_report)

    from quant209.report.truth import RunTruth

    empty_truth = RunTruth(areas=pd.DataFrame(), masses=pd.DataFrame(), recoveries=pd.DataFrame())
    monkeypatch.setattr("quant209.report.truth.read_run_truth", lambda folder: empty_truth)

    generate_called = []
    monkeypatch.setattr(
        "quant209.report.validation.generate_validation_report",
        lambda *a, **kw: generate_called.append(True),
    )

    code = cli.main(
        ["run", str(run_folder_with_cal), "--out", str(tmp_path / "out"), "--validate"]
    )

    assert code == cli.EXIT_OK  # --validate is opportunistic on `run`, not a hard requirement
    assert not generate_called
    assert "no analyst quantification workbook" in capsys.readouterr().out.lower()


# --------------------------------------------------------------------------------------
# `run --from-raw`
# --------------------------------------------------------------------------------------


def test_run_from_raw_requires_msconvert(run_folder_with_cal: Path, capsys):
    code = cli.main(["run", str(run_folder_with_cal), "--from-raw"])
    assert code == cli.EXIT_USAGE
    assert "--msconvert" in capsys.readouterr().err


def test_run_from_raw_converts_before_analyzing(tmp_path, monkeypatch, capsys):
    folder = tmp_path / "RAWRUN"
    folder.mkdir()
    (folder / "AN_Cal01_20240101.D").mkdir()
    exe = tmp_path / "msconvert.exe"
    exe.write_text("stub")

    plan = ConversionPlan(
        to_convert=[FileEntry(path=folder / "AN_Cal01_20240101.D", file_type="cal")]
    )
    monkeypatch.setattr(cli, "plan_conversion", lambda entries, out_dir: plan)

    convert_calls = []

    def fake_convert_many(sources, msconvert_exe, out_dir, extra_args=(), progress=None):
        convert_calls.append((list(sources), msconvert_exe, out_dir))
        if progress:
            progress(0, 1, "AN_Cal01_20240101.D")
        return [ConversionResult(source=sources[0], output=out_dir / "AN_Cal01_20240101.mzML",
                                  ok=True, message="ok")]

    monkeypatch.setattr(cli, "convert_many", fake_convert_many)

    monkeypatch.setattr(
        "quant209.core.pipeline.run_analysis",
        lambda entries, converted_dir, progress=None, qaqc=None, **kw: _stub_run_result(),
    )

    def fake_write_report(result, out_path, run_name=None):
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_bytes(b"x")
        return out_path

    monkeypatch.setattr("quant209.report.excel.write_report", fake_write_report)

    code = cli.main(
        ["run", str(folder), "--from-raw", "--msconvert", str(exe), "--out", str(tmp_path / "out")]
    )

    assert code == cli.EXIT_OK
    assert len(convert_calls) == 1
    assert convert_calls[0][1] == exe
    assert convert_calls[0][2] == folder / "Converted"
    assert "[convert]" in capsys.readouterr().out


def test_run_from_raw_fails_cleanly_on_conversion_error(tmp_path, monkeypatch, capsys):
    folder = tmp_path / "RAWRUN"
    folder.mkdir()
    (folder / "AN_Cal01_20240101.D").mkdir()
    exe = tmp_path / "msconvert.exe"
    exe.write_text("stub")

    plan = ConversionPlan(
        to_convert=[FileEntry(path=folder / "AN_Cal01_20240101.D", file_type="cal")]
    )
    monkeypatch.setattr(cli, "plan_conversion", lambda entries, out_dir: plan)
    monkeypatch.setattr(
        cli, "convert_many",
        lambda sources, msconvert_exe, out_dir, extra_args=(), progress=None: [
            ConversionResult(source=sources[0], output=None, ok=False, message="boom")
        ],
    )

    code = cli.main(["run", str(folder), "--from-raw", "--msconvert", str(exe)])

    assert code == cli.EXIT_FAILURE
    captured = capsys.readouterr()
    assert "boom" in captured.out or "boom" in captured.err


# --------------------------------------------------------------------------------------
# `convert`
# --------------------------------------------------------------------------------------


def test_convert_prints_plan_and_summary(tmp_path, monkeypatch, capsys):
    folder = tmp_path / "RUN"
    folder.mkdir()
    (folder / "AN_S1_20240101.D").mkdir()
    (folder / "AN_S2_20240101.D").mkdir()
    exe = tmp_path / "msconvert.exe"
    exe.write_text("stub")

    plan = ConversionPlan(
        to_convert=[FileEntry(path=folder / "AN_S1_20240101.D", file_type="sample")],
        skip=[FileEntry(path=folder / "AN_S2_20240101.D", file_type="sample")],
    )
    monkeypatch.setattr(cli, "plan_conversion", lambda entries, out_dir: plan)

    def fake_convert_many(sources, msconvert_exe, out_dir, extra_args=(), progress=None):
        if progress:
            progress(0, 1, Path(sources[0]).name)
        return [ConversionResult(source=Path(sources[0]), output=out_dir / "AN_S1_20240101.mzML",
                                  ok=True, message="ok")]

    monkeypatch.setattr(cli, "convert_many", fake_convert_many)

    code = cli.main(["convert", str(folder), "--msconvert", str(exe)])

    assert code == cli.EXIT_OK
    out = capsys.readouterr().out
    assert "plan" in out.lower()
    assert "1 to convert" in out
    assert "1 already converted" in out
    assert "[convert]" in out


def test_convert_defaults_out_dir_to_converted_subfolder(tmp_path, monkeypatch):
    folder = tmp_path / "RUN"
    folder.mkdir()
    exe = tmp_path / "msconvert.exe"
    exe.write_text("stub")

    captured = {}

    def fake_plan_conversion(entries, out_dir):
        captured["out_dir"] = out_dir
        return ConversionPlan()

    monkeypatch.setattr(cli, "plan_conversion", fake_plan_conversion)

    code = cli.main(["convert", str(folder), "--msconvert", str(exe)])

    assert code == cli.EXIT_OK
    assert captured["out_dir"] == folder / "Converted"


def test_convert_honors_out_flag(tmp_path, monkeypatch):
    folder = tmp_path / "RUN"
    folder.mkdir()
    exe = tmp_path / "msconvert.exe"
    exe.write_text("stub")
    custom_out = tmp_path / "somewhere_else"

    captured = {}

    def fake_plan_conversion(entries, out_dir):
        captured["out_dir"] = out_dir
        return ConversionPlan()

    monkeypatch.setattr(cli, "plan_conversion", fake_plan_conversion)

    cli.main(["convert", str(folder), "--msconvert", str(exe), "--out", str(custom_out)])

    assert captured["out_dir"] == custom_out


def test_convert_exits_1_on_folder_not_found(tmp_path, capsys):
    exe = tmp_path / "msconvert.exe"
    exe.write_text("stub")
    code = cli.main(["convert", str(tmp_path / "nope"), "--msconvert", str(exe)])
    assert code == cli.EXIT_FAILURE


def test_convert_exits_1_when_a_conversion_fails(tmp_path, monkeypatch, capsys):
    folder = tmp_path / "RUN"
    folder.mkdir()
    exe = tmp_path / "msconvert.exe"
    exe.write_text("stub")

    plan = ConversionPlan(to_convert=[FileEntry(path=folder / "AN_S1_20240101.D", file_type="sample")])
    monkeypatch.setattr(cli, "plan_conversion", lambda entries, out_dir: plan)
    monkeypatch.setattr(
        cli, "convert_many",
        lambda sources, msconvert_exe, out_dir, extra_args=(), progress=None: [
            ConversionResult(source=Path(sources[0]), output=None, ok=False, message="msconvert crashed")
        ],
    )

    code = cli.main(["convert", str(folder), "--msconvert", str(exe)])

    assert code == cli.EXIT_FAILURE
    captured = capsys.readouterr()
    assert "msconvert crashed" in captured.out or "msconvert crashed" in captured.err


# --------------------------------------------------------------------------------------
# `validate`
# --------------------------------------------------------------------------------------


def test_validate_writes_report_and_prints_headline(run_folder_with_cal: Path, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(
        "quant209.core.pipeline.run_analysis",
        lambda entries, converted_dir, progress=None, **kw: _stub_run_result(),
    )

    from quant209.report.truth import RunTruth

    non_empty_truth = RunTruth(
        areas=pd.DataFrame({"S1": {"1": 1000.0}}),
        masses=pd.DataFrame({"S1": {"1": 1.0}}),
        recoveries=pd.DataFrame({"S1": {"mono": 0.9}}),
    )
    monkeypatch.setattr("quant209.report.truth.read_run_truth", lambda folder: non_empty_truth)

    captured = {}

    def fake_generate_validation_report(result, truth, out_dir, run_name, **kwargs):
        captured["out_dir"] = out_dir
        return _stub_validation_report(Path(out_dir))

    monkeypatch.setattr(
        "quant209.report.validation.generate_validation_report", fake_generate_validation_report
    )

    out_dir = tmp_path / "out"
    code = cli.main(["validate", str(run_folder_with_cal), "--out", str(out_dir)])

    assert code == cli.EXIT_OK
    assert captured["out_dir"] == out_dir / "validation" / run_folder_with_cal.name
    out = capsys.readouterr().out
    assert "verdict: PASS" in out
    assert "recovery MAE" in out


def test_validate_errors_cleanly_with_no_truth_workbook(run_folder_with_cal: Path, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(
        "quant209.core.pipeline.run_analysis",
        lambda entries, converted_dir, progress=None, **kw: _stub_run_result(),
    )

    from quant209.report.truth import RunTruth

    empty_truth = RunTruth(areas=pd.DataFrame(), masses=pd.DataFrame(), recoveries=pd.DataFrame())
    monkeypatch.setattr("quant209.report.truth.read_run_truth", lambda folder: empty_truth)

    code = cli.main(["validate", str(run_folder_with_cal), "--out", str(tmp_path / "out")])

    assert code == cli.EXIT_FAILURE
    err = capsys.readouterr().err
    assert "no analyst quantification workbook" in err.lower()


def test_validate_on_nonexistent_folder_exits_1(tmp_path, capsys):
    code = cli.main(["validate", str(tmp_path / "nope")])
    assert code == cli.EXIT_FAILURE


# --------------------------------------------------------------------------------------
# --debug
# --------------------------------------------------------------------------------------


def test_debug_flag_reraises_instead_of_printing_one_line(run_folder_with_cal: Path, tmp_path, monkeypatch):
    def boom(entries, converted_dir, progress=None, qaqc=None, **kw):
        raise RuntimeError("kaboom")

    monkeypatch.setattr("quant209.core.pipeline.run_analysis", boom)

    with pytest.raises(RuntimeError, match="kaboom"):
        cli.main(["--debug", "run", str(run_folder_with_cal), "--out", str(tmp_path / "out")])


def test_without_debug_flag_prints_one_line_and_returns_1(run_folder_with_cal: Path, tmp_path, monkeypatch, capsys):
    def boom(entries, converted_dir, progress=None, qaqc=None, **kw):
        raise RuntimeError("kaboom")

    monkeypatch.setattr("quant209.core.pipeline.run_analysis", boom)

    code = cli.main(["run", str(run_folder_with_cal), "--out", str(tmp_path / "out")])

    assert code == cli.EXIT_FAILURE
    err = capsys.readouterr().err
    assert "kaboom" in err
    assert "Traceback" not in err
