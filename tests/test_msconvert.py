"""Tests for quant209.convert.msconvert: the MSConvert CLI wrapper and conversion planner."""

from __future__ import annotations

from pathlib import Path

import pytest

from quant209.convert.msconvert import (
    LAB_MSCONVERT_SETTINGS,
    ConversionPlan,
    convert_file,
    is_truncated,
    plan_conversion,
)
from quant209.core.sequence import FileEntry

# --------------------------------------------------------------------------------------
# Command construction
# --------------------------------------------------------------------------------------


def test_command_includes_lab_settings_and_output_dir(tmp_path: Path, monkeypatch):
    """The subprocess argv reproduces the lab's MSConvertGUI settings exactly: --mzML,
    --64, --zlib, the titleMaker filter, and -o pointing at the output directory."""
    exe = tmp_path / "msconvert.exe"
    exe.write_text("stub exe")
    source = tmp_path / "AN_Hex16_20231002.D"
    source.mkdir()
    out_dir = tmp_path / "out"

    captured: dict[str, list[str]] = {}

    class _FakeProc:
        returncode = 1  # never claim success. Only what was invoked matters here
        stdout = ""
        stderr = "fake failure"

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        return _FakeProc()

    import quant209.convert.msconvert as msconvert_mod

    monkeypatch.setattr(msconvert_mod.subprocess, "run", fake_run)

    convert_file(source, exe, out_dir)

    cmd = captured["cmd"]
    assert cmd[0] == str(exe)
    assert cmd[1] == str(source)
    for flag in LAB_MSCONVERT_SETTINGS:
        assert flag in cmd
    assert "-o" in cmd
    assert cmd[cmd.index("-o") + 1] == str(out_dir)


def test_lab_settings_contain_exact_flag_set():
    assert "--mzML" in LAB_MSCONVERT_SETTINGS
    assert "--64" in LAB_MSCONVERT_SETTINGS
    assert "--zlib" in LAB_MSCONVERT_SETTINGS
    assert "--filter" in LAB_MSCONVERT_SETTINGS
    idx = LAB_MSCONVERT_SETTINGS.index("--filter")
    assert "titleMaker" in LAB_MSCONVERT_SETTINGS[idx + 1]


# --------------------------------------------------------------------------------------
# convert_file with a missing exe: must not raise
# --------------------------------------------------------------------------------------


def test_convert_file_with_nonexistent_exe_fails_cleanly(tmp_path: Path):
    source = tmp_path / "AN_Hex16_20231002.D"
    source.mkdir()
    missing_exe = tmp_path / "does_not_exist" / "msconvert.exe"
    out_dir = tmp_path / "out"

    result = convert_file(source, missing_exe, out_dir)

    assert result.ok is False
    assert result.output is None
    assert "not found" in result.message.lower()


# --------------------------------------------------------------------------------------
# is_truncated
# --------------------------------------------------------------------------------------


def test_is_truncated_true_for_header_only_stub(tmp_path: Path):
    """An ~8 KB header-only mzML (interrupted MSConvert) is flagged truncated."""
    stub = tmp_path / "truncated.mzML"
    header = (
        '<?xml version="1.0" encoding="utf-8"?>\n'
        '<indexedmzML xmlns="http://psi.hupo.org/ms/mzml">\n'
        '  <mzML xmlns="http://psi.hupo.org/ms/mzml" id="stub" version="1.1.0">\n'
        '    <run id="stub_run">\n'
    )
    # Pad well past a real header but safely under the 50 KB floor.
    padding = "    <!-- pad -->\n" * 300
    stub.write_text(header + padding, encoding="utf-8")
    assert stub.stat().st_size < 50_000

    assert is_truncated(stub) is True


def test_is_truncated_false_for_a_real_converted_file(reference_run: Path):
    converted = sorted((reference_run / "Converted").glob("*.mzML"))
    assert converted, "expected at least one converted mzML in the fixture run"
    real_file = converted[0]
    assert real_file.stat().st_size >= 50_000

    assert is_truncated(real_file) is False


def test_is_truncated_true_for_missing_file(tmp_path: Path):
    assert is_truncated(tmp_path / "nope.mzML") is True


def test_is_truncated_true_for_unparseable_garbage(tmp_path: Path):
    """A large-enough file that isn't XML at all (e.g. a disk-full write of junk) is not
    just small. It's not parseable either, and must still be flagged."""
    garbage = tmp_path / "garbage.mzML"
    garbage.write_bytes(b"\x00\x01\x02 not xml at all " * 5000)
    assert garbage.stat().st_size >= 50_000

    assert is_truncated(garbage) is True


# --------------------------------------------------------------------------------------
# plan_conversion
# --------------------------------------------------------------------------------------


def _entry(name: str, suffix: str = ".D", file_type: str = "sample") -> FileEntry:
    return FileEntry(path=Path(f"{name}{suffix}"), file_type=file_type)


def test_plan_conversion_classifies_missing_truncated_and_fresh(tmp_path: Path):
    out_dir = tmp_path / "Converted"
    out_dir.mkdir()

    # Fresh, existing, large, parseable output -> skip.
    fresh_name = "AN_SMP16_20231002"
    fresh_header = (
        '<?xml version="1.0" encoding="utf-8"?>\n'
        '<indexedmzML xmlns="http://psi.hupo.org/ms/mzml">\n'
        '  <mzML xmlns="http://psi.hupo.org/ms/mzml" id="stub" version="1.1.0">\n'
        '    <run id="stub_run">\n'
    )
    fresh_padding = "      <!-- pad -->\n" * 4000  # pushes size well past the 50 KB floor
    (out_dir / f"{fresh_name}.mzML").write_text(fresh_header + fresh_padding, encoding="utf-8")

    # Existing but truncated (header-only, tiny) output -> reconvert.
    truncated_name = "AN_Ref04_20231002"
    (out_dir / f"{truncated_name}.mzML").write_bytes(b"<mzML></mzML>" + b" " * 100)

    # No output at all -> convert.
    missing_name = "AN_SMP17_20231002"

    entries = [_entry(fresh_name), _entry(truncated_name), _entry(missing_name)]
    plan = plan_conversion(entries, out_dir)

    assert isinstance(plan, ConversionPlan)
    assert [e.path.stem for e in plan.skip] == [fresh_name]
    assert [e.path.stem for e in plan.reconvert] == [truncated_name]
    assert [e.path.stem for e in plan.to_convert] == [missing_name]
    assert [e.path.stem for e in plan.pending] == [missing_name, truncated_name]


def test_plan_conversion_includes_hex_files_when_explicitly_passed(tmp_path: Path):
    """plan_conversion does not special-case file type. Hex entries the caller (GUI/CLI)
    explicitly selected are planned like any other .D entry. Only the pipeline later
    declines to quantify them."""
    out_dir = tmp_path / "Converted"
    entries = [_entry("AN_Hex01_20231002", file_type="hex")]

    plan = plan_conversion(entries, out_dir)

    assert [e.path.stem for e in plan.to_convert] == ["AN_Hex01_20231002"]


def test_plan_conversion_treats_mzml_entries_as_already_done(tmp_path: Path):
    out_dir = tmp_path / "Converted"
    entry = _entry("AN_SMP16_20231002", suffix=".mzML")

    plan = plan_conversion([entry], out_dir)

    assert plan.skip == [entry]
    assert plan.to_convert == []
    assert plan.reconvert == []
