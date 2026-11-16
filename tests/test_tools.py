"""Tests for the command-line tools under ``tools/``.

The tools are scripts, not a package, so they are loaded by path. Everything here runs on
data made up inside the test: the benchmark outputs the figure tool normally reads are
laboratory data and never enter the repository.
"""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import pytest

from quant209.core.config import load_canonical_peaks
from quant209.core.peaks import IntegrationParams

TOOLS = Path(__file__).resolve().parent.parent / "tools"


def _load(name: str):
    """Import ``tools/<name>.py`` as a module."""
    spec = importlib.util.spec_from_file_location(f"tools_{name}", TOOLS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------------------
# Parameter overrides on the two bench tools
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("tool", ["bench_areas", "bench_quant"])
def test_a_label_list_can_be_overridden_or_emptied_from_the_command_line(tool):
    """``--param forced_riders=...`` has to work, or the rules cannot be measured.

    The evidence for each named rule is a benchmark run with it and one without it. The
    override parser used to coerce everything that was not a bool or an int to a float.
    """
    apply = _load(tool)._apply_overrides
    params = apply(IntegrationParams(), ["forced_riders=83, 117,115", "forced_present=", "min_sn=2.5"])
    assert params.forced_riders == ("83", "117", "115")
    assert params.forced_present == (), "an empty value switches the rule off"
    assert params.min_sn == 2.5
    assert isinstance(params.forced_riders, tuple)


# --------------------------------------------------------------------------------------
# The validation figures, drawn from a made-up comparison
# --------------------------------------------------------------------------------------

AREA_FIELDS = ["run", "sample", "file_type", "label", "homolog", "role",
               "manual", "auto", "flags", "method", "sn", "rt_apex"]
MASS_FIELDS = ["run", "sample", "file_type", "label", "homolog", "role",
               "manual", "auto", "quantified", "flags", "method", "rt_apex", "detected"]
AROCLORS = ["1221", "1232", "1016", "1242", "1248", "1254", "1260", "1262", "1268"]
HOMOLOGS = ["mono", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona", "deca"]


def _write(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fields)
        writer.writeheader()
        writer.writerows(rows)


def _aroclor_profiles(natives) -> dict[str, dict[str, float]]:
    """Nine made-up mixtures, each a bell over the homolog series, so they are tellable apart."""
    profiles = {}
    for k, name in enumerate(AROCLORS):
        centre = 0.8 + k * 0.9
        weights = {
            a.label: float(np.exp(-0.5 * ((HOMOLOGS.index(a.homolog) - centre) / 1.1) ** 2))
            * (1.0 + 0.5 * np.sin(1.7 * i + k))
            for i, a in enumerate(natives)
        }
        total = sum(weights.values())
        profiles[name] = {label: w / total for label, w in weights.items()}
    return profiles


@pytest.fixture(scope="module")
def comparison(tmp_path_factory) -> tuple[Path, Path]:
    """A pooled comparison shaped like ``tools/master_bench.py`` output, and an Aroclor table.

    Three runs of two blanks and eight samples. Each sample is one of the mixtures at its own
    concentration. The automated value follows the manual one closely for most rows, and the
    failure modes the figures exist to show are each planted somewhere: an analyte the presence
    test missed, a rider, a forced peak, a sample the pipeline refused to quantify.
    """
    rng = np.random.default_rng(209)
    folder = tmp_path_factory.mktemp("comparison")
    analytes = load_canonical_peaks()
    natives = [a for a in analytes if a.role == "native"]
    profiles = _aroclor_profiles(natives)

    areas, masses, recoveries, runs = [], [], [], []
    for r in range(3):
        run = f"RUN_{'ABC'[r]}_2030010{r + 1}"
        files = [(f"Blank{b + 1}", "blank") for b in range(2)] + [(f"S{s + 1:02d}", "sample") for s in range(8)]
        for s, (sample, file_type) in enumerate(files):
            mixture = profiles[AROCLORS[(r * 3 + s) % len(AROCLORS)]]
            total_ng = 0.02 if file_type == "blank" else float(10 ** rng.uniform(1.0, 3.5))
            refused = file_type == "sample" and r == 2 and s == 9
            for i, a in enumerate(analytes):
                rt = float(a.cal_rt)
                if a.role != "native":
                    manual_area = float(rng.uniform(2e5, 4e5))
                    areas.append(dict(run=run, sample=sample, file_type=file_type, label=a.label,
                                      homolog=a.homolog, role=a.role, manual=manual_area,
                                      auto=manual_area * float(rng.normal(1.0, 0.01)), flags="-",
                                      method="valley_drop", sn=5000.0, rt_apex=rt))
                    continue
                manual_mass = total_ng * mixture[a.label] * float(rng.lognormal(0.0, 0.3))
                manual_area = manual_mass * 4000.0
                flags, method, ratio = "-", "valley_drop" if i % 3 else "rt_window", float(rng.lognormal(0.0, 0.04))
                draw = rng.random()
                if draw < 0.06:
                    flags, method, ratio = "absent", "rt_window", 0.003
                elif draw < 0.10:
                    flags, method, ratio = "shoulder", "deconvolved" if i % 2 else "tangent_skim", float(rng.lognormal(0.0, 0.3))
                elif draw < 0.12:
                    flags, method, ratio = "forced", "deconvolved", float(rng.lognormal(0.0, 0.25))
                elif draw < 0.16:
                    flags = "clustered"
                elif draw < 0.19:
                    flags, ratio = "low_sn", float(rng.lognormal(0.0, 0.5))
                auto_area = manual_area * ratio
                if manual_area < 30 and rng.random() < 0.5:
                    manual_area = manual_mass = 0.0
                areas.append(dict(run=run, sample=sample, file_type=file_type, label=a.label,
                                  homolog=a.homolog, role=a.role, manual=manual_area, auto=auto_area,
                                  flags=flags, method=method, sn=max(auto_area / 40.0, 0.5), rt_apex=rt))
                masses.append(dict(run=run, sample=sample, file_type=file_type, label=a.label,
                                   homolog=a.homolog, role=a.role, manual=manual_mass,
                                   auto=0.0 if refused else auto_area / 4000.0,
                                   quantified=not refused, flags=flags, method=method, rt_apex=rt,
                                   detected="absent" not in flags))
            for homolog in HOMOLOGS:
                value = float(rng.uniform(0.6, 1.0))
                recoveries.append(dict(run=run, sample=sample, homolog=homolog, manual=value,
                                       auto=value + float(rng.normal(0.0, 0.03))))
        runs.append(dict(run=run, status="ok", files=len(files) + 3, samples=8, truth_cols=len(files),
                         joined_areas=len(files) * len(analytes), joined_masses=len(files) * len(natives),
                         fit_ok=True, fit_reason="", seconds=10.0, note=""))
    runs.append(dict(run="RUN_D_20300104", status="no-truth", files="", samples="", truth_cols=0,
                     joined_areas="", joined_masses="", fit_ok="", fit_reason="", seconds=0.0,
                     note="no analyst quantification in this folder"))

    _write(folder / "areas.csv", AREA_FIELDS, areas)
    _write(folder / "masses.csv", MASS_FIELDS, masses)
    _write(folder / "recoveries.csv", ["run", "sample", "homolog", "manual", "auto"], recoveries)
    _write(folder / "runs.csv", list(runs[0]), runs)

    reference = folder / "aroclors.csv"
    with reference.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["pcb", *AROCLORS])
        for a in natives:
            writer.writerow([a.label, *(profiles[name][a.label] for name in AROCLORS)])
    return folder, reference


def _figures(out: Path) -> dict[str, str]:
    return {p.name: p.read_text(encoding="utf-8") for p in sorted(out.glob("*.svg"))}


def test_every_validation_figure_is_drawn_from_a_comparison(comparison, tmp_path, capsys):
    """The whole set draws from nothing but a comparison folder, and every file is an SVG."""
    folder, reference = comparison
    tool = _load("validation_figures")
    out = tmp_path / "figures"
    assert tool.main(["--data", str(folder), "--out", str(out), "--aroclor", str(reference),
                      "--root", str(tmp_path)]) == 0
    figures = _figures(out)
    expected = sum(len(module.FIGURES) for module in tool.MODULES)
    assert len(figures) == expected and expected >= 49
    for name, svg in figures.items():
        root = ElementTree.fromstring(svg)
        assert root.tag.endswith("svg"), name
        assert "nan" not in svg.lower().replace("penta", "").replace("nona", ""), f"{name}: a NaN was drawn"
    assert f"wrote {expected} figures" in capsys.readouterr().out


def test_the_fingerprint_figures_are_skipped_without_a_reference(comparison, tmp_path, capsys):
    """No Aroclor table is not an error. The figures that need one are left out, and said so."""
    folder, _reference = comparison
    tool = _load("validation_figures")
    out = tmp_path / "figures"
    assert tool.main(["--data", str(folder), "--out", str(out), "--aroclor", str(tmp_path / "none.csv"),
                      "--root", str(tmp_path)]) == 0
    needing = sum(1 for m in tool.MODULES for f in m.FIGURES if getattr(f, "needs_aroclors", False))
    total = sum(len(m.FIGURES) for m in tool.MODULES)
    assert needing > 0 and len(_figures(out)) == total - needing
    assert "fingerprint figures are skipped" in capsys.readouterr().out


def test_one_figure_can_be_drawn_alone_and_in_either_style(comparison, tmp_path):
    """``--only`` picks figures by number, and ``--as-first-drawn`` leaves out the later additions."""
    folder, reference = comparison
    tool = _load("validation_figures")
    common = ["--data", str(folder), "--aroclor", str(reference), "--root", str(tmp_path), "--only", "35"]
    assert tool.main([*common, "--out", str(tmp_path / "now")]) == 0
    assert tool.main([*common, "--out", str(tmp_path / "then"), "--as-first-drawn"]) == 0
    now, then = _figures(tmp_path / "now"), _figures(tmp_path / "then")
    assert len(now) == len(then) == 1 and list(now) == list(then)
    (name,) = now
    assert len(now[name]) > len(then[name]), "the legend and notes are the later additions"
