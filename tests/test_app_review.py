"""Behavioural tests for the Review dashboard and the Export page.

Same discipline as ``tests/test_app_pages.py``: ``QT_QPA_PLATFORM=offscreen`` before
pytest-qt builds the session ``QApplication``, a throwaway INI-backed ``QSettings`` in every
fixture, and **real widgets throughout**. The assertions read what a model actually returns,
what a plot actually contains and which button is actually enabled.

What is synthetic, and what is not
------------------------------------
The ``RunResult`` container is a ``SimpleNamespace`` standing in for the real one, which
wants an alignment, anchors and two parameter objects that none of this code reads.
Everything *inside* it is the real thing: real :class:`~quant209.core.peaks.PeakResult`
objects with real index bounds, real :class:`~quant209.core.config.Analyte` rows, real
``pandas`` frames, real :class:`~quant209.core.mzml.Chromatogram` traces, and a
:class:`~quant209.core.quantify.CalibrationModel` built by the pipeline's own
``build_calibration``. A fixture that faked the science would let this page render a number
no run could produce and still pass.

The run in the fixture, and why each piece is there
----------------------------------------------------
Five files (two samples, a blank, a reference and a calibration injection) over two
transitions and five analytes:

* three tetra natives on ``292.0 -> 222.0``: ``10`` ordinary, ``11`` **clustered** (so the
  overlay's second baseline hue is exercised), ``12`` **absent** (so the peak table's blank
  S/N rule and the calibration guard's refusal are both reachable).
* the tetra surrogate ``13C 52`` and the internal standard ``30D`` on ``304.0 -> 234.0``: the
  internal standard is what makes ``build_calibration`` produce real response factors rather
  than refusing every analyte for want of one.

Sample ``S2`` carries the **out-of-range tetra recovery** (42%) and one ``low_sn``
integration, which is what makes the default run amber. The knobs on :func:`build_result`
turn each of those off (green) or fail the run-fit guard (red).
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from PySide6.QtCore import QSettings, Qt

from quant209.app import charts as charts_module
from quant209.app.main import configure_app
from quant209.app.state import STEP_EXPORT, STEP_REVIEW, AppState
from quant209.app.workflow import review_page as rp
from quant209.app.workflow.export_page import (
    DEFAULT_OUT_DIRNAME,
    REPORT_SUFFIX,
    VALIDATION_DIRNAME,
    ExportPage,
    truth_workbooks,
)
from quant209.app.workflow.review_page import (
    FLAG_KINDS,
    OUT_OF_BAND_ROLE,
    TAB_CALIBRATION,
    TAB_OVERVIEW,
    TAB_QAQC,
    TAB_RECOVERIES,
    TAB_SAMPLES,
    VALUE_ROLE,
    ReviewPage,
)
from quant209.core.config import Analyte
from quant209.core.mzml import Chromatogram
from quant209.core.peaks import PeakResult
from quant209.core.quantify import build_calibration

# --------------------------------------------------------------------------------------
# The synthetic run
# --------------------------------------------------------------------------------------

T1 = (292.0, 222.0)  # three tetra natives
T2 = (304.0, 234.0)  # the tetra surrogate and the internal standard

ANALYTES = [
    Analyte("10", 292.0, 222.0, "tetra", 1, "native", "30D", 5.0, 10.00),
    Analyte("11", 292.0, 222.0, "tetra", 2, "native", "30D", 5.0, 10.30),
    Analyte("12", 292.0, 222.0, "tetra", 3, "native", "30D", 5.0, 10.60),
    Analyte("13C 52", 304.0, 234.0, "", 1, "surrogate_13C", "30D", 9.84, 10.10),
    Analyte("30D", 304.0, 234.0, "", 2, "IS_d30", "30D", 10.4125, 10.20),
]

T0, DT, NPOINTS = 9.5, 0.005, 301
TIMES = np.round(T0 + DT * np.arange(NPOINTS), 6)
BASELINE = 1000.0

FILE_IDS = ["C1", "R1", "B1", "S1", "S2"]
FILE_TYPES = {"C1": "cal", "R1": "ref", "B1": "blank", "S1": "sample", "S2": "sample"}


def _index(rt: float) -> int:
    return int(round((rt - T0) / DT))


def _trace(peaks: list[tuple[float, float]]) -> np.ndarray:
    y = np.full(NPOINTS, BASELINE)
    for centre, height in peaks:
        y = y + height * np.exp(-0.5 * ((TIMES - centre) / 0.02) ** 2)
    return y


CHROMS: dict[tuple[float, float], Chromatogram] = {
    T1: Chromatogram(*T1, TIMES, _trace([(10.00, 5.0e5), (10.30, 3.0e5), (10.60, 1.5e3)])),
    T2: Chromatogram(*T2, TIMES, _trace([(10.10, 4.0e5), (10.20, 6.0e5)])),
}


def _peak(label: str, rt: float, q1: float, q3: float, area: float, flags=()) -> PeakResult:
    apex = _index(rt)
    return PeakResult(
        label=label,
        q1=q1,
        q3=q3,
        rt_apex=rt,
        rt_start=rt - 0.06,
        rt_end=rt + 0.06,
        area=area,
        height=area / 40.0,
        baseline_start=BASELINE,
        baseline_end=BASELINE,
        # An absent peak's S/N is large and meaningless by construction (see
        # ``PeakResult.sn``). The peak table must refuse to print it anyway.
        sn=4.9e4 if "absent" in flags else 45.0,
        method="valley_drop",
        flags=list(flags),
        i_start=apex - 12,
        i_apex=apex,
        i_end=apex + 12,
        src_q1=q1,
        src_q3=q3,
    )


def _peaks_for(file_id: str, *, flagged: bool, cal_ok: bool = False) -> list[PeakResult]:
    absent = ("absent",) if file_id != "R1" else ()
    low_sn = ("low_sn",) if (flagged and file_id == "S2") else ()
    # ``12`` is absent everywhere, which is what makes the calibration guard refuse it a
    # response factor, and two tests depend on it. ``cal_ok`` places it in the
    # calibration injection instead, which is the only way to build a run this dashboard finds
    # *nothing* wrong with now that "congeners with no response factor" is one of the
    # Overview's attention classes.
    if cal_ok and file_id == "C1":
        absent, area_12 = (), 2.0e5
    else:
        area_12 = 8.0e2
    return [
        _peak("10", 10.00, *T1, 4.0e5, low_sn),
        _peak("11", 10.30, *T1, 2.2e5, ("clustered",)),
        _peak("12", 10.60, *T1, area_12, absent),
        _peak("13C 52", 10.10, *T2, 3.1e5),
        _peak("30D", 10.20, *T2, 5.4e5),
    ]


def build_result(
    tmp_path: Path, *, fit_ok=True, bad_recovery=True, flagged=True, blanks=True, cal_ok=False
):
    """A stand-in ``RunResult`` carrying real peaks, frames and a real calibration model."""
    converted = tmp_path / "run" / "Converted"
    converted.mkdir(parents=True, exist_ok=True)

    file_ids = [f for f in FILE_IDS if blanks or f != "B1"]
    files = {}
    for file_id in file_ids:
        mzml = converted / f"AN_{file_id}_20231002.mzML"
        mzml.write_text("<mzML/>", encoding="utf-8")
        files[file_id] = SimpleNamespace(
            entry=SimpleNamespace(sample_id=file_id, file_type=FILE_TYPES[file_id]),
            file_type=FILE_TYPES[file_id],
            mzml=mzml,
            peaks=_peaks_for(file_id, flagged=flagged, cal_ok=cal_ok),
            skipped_analytes=[],
        )

    sample_ids = ["S1", "S2"]
    blank_ids = ["B1"] if blanks else []
    reported = sample_ids + blank_ids

    recoveries = pd.DataFrame(
        index=list(rp.HOMOLOG_SURROGATE), columns=reported, dtype=float
    )
    tetra = [1.05, 0.42 if bad_recovery else 1.02] + ([0.98] if blanks else [])
    recoveries.loc["tetra", reported] = tetra

    labels = [a.label for a in ANALYTES]
    masses = pd.DataFrame(3.0, index=labels, columns=file_ids)
    corrected = masses.copy()
    censored = corrected.copy()
    censored.loc["10", "S2"] = 0.0  # the one cell the LOQ censoring zeroed
    loq = pd.Series(0.5, index=labels)

    calibration = build_calibration(files["C1"].peaks, ANALYTES, name="C1")

    return SimpleNamespace(
        analytes=ANALYTES,
        files=files,
        sample_ids=sample_ids,
        blank_ids=blank_ids,
        cal_ids=["C1"],
        ref_ids=["R1"],
        reported_ids=reported,
        calibration=calibration,
        recoveries=recoveries,
        masses=masses,
        masses_corrected=corrected,
        masses_censored=censored,
        loq=loq,
        flags={fid: [] for fid in file_ids},
        status="OK" if fit_ok else "FAILED-QC",
        status_reason="" if fit_ok else "only 12% of the labeled standards were found",
        fit=SimpleNamespace(
            fit_ok=fit_ok,
            reason="" if fit_ok else "only 12% of the labeled standards were found",
            summary="native 3/3 (100%), labeled 2/2 (100%)",
        ),
        notes=[],
        reference_for={fid: "R1" for fid in reported},
    )


# --------------------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------------------


@pytest.fixture
def settings(tmp_path: Path) -> QSettings:
    return QSettings(str(tmp_path / "quant209.ini"), QSettings.Format.IniFormat)


@pytest.fixture
def state(settings: QSettings, tmp_path: Path) -> AppState:
    app_state = AppState(settings=settings)
    app_state.set_run_folder(tmp_path / "run")
    return app_state


@pytest.fixture
def styled(qapp):
    configure_app(qapp)
    return qapp


@pytest.fixture
def chroms(monkeypatch):
    """Serve every file the same two synthetic transitions, and count the reads."""
    reads: list[Path] = []

    def read_srm_mzml(path):
        reads.append(Path(path))
        return dict(CHROMS)

    monkeypatch.setattr("quant209.core.mzml.read_srm_mzml", read_srm_mzml)
    return reads


@pytest.fixture
def result(tmp_path: Path):
    return build_result(tmp_path)


def review(state: AppState, qtbot, result=None) -> ReviewPage:
    if result is not None:
        state.set_result(result)
    page = ReviewPage(state)
    qtbot.addWidget(page)
    page.on_activated()
    return page


# --------------------------------------------------------------------------------------
# Overview: the verdict
# --------------------------------------------------------------------------------------


def test_a_clean_run_reads_green_and_says_what_was_checked(tmp_path: Path):
    summary = rp.summarise(
        build_result(tmp_path, bad_recovery=False, flagged=False, cal_ok=True)
    )
    tone, headline, why = rp.verdict(summary)

    assert (tone, headline) == ("success", "No problems found")
    assert "retention-time table fits" in why
    assert "recoveries within 50%–150%" in why
    assert summary.failing == []


def test_an_out_of_band_recovery_reads_amber_and_names_the_sample(tmp_path: Path):
    summary = rp.summarise(build_result(tmp_path))
    tone, headline, why = rp.verdict(summary)

    assert (tone, headline) == ("warning", "Review required")
    assert "1 of 3 recoveries outside 50%–150%" in why
    rows = [(item.where, item.what, item.value) for item in summary.failing]
    assert ("S2", "tetra surrogate recovery", "42%") in rows


def test_a_failed_run_fit_reads_red_and_quotes_the_reason_verbatim(tmp_path: Path):
    reason = "only 12% of the labeled standards were found"
    summary = rp.summarise(build_result(tmp_path, fit_ok=False))
    tone, headline, why = rp.verdict(summary)

    assert tone == "error"
    assert headline == "Run failed the retention-time fit check"
    assert why == reason, "the guard's own words, not a paraphrase"
    assert summary.failing[0].what == "retention-time fit"


def test_a_run_with_no_blanks_is_amber_because_nothing_was_censored(tmp_path: Path):
    summary = rp.summarise(build_result(tmp_path, bad_recovery=False, flagged=False, blanks=False))
    tone, _headline, why = rp.verdict(summary)

    assert tone == "warning"
    assert "no blanks" in why
    assert ("this run", "blanks") in [(i.where, i.what) for i in summary.failing]


def test_flagged_integrations_past_the_limit_are_amber(tmp_path: Path):
    summary = rp.summarise(build_result(tmp_path, bad_recovery=False))

    assert summary.flagged_integrations == 1
    assert summary.flagged_fraction > rp.FLAGGED_FRACTION_LIMIT
    assert rp.verdict(summary)[0] == "warning"
    assert ("S2", "flagged integrations", "1 of 5 peaks") in [
        (i.where, i.what, i.value) for i in summary.failing
    ]


# --------------------------------------------------------------------------------------
# Overview: the attention classes
#
# A "What" column that always reads the same value carries no information. These pin each
# class the table distinguishes, the scoping decision behind the two that could otherwise
# flood it, and the severity bands that decide what an analyst reads first.
# --------------------------------------------------------------------------------------


def _add_second_calibration(result, *, disagree_on: str = "10") -> None:
    """Give the run a second calibration injection whose RRF for one analyte is halved."""
    peaks = []
    for peak in result.files["C1"].peaks:
        area = peak.area * 2.0 if peak.label == disagree_on else peak.area
        peaks.append(_peak(peak.label, peak.rt_apex, peak.q1, peak.q3, area, tuple(peak.flags)))
    result.files["C2"] = SimpleNamespace(
        entry=SimpleNamespace(sample_id="C2", file_type="cal"),
        file_type="cal",
        mzml=result.files["C1"].mzml,
        peaks=peaks,
        skipped_analytes=[],
    )
    result.cal_ids = ["C1", "C2"]


def test_the_attention_table_names_more_than_one_kind_of_trouble(tmp_path: Path):
    """One What value over every row is the failure mode this table exists to avoid."""
    result = build_result(tmp_path)
    result.recoveries.loc["tetra", "B1"] = float("nan")
    result.masses_censored.loc[["10", "11", "12"], "S1"] = 0.0
    _add_second_calibration(result)

    summary = rp.summarise(result)
    kinds = {item.what for item in summary.failing}

    assert len(kinds) > 1, "a What column with one value in it is not a column"
    assert {
        "response factors",
        "response-factor spread",
        "no congener above the limit",
        "tetra surrogate recovery",
        "tetra recovery gap",
        "flagged integrations",
    } <= kinds


def test_congeners_with_no_response_factor_are_one_row_carrying_the_count(tmp_path: Path):
    """One row, not one per analyte: 185 identical What cells is the defect, not the fix."""
    summary = rp.summarise(build_result(tmp_path, bad_recovery=False, flagged=False))

    rows = [i for i in summary.failing if i.what == "response factors"]
    assert len(rows) == 1
    assert rows[0].where == "this run"
    assert "1 of 5 congeners have none" in rows[0].value
    # ...and it is gone entirely once the guard accepts every analyte.
    clean = rp.summarise(build_result(tmp_path, bad_recovery=False, flagged=False, cal_ok=True))
    assert [i for i in clean.failing if i.what == "response factors"] == []


def test_response_factors_that_disagree_between_injections_are_reported(tmp_path: Path):
    result = build_result(tmp_path, bad_recovery=False, flagged=False)
    assert rp.rrf_spread(rp.per_injection_models(result)) == [], "one injection has no spread"

    _add_second_calibration(result, disagree_on="10")
    models = rp.per_injection_models(result)

    assert set(models) == {"C1", "C2"}
    assert "10" in rp.rrf_spread(models), "halving one injection's area is a 47% RSD"
    rows = [i for i in rp.summarise(result).failing if i.what == "response-factor spread"]
    assert len(rows) == 1
    assert "2 calibration injections" in rows[0].value


def test_a_sample_that_reported_nothing_is_flagged_and_a_clean_blank_is_not(tmp_path: Path):
    """A blank censoring to zero is the limit working, not a problem. See ``summarise``."""
    result = build_result(tmp_path, bad_recovery=False, flagged=False)
    result.masses_censored.loc[["10", "11", "12"], ["S2", "B1"]] = 0.0

    rows = [i for i in rp.summarise(result).failing if i.what == "no congener above the limit"]

    assert [i.where for i in rows] == ["S2"]
    assert rows[0].value == "0 of 3 congeners cleared the limit of quantification"


def test_a_homolog_measured_elsewhere_and_missing_here_is_a_gap(tmp_path: Path):
    result = build_result(tmp_path, bad_recovery=False, flagged=False)
    result.recoveries.loc["tetra", "B1"] = float("nan")

    gaps = rp.recovery_gaps(result.recoveries, ["S1", "S2", "B1"])

    assert gaps == {"B1": ["tetra"]}, "the nine homologs this run never measured are not gaps"
    rows = [i for i in rp.summarise(result).failing if i.what == "tetra recovery gap"]
    assert [i.where for i in rows] == ["B1"]
    assert "13C 52" in rows[0].value and "uncorrected" in rows[0].value


def test_a_run_with_no_usable_reference_says_so_once(tmp_path: Path):
    result = build_result(tmp_path, bad_recovery=False, flagged=False)
    result.reference_for = {}

    rows = [i for i in rp.summarise(result).failing if i.what == "reference injection"]

    assert len(rows) == 1, "one run-level row, not one per reported file"
    assert rows[0].where == "this run"
    assert rows[0].tone == "error"


def test_an_unpaired_file_is_named_when_the_rest_of_the_run_was_paired(tmp_path: Path):
    result = build_result(tmp_path, bad_recovery=False, flagged=False)
    del result.reference_for["S2"]

    rows = [i for i in rp.summarise(result).failing if i.what == "reference pairing"]

    assert [i.where for i in rows] == ["S2"]


def test_attention_rows_sort_worst_class_first(tmp_path: Path):
    """Between classes the order is the policy. Magnitude only breaks ties inside one."""
    result = build_result(tmp_path, fit_ok=False)
    result.recoveries.loc["tetra", "B1"] = float("nan")
    result.masses_censored.loc[["10", "11", "12"], "S1"] = 0.0
    _add_second_calibration(result)

    order = [item.what for item in rp.summarise(result).failing]

    assert order.index("retention-time fit") == 0
    assert order.index("response factors") < order.index("response-factor spread")
    assert order.index("response-factor spread") < order.index("no congener above the limit")
    assert order.index("no congener above the limit") < order.index("tetra surrogate recovery")
    assert order.index("tetra surrogate recovery") < order.index("tetra recovery gap")
    assert order.index("tetra recovery gap") < order.index("flagged integrations")


def test_a_severe_item_never_outranks_the_class_above_it(tmp_path: Path):
    """The within-band tie-break is capped, so bands cannot interleave."""
    result = build_result(tmp_path)
    result.recoveries.loc["tetra", "S2"] = 12.0  # a wildly out-of-band recovery
    result.masses_censored.loc[["10", "11", "12"], "S1"] = 0.0

    summary = rp.summarise(result)
    by_what = {i.what: i.weight for i in summary.failing}

    assert by_what["tetra surrogate recovery"] < rp.W_NOTHING_ABOVE_LOQ
    assert by_what["no congener above the limit"] == rp.W_NOTHING_ABOVE_LOQ
    assert max(i.weight for i in summary.failing if i.what == "flagged integrations") < (
        rp.W_RECOVERY_GAP
    )


def test_the_stat_cards_count_the_run(state, qtbot, result, styled):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_OVERVIEW)

    assert tab.cards["samples"].value_label.text() == "2"
    assert tab.cards["blanks"].value_label.text() == "1"
    assert tab.cards["recovery"].value_label.text() == "42–105%"
    assert tab.cards["flagged"].value_label.text() == "1"
    assert tab.cards["censored"].value_label.text() == "1 of 9"
    assert tab.banner.tone == "warning"
    assert tab.table.isVisibleTo(tab)


def test_the_failing_table_is_empty_state_when_nothing_is_wrong(state, qtbot, tmp_path):
    state.set_result(build_result(tmp_path, bad_recovery=False, flagged=False, cal_ok=True))
    page = review(state, qtbot)
    tab = page.ensure_tab(TAB_OVERVIEW)

    assert tab.model.rowCount() == 0
    assert tab.empty.isVisibleTo(tab)
    assert not tab.table.isVisibleTo(tab)
    assert tab.empty.headline_label.text() == "Nothing needs attention"


def test_the_header_chip_carries_the_verdict(state, qtbot, result):
    page = review(state, qtbot, result)
    assert page.page_status() == ("Review required", "warning")

    state.set_result(build_result(Path(state.run_folder).parent, fit_ok=False))
    page.refresh()
    assert page.page_status() == ("Fit failed", "error")


# --------------------------------------------------------------------------------------
# Samples: the chromatogram browser
# --------------------------------------------------------------------------------------


def test_the_rail_groups_the_files_and_selects_the_first_sample(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_SAMPLES)

    texts = [tab.file_list.item(i).text() for i in range(tab.file_list.count())]
    assert texts == ["SAMPLES", "S1", "S2", "BLANKS", "B1", "REFERENCE", "R1",
                     "CALIBRATION", "C1"]
    assert tab.current_file_id() == "S1"
    assert tab.transition_combo.count() == 2
    assert tab.transition_combo.itemText(0).startswith("tetra")


def test_selecting_a_transition_draws_every_integration_on_it(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_SAMPLES)
    view = tab.chart

    assert view.trace is not None, "one acquired trace"
    assert len(view.region_items) == 3, "one shaded region per analyte on this transition"
    assert len(view.baseline_items) == 3
    assert len(view.apex_scatter.points()) == 3
    assert len(view.labels) == 3
    assert tab.model.rowCount() == 3
    assert [tab.model.row_at(i).label for i in range(3)] == ["10", "11", "12"]


def test_switching_transitions_replots_the_other_channel(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_SAMPLES)

    tab.transition_combo.setCurrentIndex(1)

    assert len(tab.chart.region_items) == 2
    assert [tab.model.row_at(i).label for i in range(2)] == ["13C 52", "30D"]


def test_switching_files_reads_once_and_then_serves_from_the_cache(
    state, qtbot, result, chroms
):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_SAMPLES)
    tab.file_list.setCurrentRow(2)  # S2
    reads_after_s2 = len(chroms)

    tab.file_list.setCurrentRow(1)  # back to S1, already cached
    tab.transition_combo.setCurrentIndex(1)

    assert len(chroms) == reads_after_s2, "a cached file is never re-read"
    assert len(page.cache) == 2


def test_the_absent_row_shows_no_signal_to_noise(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_SAMPLES)
    model = tab.model

    absent_row = [i for i in range(model.rowCount()) if model.row_at(i).label == "12"][0]
    assert model.data(model.index(absent_row, 3)) == ""
    assert "not a measurement" in model.data(
        model.index(absent_row, 3), Qt.ItemDataRole.ToolTipRole
    )
    ordinary_row = [i for i in range(model.rowCount()) if model.row_at(i).label == "10"][0]
    assert model.data(model.index(ordinary_row, 3)) == "45"


def test_a_flagged_row_is_tinted_and_carries_its_flags_as_chips(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_SAMPLES)
    tab.file_list.setCurrentRow(2)  # S2, the file with the low_sn integration
    model = tab.model

    row = [i for i in range(model.rowCount()) if model.row_at(i).label == "10"][0]
    assert model.data(model.index(row, 5), rp.FLAGS_ROLE) == ["low_sn"]
    assert model.data(model.index(row, 0), Qt.ItemDataRole.BackgroundRole) is not None
    clustered = [i for i in range(model.rowCount()) if model.row_at(i).label == "11"][0]
    assert model.data(model.index(clustered, 5), rp.FLAGS_ROLE) == ["clustered"]


def test_clicking_a_peak_row_zooms_the_view_to_that_peak(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_SAMPLES)
    view = tab.chart
    before = view.view_box.viewRange()[0]
    peak = tab.model.row_at(1)  # analyte "11" at 10.30 min

    tab.table.clicked.emit(tab.model.index(1, 0))

    (x0, x1), _y = view.zoom_target
    assert x0 == pytest.approx(peak.rt_start - charts_module.ZOOM_PAD_MIN)
    assert x1 == pytest.approx(peak.rt_end + charts_module.ZOOM_PAD_MIN)
    qtbot.waitUntil(
        lambda: (view.view_box.viewRange()[0][1] - view.view_box.viewRange()[0][0])
        < (before[1] - before[0]),
        timeout=3000,
    )


def test_labels_thin_out_instead_of_overlapping(state, qtbot, result, chroms, styled):
    """Every label is drawn, but only the ones that clear their neighbour are shown."""
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_SAMPLES)
    view = tab.chart
    view.resize(400, 240)
    view.view_box.setRange(xRange=(9.5, 11.0), padding=0)

    assert view.visible_label_count() <= len(view.labels)
    view.view_box.setRange(xRange=(10.25, 10.35), padding=0)
    assert view.visible_label_count() >= 1


def _crowded_view(qtbot):
    """Two analytes 0.03 min apart, closer together than one label is wide, on one trace.

    The rejected one elutes *first*, which is the arrangement a left-to-right thinning pass
    gets wrong: on real calibration data a near-zero ``absent`` estimate two seconds ahead of
    the dominant congener took the only label the window had room for.
    """
    view = charts_module.ChromatogramView()
    qtbot.addWidget(view)
    view.resize(200, 240)
    rejected = _peak("12", 10.20, *T1, 8.0e2, ("absent",))
    dominant = _peak("10", 10.23, *T1, 4.0e5)
    view.show_trace(TIMES, CHROMS[T1].intensities, [rejected, dominant])
    # The whole acquisition in view: 0.03 min is then two per cent of the width, well inside
    # one label's own box, so exactly one of the two names can be drawn.
    view.view_box.setRange(xRange=(TIMES[0], TIMES[-1]), padding=0)
    return view, rejected, dominant


def test_a_rejected_analyte_does_not_take_its_dominant_neighbours_label(qtbot, styled):
    view, _rejected, _dominant = _crowded_view(qtbot)

    shown = {mark.text for mark in view.labels if mark.item.isVisible()}

    assert shown == {"10"}, "the detected analyte is named; the absent estimate gives way"


def test_the_analyte_the_view_was_taken_to_keeps_its_label(qtbot, styled):
    """Whatever the ranking says, the analyte the analyst clicked is the one named."""
    view, rejected, _dominant = _crowded_view(qtbot)
    assert "12" not in {m.text for m in view.labels if m.item.isVisible()}

    view.zoom_to_peak(rejected, animate=False)

    assert "12" in {mark.text for mark in view.labels if mark.item.isVisible()}
    view.reset_zoom()
    assert view.pinned is None


def test_every_region_is_bounded_by_its_perpendicular_drops(state, qtbot, result, chroms):
    """Two hairlines per region, so abutting cluster baselines are separable by eye."""
    page = review(state, qtbot, result)
    view = page.ensure_tab(TAB_SAMPLES).chart

    points = sum(len(item.getData()[0]) for item in view.boundary_items)

    assert view.boundary_items
    assert points == 4 * len(view.region_items), "two two-point segments per region"


# --------------------------------------------------------------------------------------
# Calibration
# --------------------------------------------------------------------------------------


def test_the_calibration_table_lists_every_analyte_with_its_injections(
    state, qtbot, result, chroms
):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_CALIBRATION)

    assert tab.model.injections == ("C1",)
    assert tab.model.rowCount() == len(ANALYTES)
    # Answer first, evidence after: the response factor every mass is multiplied by comes
    # before the injections it was averaged from, so it survives a narrow window.
    assert tab.model.HEADERS == (
        "Analyte", "Homolog", "Transition", "Combined RRF", "C1", "Flag",
    )
    assert "RRFs from 1 injection" in tab.summary_label.text()


def test_the_guards_refusal_appears_in_the_flag_column(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_CALIBRATION)
    model = tab.model

    refused = [i for i in range(model.rowCount()) if model.row_at(i).label == "12"][0]
    assert model.display(model.row_at(refused), 3) == "—", "no combined RRF"
    assert "not detected" in model.display(model.row_at(refused), 5)
    ordinary = [i for i in range(model.rowCount()) if model.row_at(i).label == "10"][0]
    assert model.display(model.row_at(ordinary), 5) == ""
    assert model.display(model.row_at(ordinary), 3) != "—"
    assert "1 analyte flagged" in tab.summary_label.text()


def test_the_calibration_search_box_filters_rows(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_CALIBRATION)

    tab.search.setText("13C")

    assert tab.model.rowCount() == 1
    assert tab.model.row_at(0).label == "13C 52"


def test_selecting_a_calibration_row_previews_its_peak(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_CALIBRATION)

    row = [i for i in range(tab.model.rowCount()) if tab.model.row_at(i).label == "11"][0]
    tab.table.clicked.emit(tab.model.index(row, 0))

    assert tab.chart.trace is not None
    assert len(tab.chart.region_items) == 3, "the whole channel, so neighbours are visible"
    assert "11 in C1" in tab.preview_caption.text()
    assert tab.chart.zoom_target is not None


# --------------------------------------------------------------------------------------
# Recoveries
# --------------------------------------------------------------------------------------


def test_the_heatmap_model_matches_the_result_frame(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_RECOVERIES)
    model = tab.model

    assert model.homologs == list(rp.HOMOLOG_SURROGATE)
    assert model.files == ["S1", "S2", "B1"]
    tetra = model.homologs.index("tetra")
    for column, file_id in enumerate(model.files):
        expected = result.recoveries.loc["tetra", file_id]
        assert model.data(model.index(tetra, column), VALUE_ROLE) == pytest.approx(expected)
    assert model.data(model.index(tetra, 1)) == "42%"
    assert model.data(model.index(0, 0), VALUE_ROLE) is None, "an unmeasured homolog"
    assert model.data(model.index(0, 0)) == "—"


def test_out_of_band_cells_are_flagged_in_the_model(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    model = page.ensure_tab(TAB_RECOVERIES).model
    tetra = model.homologs.index("tetra")

    assert model.data(model.index(tetra, 1), OUT_OF_BAND_ROLE) is True  # 42%
    assert model.data(model.index(tetra, 0), OUT_OF_BAND_ROLE) is False  # 105%
    assert model.data(model.index(0, 0), OUT_OF_BAND_ROLE) is False  # no measurement


def test_the_heatmap_colours_diverge_around_one_hundred_percent():
    low, mid, high = (
        charts_module.recovery_color(0.2),
        charts_module.recovery_color(1.0),
        charts_module.recovery_color(1.9),
    )
    assert low.red() > low.blue(), "below the band leans amber"
    assert high.green() > high.red(), "above the band leans teal"
    assert mid != low and mid != high
    assert charts_module.recovery_in_band(1.0)
    assert not charts_module.recovery_in_band(0.42)


def test_clicking_a_recovery_cell_previews_the_surrogate_peak(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_RECOVERIES)
    tetra = tab.model.homologs.index("tetra")

    tab.table.clicked.emit(tab.model.index(tetra, 1))

    assert tab.chart.trace is not None
    assert "13C 52 in S2" in tab.preview_caption.text()
    assert "42%" in tab.preview_caption.text()


def test_the_surrogate_preview_is_a_strip_until_a_cell_is_chosen(state, qtbot, result, chroms):
    """Before any click the preview must not be half the screen. See ``RecoveriesTab``.

    An expanding placeholder in a card added at stretch 1 would let every *first* view of
    this tab become a 45%-of-screen empty white card carrying one grey sentence.
    """
    from PySide6.QtWidgets import QSizePolicy

    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_RECOVERIES)
    layout = tab._outer

    assert not tab.chart.isVisibleTo(tab), "no chart until there is a trace"
    assert tab.placeholder.sizePolicy().verticalPolicy() == QSizePolicy.Policy.Fixed
    assert layout.stretch(tab._preview_row) == 0, "the card hugs its one line"
    assert layout.stretch(tab._filler_row) == 1, "the surplus goes below it"
    strip = tab.preview_card.sizeHint().height()
    assert strip < 200, f"a title and one line, not a panel (got {strip} px)"

    tab.table.clicked.emit(tab.model.index(tab.model.homologs.index("tetra"), 1))
    tab.preview_card.layout().activate()  # size hints are cached until the layout re-runs

    assert tab.chart.isVisibleTo(tab)
    assert layout.stretch(tab._preview_row) == 1, "and takes the height once it has a trace"
    assert layout.stretch(tab._filler_row) == 0
    assert tab.preview_card.sizeHint().height() > strip, "the chart's own minimum, at least"


def test_a_preview_with_nothing_to_draw_collapses_back_to_the_strip(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_RECOVERIES)
    tab.table.clicked.emit(tab.model.index(tab.model.homologs.index("tetra"), 1))
    assert tab._outer.stretch(tab._preview_row) == 1

    # ``mono`` is tracked by a surrogate no file in this run carries a peak for.
    tab.table.clicked.emit(tab.model.index(tab.model.homologs.index("mono"), 0))

    assert not tab.chart.isVisibleTo(tab)
    assert tab.placeholder.isVisibleTo(tab)
    assert "13C 3" in tab.placeholder.text()
    assert tab._outer.stretch(tab._preview_row) == 0, "a message is a strip, not a panel"


def test_the_calibration_preview_still_fills_its_column(state, qtbot, result, chroms):
    """The other preview is a side panel beside a full-height table and must not collapse."""
    from PySide6.QtWidgets import QSizePolicy

    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_CALIBRATION)

    assert tab.placeholder.sizePolicy().verticalPolicy() == QSizePolicy.Policy.Expanding


def test_the_recovery_grid_asks_for_its_rows_and_no_more(state, qtbot, result, chroms):
    """The grid is a banner, not a pane: ten rows tall whatever the height available.

    Without a real size hint, an expanding grid leaves half the tab blank and strands the
    legend at the bottom of the emptiness, and resizing it from a ``resizeEvent`` instead
    paints the legend over its last two rows.
    """
    from PySide6.QtWidgets import QSizePolicy

    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_RECOVERIES)
    rows = tab.model.rowCount()

    assert tab.table.sizePolicy().verticalPolicy() == QSizePolicy.Policy.Fixed
    hint = tab.table.sizeHint().height()
    assert hint == tab.table.minimumSizeHint().height()
    # Its rows, a header and a reserved scrollbar, not a screenful.
    assert rows * rp.ROW_HEIGHT <= hint <= rows * rp.ROW_HEIGHT + 100


def test_a_run_with_no_recoveries_shows_an_empty_state(state, qtbot, tmp_path, chroms):
    result = build_result(tmp_path)
    result.recoveries = pd.DataFrame()
    state.set_result(result)
    page = review(state, qtbot)
    tab = page.ensure_tab(TAB_RECOVERIES)

    assert tab.empty.isVisibleTo(tab)
    assert not tab.table.isVisibleTo(tab)


# --------------------------------------------------------------------------------------
# QA/QC
# --------------------------------------------------------------------------------------


def test_every_flag_in_the_run_becomes_one_event(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_QAQC)
    kinds = [event.kind for event in tab._events]

    assert kinds.count("recovery") == 1
    assert kinds.count("clustered") == 5, "every file's analyte 11"
    assert kinds.count("absent") == 4, "every file but the reference"
    assert kinds.count("low_sn") == 1
    assert kinds.count("censored") == 1
    assert set(kinds) <= set(FLAG_KINDS)


def test_a_forced_peak_says_in_words_that_a_rule_produced_it():
    """``forced`` is the one flag that means "no test found this peak in this file".

    The QA/QC table is where a reviewer meets it, so its context has to say that in a
    sentence and not leave the cell empty: the number beside it was put there by a rule
    that names the congener, and that is what they need to know before trusting it.
    """
    peak = SimpleNamespace(rt_apex=42.41, sn=36.2, i_start=10, i_end=40)
    context = rp._flag_context("forced", peak)
    assert "rule" in context and "detection test" in context
    assert rp.FLAG_TONES["forced"] == "info", "how it was measured, not a defect in it"
    assert "forced" not in rp.QUALITY_FLAGS, "a forced peak is not counted as a flagged integration"


def test_the_filter_chips_reduce_the_rows(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_QAQC)
    everything = tab.model.rowCount()

    tab.kind_buttons["absent"].setChecked(True)
    tab._on_kind_clicked("absent")

    assert tab.model.rowCount() == 4
    assert tab.model.rowCount() < everything
    assert all(row.kind == "absent" for row in tab.model.rows)
    assert not tab.kind_buttons["All"].isChecked()
    assert f"4 of {everything} flag events" == tab.count_label.text()

    tab.kind_buttons["All"].setChecked(True)
    tab._on_kind_clicked("All")
    assert tab.model.rowCount() == everything


def test_the_sample_filter_narrows_to_one_file(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_QAQC)

    tab.file_combo.setCurrentIndex(tab.file_combo.findData("S2"))

    assert tab.model.rowCount() > 0
    assert {row.file_id for row in tab.model.rows} == {"S2"}


def test_turning_off_the_last_chip_falls_back_to_all(state, qtbot, result, chroms):
    """An empty filter row over an empty table is a dead end. All is the floor."""
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_QAQC)

    tab.kind_buttons["absent"].setChecked(True)
    tab._on_kind_clicked("absent")
    tab.kind_buttons["absent"].setChecked(False)
    tab._on_kind_clicked("absent")

    assert tab.kind_buttons["All"].isChecked()
    assert tab.model.rowCount() == len(tab._events)


def test_exporting_writes_the_filtered_rows_as_csv(state, qtbot, result, chroms, tmp_path):
    page = review(state, qtbot, result)
    tab = page.ensure_tab(TAB_QAQC)
    tab.kind_buttons["recovery"].setChecked(True)
    tab._on_kind_clicked("recovery")

    written = tab.export_csv(tmp_path / "flags.csv")

    lines = written.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0] == "file,analyte,flag,value"
    assert len(lines) == 2
    assert lines[1].startswith("S2,13C 52,recovery,")


# --------------------------------------------------------------------------------------
# The page itself
# --------------------------------------------------------------------------------------


def test_tabs_are_built_only_when_they_are_opened(state, qtbot, result, chroms):
    page = review(state, qtbot, result)
    assert sorted(page._tabs) == [TAB_OVERVIEW], "opening Review builds the verdict, nothing else"

    page.tabs.setCurrentIndex(TAB_RECOVERIES)

    assert sorted(page._tabs) == [TAB_OVERVIEW, TAB_RECOVERIES]


def test_a_new_result_discards_the_tabs_built_for_the_old_one(state, qtbot, result, chroms, tmp_path):
    page = review(state, qtbot, result)
    page.tabs.setCurrentIndex(TAB_SAMPLES)
    stale_tab = page.ensure_tab(TAB_SAMPLES)
    assert page.ensure_tab(TAB_OVERVIEW).banner.tone == "warning"

    state.set_result(build_result(tmp_path, bad_recovery=False, flagged=False))
    page.refresh()

    assert page._tabs[TAB_SAMPLES] is not stale_tab, "no tab survives its own run"
    assert TAB_OVERVIEW not in page._tabs, "an unopened tab is not rebuilt eagerly either"
    assert page.ensure_tab(TAB_OVERVIEW).banner.tone == "success"


def test_no_result_shows_the_empty_state(state, qtbot):
    page = review(state, qtbot)

    assert page.empty.isVisibleTo(page)
    assert not page.tabs.isVisibleTo(page)
    assert page.page_status() is None


# --------------------------------------------------------------------------------------
# Export
# --------------------------------------------------------------------------------------


@pytest.fixture
def export(state, qtbot, result) -> ExportPage:
    state.set_result(result)
    page = ExportPage(state)
    qtbot.addWidget(page)
    page.on_activated()
    return page


def test_the_output_folder_defaults_beside_the_run_and_is_remembered(export, state, tmp_path):
    assert export.output_dir() == Path(state.run_folder) / DEFAULT_OUT_DIRNAME
    assert export.report_path().name == f"run{REPORT_SUFFIX}"

    export.set_output_dir(tmp_path / "elsewhere")

    assert export.output_dir() == tmp_path / "elsewhere"
    assert str(tmp_path / "elsewhere") in export.folder_label.text()
    assert ExportPage(state).output_dir() == tmp_path / "elsewhere", "persisted"


def test_writing_the_excel_report_calls_the_writer_and_shows_the_path(
    export, state, qtbot, monkeypatch, tmp_path
):
    calls: dict[str, object] = {}

    def write_report(result, out_path, run_name=None):
        calls.update(result=result, out_path=Path(out_path), run_name=run_name)
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        Path(out_path).write_bytes(b"xlsx")
        return Path(out_path)

    monkeypatch.setattr("quant209.report.excel.write_report", write_report)

    export.report_card.button.click()
    qtbot.waitUntil(lambda: state.busy is False, timeout=10_000)

    assert calls["out_path"] == export.report_path()
    assert calls["run_name"] == "run"
    assert export.report_card.banner.tone == "success"
    assert str(export.report_path()) in export.report_card.banner.body_label.text()
    assert export.report_card.open_button.isVisibleTo(export.report_card)
    assert state.export_path == export.report_path()
    assert export.page_status() == ("Written", "success")


def test_a_failing_writer_is_reported_in_the_card_not_a_dialog(
    export, state, qtbot, monkeypatch
):
    monkeypatch.setattr(
        "PySide6.QtWidgets.QMessageBox.critical",
        lambda *a, **k: pytest.fail("a modal dialog was opened for an in-page error"),
    )

    def explode(*args, **kwargs):
        raise OSError("the output folder is read-only")

    monkeypatch.setattr("quant209.report.excel.write_report", explode)

    export.report_card.button.click()
    qtbot.waitUntil(lambda: state.busy is False, timeout=10_000)

    assert export.report_card.banner.tone == "error"
    assert "read-only" in export.report_card.banner.body_label.text()
    assert state.export_path is None


def test_the_validation_button_is_disabled_and_explains_why_without_a_workbook(export):
    assert not export.validation_card.button.isEnabled()
    text = export.validation_card.description.text()
    assert "no analyst quantification workbook was found" in text
    assert "run" in text


def test_the_validation_button_unlocks_when_a_workbook_is_present(export, state):
    (Path(state.run_folder) / "SMP16-20_20231003.xlsx").write_bytes(b"x")
    export.refresh()

    assert export.validation_card.button.isEnabled()
    assert "SMP16-20_20231003.xlsx" in export.validation_card.description.text()


def test_writing_the_validation_report_reports_the_verdict_and_the_folder(
    export, state, qtbot, monkeypatch
):
    (Path(state.run_folder) / "SMP16-20_20231003.xlsx").write_bytes(b"x")
    export.refresh()
    calls: dict[str, object] = {}

    monkeypatch.setattr(
        "quant209.report.truth.read_run_truth",
        lambda folder: SimpleNamespace(areas=SimpleNamespace(empty=False)),
    )

    def generate(result, truth, out_dir, run_name, **kwargs):
        calls.update(out_dir=Path(out_dir), run_name=run_name)
        return SimpleNamespace(passed=True)

    monkeypatch.setattr("quant209.report.validation.generate_validation_report", generate)

    export.validation_card.button.click()
    qtbot.waitUntil(lambda: state.busy is False, timeout=30_000)

    assert calls["out_dir"] == export.output_dir() / VALIDATION_DIRNAME / "run"
    assert export.validation_card.banner.tone == "success"
    assert "PASS" in export.validation_card.banner.title_label.text()


def test_a_run_folder_with_no_usable_truth_fails_in_the_card(
    export, state, qtbot, monkeypatch
):
    (Path(state.run_folder) / "SMP16-20_20231003.xlsx").write_bytes(b"x")
    export.refresh()
    monkeypatch.setattr(
        "quant209.report.truth.read_run_truth",
        lambda folder: SimpleNamespace(areas=SimpleNamespace(empty=True)),
    )

    export.validation_card.button.click()
    qtbot.waitUntil(lambda: state.busy is False, timeout=30_000)

    assert export.validation_card.banner.tone == "error"
    assert "nothing to validate against" in export.validation_card.banner.body_label.text()


def test_both_buttons_are_gated_while_a_worker_runs(export, state):
    state.set_busy(True)
    export.refresh()

    assert not export.report_card.button.isEnabled()
    assert not export.validation_card.button.isEnabled()
    assert not export.choose_button.isEnabled()
    assert export.page_status() == ("Writing", "accent")


def test_export_shows_an_empty_state_before_a_run(state, qtbot):
    state.set_result(None)
    page = ExportPage(state)
    qtbot.addWidget(page)

    assert page.empty.isVisibleTo(page)
    assert not page.report_card.isVisibleTo(page)
    assert page.page_status() is None


# --------------------------------------------------------------------------------------
# The constants these two pages duplicate rather than import
# --------------------------------------------------------------------------------------


def test_the_recovery_band_matches_the_qaqc_defaults():
    from quant209.core.qaqc import QaqcParams

    assert (rp.RECOVERY_LOW, rp.RECOVERY_HIGH) == (
        QaqcParams().recovery_low,
        QaqcParams().recovery_high,
    )
    assert (charts_module.RECOVERY_LOW, charts_module.RECOVERY_HIGH) == (
        rp.RECOVERY_LOW,
        rp.RECOVERY_HIGH,
    )


def test_the_homolog_surrogates_match_the_standards_table():
    from quant209.core.standards import HOMOLOG_SURROGATE

    assert rp.HOMOLOG_SURROGATE == HOMOLOG_SURROGATE
    assert list(rp.HOMOLOG_SURROGATE) == list(HOMOLOG_SURROGATE), "chlorination order"


def test_the_flag_vocabulary_matches_the_integration_engine():
    from quant209.core.integrate import KNOWN_FLAGS

    assert KNOWN_FLAGS <= set(rp.FLAG_KINDS)
    assert set(rp.QUALITY_FLAGS) <= KNOWN_FLAGS
    assert set(rp.FLAG_TONES) == set(rp.FLAG_KINDS)


def test_the_engines_identifiers_are_written_as_words_on_screen():
    """Stored as ``low_sn`` and ``valley_drop``, read as "low S/N" and "valley drop".

    The kind an event carries is the engine's, and the CSV export and every filter compare
    against it. Only the painted text is translated. A screen that shows an analyst the
    integrator's variable names is a screen no analyst read.
    """
    from quant209.core.integrate import KNOWN_METHODS

    assert rp.flag_text("low_sn") == "low S/N"
    assert rp.flag_text("narrow") == "narrow", "no entry means the kind is already a word"
    assert set(rp.FLAG_LABELS) <= set(rp.FLAG_KINDS)
    assert set(rp.METHOD_LABELS) == KNOWN_METHODS


def test_the_peak_table_prints_the_method_in_words_and_keeps_the_stored_one(
    state, qtbot, result, chroms
):
    page = review(state, qtbot, result)
    model = page.ensure_tab(TAB_SAMPLES).model

    assert model.row_at(0).method == "valley_drop"
    assert model.data(model.index(0, 4)) == "valley drop"


def test_the_canonical_transition_list_is_the_twenty_one_acquired_channels():
    from quant209.core.config import load_canonical_peaks

    transitions = rp.canonical_transitions(load_canonical_peaks())

    assert len(transitions) == 21
    assert sum(len(t.labels) for t in transitions) == 185
    assert [t.group for t in transitions][:3] == ["mono", "13C 3", "di"]
    assert "325.9 → 255.9" in [t for t in transitions if t.group == "penta"][0].text


def test_the_truth_workbook_scan_matches_the_truth_readers(tmp_path: Path):
    from quant209.report.truth import _find_quant_workbooks

    for name in (
        "SMP16-20_20231003.xlsx",
        "SMP21-25_20231003.xlsx",
        "~$SMP16-20_20231003.xlsx",
        "SMP16-30_20231002_sequence.xlsx",
    ):
        (tmp_path / name).write_bytes(b"x")

    assert truth_workbooks(tmp_path) == _find_quant_workbooks(tmp_path)
    assert len(truth_workbooks(tmp_path)) == 2


def test_the_export_paths_match_the_command_line_defaults():
    from quant209.cli import DEFAULT_OUT_DIR

    assert DEFAULT_OUT_DIRNAME == DEFAULT_OUT_DIR
    assert REPORT_SUFFIX == "_report.xlsx"


# --------------------------------------------------------------------------------------
# Cold start: these two pages stay off the heavy import path
# --------------------------------------------------------------------------------------


def test_the_review_and_export_modules_import_no_scientific_stack():
    """Importing the two pages must not reach pyqtgraph, pandas, openpyxl or the pipeline.

    ``quant209.app.charts`` is in the banned list on purpose: it is where pyqtgraph lives,
    and the whole lazy-chart arrangement is worth nothing if importing the Review module
    pulls it in anyway.
    """
    probe = (
        "import sys;"
        "import quant209.app.workflow.review_page, quant209.app.workflow.export_page;"
        "banned=('numpy','pandas','pyqtgraph','matplotlib','openpyxl',"
        "'quant209.app.charts','quant209.core.pipeline','quant209.core.qaqc',"
        "'quant209.core.mzml','quant209.report.excel','quant209.report.truth',"
        "'quant209.report.validation');"
        "print(','.join(m for m in banned if m in sys.modules))"
    )
    proc = subprocess.run(
        [sys.executable, "-c", probe], capture_output=True, text=True, timeout=120
    )

    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "", f"eagerly imported: {proc.stdout.strip()}"


def test_the_window_still_reaches_review_and_export_without_them_being_built(
    qtbot, qapp, state, result
):
    """The page registry must still hand back the two pages, lazily, through the window."""
    from quant209.app.main import MainWindow

    configure_app(qapp)
    state.set_result(result)
    window = MainWindow(state)
    qtbot.addWidget(window)

    assert window.built_steps() == [0]
    assert window.goto_step(STEP_REVIEW)
    assert isinstance(window.page(STEP_REVIEW), ReviewPage)
    assert window.goto_step(STEP_EXPORT)
    assert isinstance(window.page(STEP_EXPORT), ExportPage)
    assert window.header.title_label.text() == window.page(STEP_EXPORT).title


def test_a_file_without_an_internal_standard_is_an_error_not_a_gap(tmp_path):
    """No standard, no mass, and the analyst is told so in red: in the banner, in the
    attention table and as a QA/QC event. Silence would read as a sample with no PCBs."""
    result = build_result(tmp_path, bad_recovery=False, flagged=False, cal_ok=True)
    result.flags["S2"] = ["missing_internal_standard:30D=absent"]

    summary = rp.summarise(result)
    tone, headline, why = rp.verdict(summary)
    assert (tone, headline) == ("error", "Internal standard missing") and "S2" in why

    item = next(i for i in summary.failing if i.what == "internal standard missing")
    assert (item.where, item.tone) == ("S2", "error")
    assert "not in this injection" in item.value and "no mass reported for this file" in item.value

    events = [e for e in rp.flag_events(result) if e.kind == "missing_is"]
    assert [(e.file_id, e.analyte) for e in events] == [("S2", "30D")]


def test_an_excluded_calibration_is_named_in_the_attention_table(tmp_path):
    result = build_result(tmp_path, bad_recovery=False, flagged=False, cal_ok=True)
    result.excluded_calibrations = {"Cal01.mzML": "only 31% of the native analytes were found"}

    item = next(i for i in rp.summarise(result).failing if i.what == "calibration injection excluded")
    assert item.where == "Cal01.mzML" and "anchors no file" in item.value
