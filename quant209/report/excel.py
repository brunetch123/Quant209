"""Write the analyst-format Excel workbook from a run result.

The layout follows a manually quantified workbook address for address, so a reviewer who
navigates those by muscle memory can navigate this one the same way. It was established by
reading a real workbook cell by cell with ``openpyxl``. ``report.truth`` documents the same
layout from the reading side. This module is its writing-side twin and reuses its column
vocabulary, with ``area`` at column P and ``RRF`` at column H of the ``cal`` sheet.

What is mirrored, address for address
--------------------------------------
* **Per-file sheets**, one per :attr:`RunResult.reported_ids`, samples then blanks, keeping
  the analyst's own inclusion of blank sheets. Header row 15, data from row 16, one row per
  analyte in **canonical order** (``result.analytes``, which is
  ``core.config.load_canonical_peaks`` order unless the caller overrode it). That row order
  was verified against a real workbook: row 16 = congener '1', row 20 = congener '4', row 200
  = the last of 185 analytes. Columns: B "mass transition", C label, D "Cal RT", F "RT" (this
  file's measured apex), G and P "Area", H "IS", I "RRF", J "Mass (ng)" (uncorrected), K
  "Corrected Mass", N "Peak", O "RT" (duplicate), Q "Height". Both area columns are written
  because the analyst's sheet duplicates the area into a second N-Q block pasted from the raw
  integration results, and P is the column ``report.truth`` reads as ground truth.
* **cal sheet**: header row 5, data from row 6, same canonical row order. C label, D "RT" (the
  analyte's ``cal_rt``), E "Area" (this run's calibration area), F "IS mass/area", G "mass"
  (the certified calibration mass), H "RRF", J "Peak". That is the exact set
  ``report.truth`` reads, C, E and H, plus the analyst's other columns.
* **Summary sheet**: the Recoveries block at the analyst's own row addressing, header row 3,
  mono at row 4 through deca at row 13 (``core.standards.HOMOLOG_SURROGATE`` order, matching
  both the analyst's block and ``report.truth``'s reader).
* **ref sheet**: the analyst's ``ref`` sheet summary block at its own row addressing, header
  row 2, mono at row 3 through deca at row 12.

Three deliberate differences from the analyst's literal cell contents, documented in place
below and in the ``run_info`` sheet:

1. **Column H ("IS") carries the *normalizing internal-standard label*** ("30D"/"204") on
   every native/surrogate row, rather than the analyst's convention of leaving it blank except
   on the two internal-standard rows themselves (where it holds a raw ``IS_spike/IS_area``
   factor). ``report.truth`` never reads this column, so there is no ground truth to
   contradict, and a label per row is more useful in an automated report than a factor on two
   rows a reader has to already know to look for.
2. **Column K ("Corrected Mass") holds the *uncensored* recovery-corrected mass**, which is
   ``RunResult.masses_corrected`` rather than ``masses_censored``. Inspection of a real
   workbook, both its per-sample sheets and its ``Uncorrected_Mass`` and
   ``SS_Corrected_Mass`` pivots, shows the raw computed mass reported everywhere, with no
   evidence of LOQ censoring anywhere in it. Column K therefore mirrors what the analyst
   shows, and the censored figure goes in a separate, clearly labelled additive column.
3. **The Summary sheet carries no single "REF" column.** An analyst keeps one workbook per
   batch, so one reference injection is unambiguous there. This pipeline brackets each file to
   whichever reference injection preceded it (:func:`core.pipeline.bracket_references`), so a
   run spanning several batches has several references and a single REF column would be a
   fiction.
   Each file's own reference is on the ``run_info`` sheet, and the surrogate masses that define
   100% recovery for each reference are on the ``ref`` sheet, keyed by reference id.

Four kinds of analyst sheet are not reproduced. ``cal_info`` holds nothing but the certified
calibration masses, which are already this workbook's ``cal`` sheet column G. ``QuickView``,
``Uncorrected_Mass`` and ``SS_Corrected_Mass`` are pivots that inspection showed to be
cell-for-cell identical to the per-sample sheets' own J and K columns, and
``RunResult.masses`` and ``masses_corrected`` already are those pivots as DataFrames.
``Profile`` is a total-mass-normalised homolog-fraction pivot, likewise pure arithmetic on
numbers already in this workbook. The numbered template sheets '8', '9' and '10' are the
workbook's unused blank slots, described in ``report.truth``'s docstring. Every one of those
pivots also carried stale ``#REF!`` and ``#DIV/0!`` cells in the real workbook, which is
further reason to derive from ``RunResult`` directly.

One sheet is wholly new. ``run_info`` is an audit trail with no analyst equivalent. It
carries the run identity, ``IntegrationParams``, the RT-alignment verdict and fit report, the
reference pairing policy and each file's own pairing, the LOQ parameters, and every file's
skipped analytes. Its own first cell says it is new.
"""

from __future__ import annotations

from dataclasses import fields
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.worksheet import Worksheet

from quant209 import __version__
from quant209.core.pipeline import RunResult
from quant209.core.quantify import IS_ROLES
from quant209.core.standards import HOMOLOG_SURROGATE, IS_SPIKE_NG, SURROGATE_SPIKE_NG

# --------------------------------------------------------------------------------------
# Per-file sheet: 1-based column indices, mirroring the analyst's own layout (see the
# module docstring). A/E/L/M are the analyst's own blank spacer columns and are left blank
# here too, for the same "same address, same meaning" reason the filled columns matter.
# --------------------------------------------------------------------------------------
COL_TRANSITION = 2  # B: "mass transition"
COL_LABEL = 3  # C: analyte label (no header text, matching the analyst sheet)
COL_CAL_RT = 4  # D: "Cal RT"
COL_RT = 6  # F: "RT" (this file's measured apex)
COL_AREA = 7  # G: "Area"
COL_IS = 8  # H: "IS". See module docstring, difference (1).
COL_RRF = 9  # I: "RRF"
COL_MASS = 10  # J: "Mass (ng)" (uncorrected)
COL_CORRECTED = 11  # K: "Corrected Mass". See module docstring, difference (2).
COL_PEAK = 14  # N: "Peak"
COL_RT_DUP = 15  # O: "RT" (duplicate)
COL_AREA_DUP = 16  # P: "Area", the duplicate report.truth reads as ground truth
COL_HEIGHT = 17  # Q: "Height"
# Additive: a spacer at R, then three columns with no analyst counterpart.
COL_CENSORED = 19  # S
COL_DETECTED = 20  # T
COL_FLAGS = 21  # U

SAMPLE_HEADER_ROW = 15
SAMPLE_DATA_START_ROW = 16

_SAMPLE_HEADERS: dict[int, str] = {
    COL_TRANSITION: "mass transition",
    COL_CAL_RT: "Cal RT",
    COL_RT: "RT",
    COL_AREA: "Area",
    COL_IS: "IS",
    COL_RRF: "RRF",
    COL_MASS: "Mass (ng)",
    COL_CORRECTED: "Corrected Mass",
    COL_PEAK: "Peak",
    COL_RT_DUP: "RT",
    COL_AREA_DUP: "Area",
    COL_HEIGHT: "Height",
    COL_CENSORED: "Censored Mass (ng) [quant209]",
    COL_DETECTED: "Detected [quant209]",
    COL_FLAGS: "Flags [quant209]",
}

# --------------------------------------------------------------------------------------
# cal sheet
# --------------------------------------------------------------------------------------
CAL_COL_TRANSITION = 2  # B
CAL_COL_LABEL = 3  # C, report.truth's ground-truth label column
CAL_COL_RT = 4  # D: "RT" (the analyte's own cal_rt)
CAL_COL_AREA = 5  # E, report.truth's ground-truth area column
CAL_COL_IS_RESPONSE = 6  # F: "IS mass/area"
CAL_COL_MASS = 7  # G: "mass" (certified calibration mass)
CAL_COL_RRF = 8  # H, report.truth's ground-truth RRF column
CAL_COL_PEAK = 10  # J: "Peak"
CAL_COL_INVALID = 12  # L: additive, why this analyte has no RRF

CAL_HEADER_ROW = 5
CAL_DATA_START_ROW = 6

_CAL_HEADERS: dict[int, str] = {
    CAL_COL_TRANSITION: "mass transition",
    CAL_COL_RT: "RT",
    CAL_COL_AREA: "Area",
    CAL_COL_IS_RESPONSE: "IS mass/area",
    CAL_COL_MASS: "mass",
    CAL_COL_RRF: "RRF",
    CAL_COL_PEAK: "Peak",
    CAL_COL_INVALID: "No RRF because... [quant209]",
}

# --------------------------------------------------------------------------------------
# ref sheet: the analyst's summary-block row addressing (header row 2, mono..deca at
# rows 3..12), one column per reference injection instead of per workbook.
# --------------------------------------------------------------------------------------
REF_COL_LABEL = 1
REF_HEADER_ROW = 2
REF_DATA_START_ROW = 3

# --------------------------------------------------------------------------------------
# Summary sheet: the Recoveries block's row addressing (header row 3, mono..deca at
# rows 4..13), one column per reported file instead of the analyst's single-workbook set.
# --------------------------------------------------------------------------------------
SUMMARY_COL_LABEL = 1
SUMMARY_HEADER_ROW = 3
SUMMARY_DATA_START_ROW = 4

_HEADER_FONT = Font(bold=True)

_NOT_REPRODUCED = (
    "Analyst sheets not reproduced here (see quant209.report.excel's module docstring for the "
    "full reasoning): cal_info (redundant with this workbook's cal sheet, column 'mass'); "
    "QuickView, Uncorrected_Mass, SS_Corrected_Mass (pivots verified cell-for-cell identical "
    "to the per-sample sheets' own J/K columns); Profile (a derived homolog-fraction pivot, "
    "likewise pure arithmetic on numbers already here); and the unused numbered template "
    "sheets ('8', '9', '10')."
)


# --------------------------------------------------------------------------------------
# Small shared helpers
# --------------------------------------------------------------------------------------


def _style_header(ws: Worksheet, row: int, cols: range) -> None:
    for c in cols:
        ws.cell(row, c).font = _HEADER_FONT


def _autosize(ws: Worksheet, widths: dict[int, int], default: int = 12, max_col: int = 25) -> None:
    for col in range(1, max_col + 1):
        ws.column_dimensions[get_column_letter(col)].width = widths.get(col, default)


def _value(frame: pd.DataFrame, label: str, column: str) -> float | None:
    """A finite float at ``frame.at[label, column]``, or ``None`` for anything else:
    a missing index or column, NaN, or a non-numeric cell. ``None`` leaves the Excel cell blank
    rather than writing a literal "nan", which openpyxl would otherwise happily do."""
    if label not in frame.index or column not in frame.columns:
        return None
    raw = frame.at[label, column]
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def _transition_text(q1: float, q3: float) -> str:
    return f"{q1:.1f} -> {q3:.1f}"


# --------------------------------------------------------------------------------------
# Per-file sheet
# --------------------------------------------------------------------------------------


def _write_sample_sheet(wb: Workbook, result: RunResult, file_id: str) -> None:
    record = result.files[file_id]
    ws = wb.create_sheet(file_id[:31])
    ws.cell(1, 1, f"{file_id} ({record.file_type})")

    for col, text in _SAMPLE_HEADERS.items():
        ws.cell(SAMPLE_HEADER_ROW, col, text)
    _style_header(ws, SAMPLE_HEADER_ROW, range(COL_TRANSITION, COL_FLAGS + 1))
    ws.freeze_panes = f"A{SAMPLE_DATA_START_ROW}"

    peaks_by_label = {p.label: p for p in record.peaks}

    for i, analyte in enumerate(result.analytes):
        row = SAMPLE_DATA_START_ROW + i
        label = analyte.label
        ws.cell(row, COL_TRANSITION, _transition_text(analyte.q1, analyte.q3))
        ws.cell(row, COL_LABEL, label)
        ws.cell(row, COL_CAL_RT, analyte.cal_rt).number_format = "0.000"

        peak = peaks_by_label.get(label)
        if peak is not None:
            ws.cell(row, COL_RT, peak.rt_apex).number_format = "0.000"
            ws.cell(row, COL_RT_DUP, peak.rt_apex).number_format = "0.000"
            ws.cell(row, COL_HEIGHT, float(peak.height)).number_format = "0"
            ws.cell(row, COL_PEAK, 1)
            if peak.flags:
                ws.cell(row, COL_FLAGS, ",".join(peak.flags))

        area = _value(result.areas, label, file_id)
        if area is not None:
            ws.cell(row, COL_AREA, area).number_format = "0"
            ws.cell(row, COL_AREA_DUP, area).number_format = "0"

        if analyte.role not in IS_ROLES:
            ws.cell(row, COL_IS, analyte.is_norm)

        rrf = result.calibration.rrf.get(label)
        if rrf is not None:
            ws.cell(row, COL_RRF, rrf).number_format = "0.0000"

        mass = _value(result.masses, label, file_id)
        if mass is not None:
            ws.cell(row, COL_MASS, mass).number_format = "0.0000"

        corrected = _value(result.masses_corrected, label, file_id)
        if corrected is not None:
            ws.cell(row, COL_CORRECTED, corrected).number_format = "0.0000"

        censored = _value(result.masses_censored, label, file_id)
        if censored is not None:
            ws.cell(row, COL_CENSORED, censored).number_format = "0.0000"

        if label in result.detected.index and file_id in result.detected.columns:
            ws.cell(row, COL_DETECTED, bool(result.detected.at[label, file_id]))

    _autosize(
        ws,
        {COL_TRANSITION: 16, COL_LABEL: 10, COL_FLAGS: 24, COL_CORRECTED: 16,
         COL_CENSORED: 22, COL_DETECTED: 14},
    )


# --------------------------------------------------------------------------------------
# cal sheet
# --------------------------------------------------------------------------------------


def _write_cal_sheet(wb: Workbook, result: RunResult) -> None:
    ws = wb.create_sheet("cal")
    calibration = result.calibration
    injections = ", ".join(calibration.injections) or "(none)"
    ws.cell(1, 1, f"Calibration model: RRFs averaged over {injections}.")
    ws.cell(
        2, 1,
        "Area / IS mass/area shown are from the first calibration injection "
        f"({calibration.injections[0]})." if calibration.injections
        else "Area / IS mass/area: no usable calibration injection.",
    )
    ws.cell(3, 1, f"Native calibration-area floor applied: {calibration.area_floor:.0f} counts*s.")

    for col, text in _CAL_HEADERS.items():
        ws.cell(CAL_HEADER_ROW, col, text)
    _style_header(ws, CAL_HEADER_ROW, range(CAL_COL_TRANSITION, CAL_COL_INVALID + 1))
    ws.freeze_panes = f"A{CAL_DATA_START_ROW}"

    for i, analyte in enumerate(result.analytes):
        row = CAL_DATA_START_ROW + i
        label = analyte.label
        ws.cell(row, CAL_COL_TRANSITION, _transition_text(analyte.q1, analyte.q3))
        ws.cell(row, CAL_COL_LABEL, label)
        ws.cell(row, CAL_COL_RT, analyte.cal_rt).number_format = "0.000"

        area = calibration.cal_area.get(label)
        if area is not None and np.isfinite(area):
            ws.cell(row, CAL_COL_AREA, float(area)).number_format = "0"

        if analyte.role in IS_ROLES:
            response = calibration.is_response.get(label)
            if response is not None and np.isfinite(response):
                ws.cell(row, CAL_COL_IS_RESPONSE, float(response)).number_format = "0.000000"

        if np.isfinite(analyte.cal_mass):
            ws.cell(row, CAL_COL_MASS, float(analyte.cal_mass)).number_format = "0.0000"

        rrf = calibration.rrf.get(label)
        if rrf is not None:
            ws.cell(row, CAL_COL_RRF, rrf).number_format = "0.0000"
            ws.cell(row, CAL_COL_PEAK, 1)

        reason = calibration.invalid.get(label)
        if reason:
            ws.cell(row, CAL_COL_INVALID, reason)

    _autosize(ws, {CAL_COL_TRANSITION: 16, CAL_COL_LABEL: 10, CAL_COL_INVALID: 44})


# --------------------------------------------------------------------------------------
# ref sheet
# --------------------------------------------------------------------------------------


def _write_ref_sheet(wb: Workbook, result: RunResult) -> None:
    ws = wb.create_sheet("ref")
    ws.cell(
        1, 1,
        "Reference-injection surrogate masses (ng), the basis for recovery correction. One "
        "column per reference injection this run used. The run_info sheet records which "
        "reference each sample and blank was paired with. Row addressing follows the "
        "analyst ref sheet's own summary block, with the header on row 2 and mono to deca "
        "on rows 3-12. The per-analyte table below that block is not reproduced, because it "
        "duplicates data already on this workbook's other sheets.",
    )
    ws.cell(REF_HEADER_ROW, REF_COL_LABEL, "Recoveries")
    ref_ids = result.ref_ids
    for j, ref_id in enumerate(ref_ids):
        ws.cell(REF_HEADER_ROW, REF_COL_LABEL + 1 + j, ref_id)
    _style_header(ws, REF_HEADER_ROW, range(REF_COL_LABEL, REF_COL_LABEL + 1 + len(ref_ids)))
    ws.freeze_panes = f"A{REF_DATA_START_ROW}"

    for i, (homolog, surrogate) in enumerate(HOMOLOG_SURROGATE.items()):
        row = REF_DATA_START_ROW + i
        ws.cell(row, REF_COL_LABEL, homolog)
        for j, ref_id in enumerate(ref_ids):
            mass = _value(result.masses, surrogate, ref_id)
            if mass is not None:
                ws.cell(row, REF_COL_LABEL + 1 + j, mass).number_format = "0.0000"

    _autosize(ws, {REF_COL_LABEL: 14}, default=14)


# --------------------------------------------------------------------------------------
# Summary sheet: the Recoveries block
# --------------------------------------------------------------------------------------


def _write_summary_sheet(wb: Workbook, result: RunResult) -> None:
    ws = wb.create_sheet("Summary")
    ws.cell(
        1, 1,
        "Surrogate recovery by homolog, as a fraction of the bracketing reference "
        "injection's own surrogate mass. There is no single 'REF' column. This run brackets "
        "each file to whichever reference injection preceded it, so there are several "
        "references rather than one. The ref sheet holds their surrogate masses and "
        "run_info holds each file's own pairing.",
    )
    ws.cell(SUMMARY_HEADER_ROW, SUMMARY_COL_LABEL, "Recoveries")
    ids = result.reported_ids
    for j, file_id in enumerate(ids):
        ws.cell(SUMMARY_HEADER_ROW, SUMMARY_COL_LABEL + 1 + j, file_id)
    _style_header(ws, SUMMARY_HEADER_ROW, range(SUMMARY_COL_LABEL, SUMMARY_COL_LABEL + 1 + len(ids)))
    ws.freeze_panes = f"A{SUMMARY_DATA_START_ROW}"

    for i, homolog in enumerate(HOMOLOG_SURROGATE):
        row = SUMMARY_DATA_START_ROW + i
        ws.cell(row, SUMMARY_COL_LABEL, homolog)
        for j, file_id in enumerate(ids):
            value = _value(result.recoveries, homolog, file_id)
            if value is not None:
                ws.cell(row, SUMMARY_COL_LABEL + 1 + j, value).number_format = "0.00%"

    _autosize(ws, {SUMMARY_COL_LABEL: 14}, default=12)


# --------------------------------------------------------------------------------------
# SamplesList_Standards: file roster + the spike-standard amounts
# --------------------------------------------------------------------------------------


def _write_samples_list_sheet(wb: Workbook, result: RunResult) -> None:
    ws = wb.create_sheet("SamplesList_Standards")
    ws.cell(
        1, 1,
        "Run file list. Mirrors the analyst SamplesList_Standards sheet's Batch ID / Sample "
        "Type columns; Sample ID, Material, Deployment and Collection are lab metadata this "
        "pipeline never receives and are not reproduced.",
    )
    header_row = 3
    headers = ["#", "Batch ID", "File type", "Acquired (UTC)", "Reference used", "Flags"]
    for col, text in enumerate(headers, start=1):
        ws.cell(header_row, col, text)
    _style_header(ws, header_row, range(1, len(headers) + 1))

    file_ids = result.reported_ids + [fid for fid in result.cal_ids + result.ref_ids]
    row = header_row
    for i, file_id in enumerate(file_ids, start=1):
        row = header_row + i
        record = result.files[file_id]
        ws.cell(row, 1, i)
        ws.cell(row, 2, file_id)
        ws.cell(row, 3, record.file_type)
        ws.cell(row, 4, record.acquired.isoformat() if record.acquired else None)
        ws.cell(row, 5, result.reference_for.get(file_id, ""))
        ws.cell(row, 6, "; ".join(result.flags.get(file_id, [])))
    ws.freeze_panes = f"A{header_row + 1}"

    standards_row = row + 3
    ws.cell(standards_row, 1, "Standards spiked per injection (ng)")
    ws.cell(standards_row + 1, 1, "Internal standards (d-PCB30 '30D', PCB204 '204')")
    ws.cell(standards_row + 1, 2, IS_SPIKE_NG).number_format = "0.0000"
    ws.cell(standards_row + 2, 1, "13C surrogates (each of the ten homolog surrogates)")
    ws.cell(standards_row + 2, 2, SURROGATE_SPIKE_NG).number_format = "0.0000"

    _autosize(ws, {2: 16, 5: 16, 6: 30}, default=13)


# --------------------------------------------------------------------------------------
# run_info: additive audit-trail sheet, no analyst equivalent
# --------------------------------------------------------------------------------------


def _write_run_info_sheet(wb: Workbook, result: RunResult, run_name: str) -> None:
    ws = wb.create_sheet("run_info")
    ws.cell(
        1, 1,
        "quant209 run_info, the audit trail for this run. This sheet is new and has no "
        "counterpart in a manually quantified workbook.",
    ).font = _HEADER_FONT
    ws.cell(2, 1, _NOT_REPRODUCED)

    kv_rows: list[tuple[str, object]] = [
        ("Run name", run_name),
        ("Generated", datetime.now().isoformat(timespec="seconds")),
        ("quant209 version", __version__),
        ("Run status", result.status),
        ("Status reason", result.status_reason),
        ("Notes", " | ".join(result.notes) if result.notes else ""),
        ("Calibration injections", ", ".join(result.cal_files)),
        ("Reference injections", ", ".join(result.ref_files)),
        ("Reference pairing policy", result.reference_policy),
        ("RT alignment is identity (canonical GC method)", result.alignment.is_identity),
        ("RT-fit guard", result.fit.summary),
        ("LOQ k (mean + k*SD of usable blanks)", result.qaqc.k_sd),
        ("Censor at/below LOQ", result.qaqc.censor_below_loq),
        ("Recovery flag band (low, high)",
         f"{result.qaqc.recovery_low:.0%}-{result.qaqc.recovery_high:.0%}"),
        ("Min usable blanks required for an LOQ", result.qaqc.min_blanks_for_loq),
    ]
    row = 4
    for key, value in kv_rows:
        ws.cell(row, 1, key).font = _HEADER_FONT
        ws.cell(row, 2, value)
        row += 1

    row += 1
    ws.cell(row, 1, "IntegrationParams").font = _HEADER_FONT
    row += 1
    for field in fields(result.params):
        ws.cell(row, 1, field.name)
        ws.cell(row, 2, getattr(result.params, field.name))
        row += 1

    row += 1
    ws.cell(row, 1, "Reference used, per reported file").font = _HEADER_FONT
    row += 1
    ws.cell(row, 1, "file id")
    ws.cell(row, 2, "reference")
    _style_header(ws, row, range(1, 3))
    row += 1
    for file_id in result.reported_ids:
        ws.cell(row, 1, file_id)
        ws.cell(row, 2, result.reference_for.get(file_id, ""))
        row += 1

    row += 1
    ws.cell(row, 1, "Skipped analytes, per file").font = _HEADER_FONT
    row += 1
    ws.cell(row, 1, "file id")
    ws.cell(row, 2, "skipped analytes (no RT anchor / no peak placed)")
    _style_header(ws, row, range(1, 3))
    row += 1
    for file_id, labels in result.skipped_analytes.items():
        if not labels:
            continue
        ws.cell(row, 1, file_id)
        ws.cell(row, 2, ", ".join(labels))
        row += 1

    _autosize(ws, {1: 46, 2: 70}, default=16)


# --------------------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------------------


def write_report(result: RunResult, out_path: str | Path, run_name: str | None = None) -> Path:
    """Write the analyst-format workbook for one run and return the path written.

    ``run_name`` labels the run on the ``run_info`` sheet. It defaults to ``out_path``'s stem
    when not given. One sheet per :attr:`RunResult.reported_ids` (samples then blanks, as the
    analyst's own workbook does), plus ``cal``, ``ref``, ``Summary``, ``SamplesList_Standards``
    and the additive ``run_info`` sheet. See the module docstring for the full layout and
    every deliberate difference from the analyst's own workbook.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    run_name = run_name or out_path.stem

    wb = Workbook()
    wb.remove(wb.active)

    _write_samples_list_sheet(wb, result)
    _write_cal_sheet(wb, result)
    _write_ref_sheet(wb, result)
    for file_id in result.reported_ids:
        _write_sample_sheet(wb, result, file_id)
    _write_summary_sheet(wb, result)
    _write_run_info_sheet(wb, result, run_name)

    wb.save(out_path)
    return out_path
