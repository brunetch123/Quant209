"""Figures 33 to 41: the run table, and how closely automated areas, recoveries and masses track
the analyst's.

Figures 34, 35, 37 and 39 are drawn on *detected* peaks: both areas positive and the automated
peak not flagged absent. That population follows the pipeline's own presence decision, so a
change in what the pipeline calls absent changes what these figures describe.
"""

from __future__ import annotations

import math
import re
from collections import Counter, defaultdict

import numpy as np

from .data import Dataset
from .svg import (ACCENT, BLUE, GOOD, HOMOLOGS, INK, INK2, INK3, MARGINAL, BAD, METHOD_COLORS, MONO, RED,
                  RULE, SANS, WARN, header, swatch, txt)

ROW_ALT = "#F0F2F5"
HOM_SHORT = {"mono": "1", "di": "2", "tri": "3", "tetra": "4", "penta": "5",
             "hexa": "6", "hepta": "7", "octa": "8", "nona": "9", "deca": "10"}
SUPERSCRIPT = {0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴",
               5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹"}

DATE8_RE = re.compile(r"(\d{4})(\d{2})(\d{2})")
DATE6_RE = re.compile(r"_(\d{4})(\d{2})$")


def power_of_ten(exponent: int) -> str:
    """``10`` raised to ``exponent``, written with a superscript."""
    return "10" + "".join(SUPERSCRIPT[int(d)] for d in str(abs(exponent)))


def detected_peaks(ds: Dataset) -> list[dict]:
    """Peaks both parties integrated: positive areas, and not flagged absent."""
    return [r for r in ds.areas if r["manual_f"] > 0 and r["auto_f"] > 0 and "absent" not in r["flags"]]


# ---------------------------------------------------------------------------------------------
# 33  The runs that were scored
# ---------------------------------------------------------------------------------------------
def _run_date(name: str) -> str:
    m8 = DATE8_RE.findall(name)
    if m8:
        y, m, d = m8[-1]
        return f"{y}-{m}-{d}"
    m6 = DATE6_RE.search(name)
    if m6:
        return f"{m6.group(1)}-{m6.group(2)}"
    return ""


def _run_name(name: str) -> str:
    s = DATE8_RE.sub("", name)
    s = DATE6_RE.sub("", s)
    s = re.sub(r"_{2,}", "_", s)
    return s.strip("_")


def fig_33(ds: Dataset) -> tuple[str, str]:
    rows = [r for r in ds.runs if r["status"] == "ok"]
    rows.sort(key=lambda r: _run_date(r["run"]))

    col_spec = [
        ("Run", "start", MONO, 240),
        ("Date", "middle", SANS, 78),
        ("Cal", "end", SANS, 30),
        ("Ref", "end", SANS, 30),
        ("Blk", "end", SANS, 30),
        ("Smp", "end", SANS, 35),
        ("Peaks", "end", SANS, 55),
        ("Time (s)", "end", SANS, 48),
    ]
    PAD_L = 14
    PAD_R = 10
    COL_GAP = 8
    HDR_H = 32
    ROW_H = 19
    TITLE_H = 36

    col_x = []
    x = PAD_L
    for _, _, _, width in col_spec:
        col_x.append(x)
        x += width + COL_GAP
    W = x - COL_GAP + PAD_R
    H = TITLE_H + HDR_H + ROW_H * len(rows) + 12

    def text_x(ci):
        _, anchor, _, width = col_spec[ci]
        if anchor == "start":
            return col_x[ci]
        if anchor == "end":
            return col_x[ci] + width
        return col_x[ci] + width // 2

    svg = header(W, H)
    svg += txt(PAD_L, TITLE_H - 10, "Validation runs", size=16, weight="700")
    total_peaks = sum(int(r["joined_areas"]) for r in rows)
    n_scored = sum(1 for r in rows if int(r["joined_areas"]) > 0)
    svg += txt(W - PAD_R, TITLE_H - 10,
               f"{len(rows)} runs, {n_scored} scored, {total_peaks:,} peak comparisons",
               size=11, color=INK3, anchor="end")

    hy = TITLE_H
    svg += f'<rect x="0" y="{hy}" width="{W}" height="{HDR_H}" fill="{RULE}" fill-opacity="0.4"/>\n'
    for ci, (label, anchor, _, _) in enumerate(col_spec):
        svg += txt(text_x(ci), hy + 21, label, size=11, weight="600", color=INK2, anchor=anchor)
    svg += f'<line x1="0" y1="{hy + HDR_H}" x2="{W}" y2="{hy + HDR_H}" stroke="{RULE}" stroke-width="1"/>\n'

    for ri, r in enumerate(rows):
        ry = TITLE_H + HDR_H + ri * ROW_H
        ty = ry + 14
        if ri % 2 == 1:
            svg += f'<rect x="0" y="{ry}" width="{W}" height="{ROW_H}" fill="{ROW_ALT}"/>\n'

        peaks = int(r["joined_areas"])
        tc = ds.file_types.get(r["run"], {})
        clr = INK if peaks > 0 else INK3

        svg += txt(text_x(0), ty, _run_name(r["run"]), size=10, color=clr, anchor="start", fam=MONO)
        svg += txt(text_x(1), ty, _run_date(r["run"]), size=10, color=clr, anchor="middle")
        svg += txt(text_x(2), ty, tc.get("cal", ""), size=10, color=clr, anchor="end")
        svg += txt(text_x(3), ty, tc.get("ref", ""), size=10, color=clr, anchor="end")
        svg += txt(text_x(4), ty, tc.get("blank", ""), size=10, color=clr, anchor="end")
        svg += txt(text_x(5), ty, tc.get("sample", ""), size=10, color=clr, anchor="end")
        svg += txt(text_x(6), ty, f"{peaks:,}" if peaks else "—", size=10, color=clr, anchor="end")
        svg += txt(text_x(7), ty, f"{float(r['seconds']):.0f}", size=10, color=clr, anchor="end")

    svg += f'<line x1="0" y1="{H - 12}" x2="{W}" y2="{H - 12}" stroke="{RULE}" stroke-width="0.7"/>\n'
    svg += txt(PAD_L, H - 1, "Grey rows: processed but no area truth available for scoring", size=9, color=INK3)
    svg += "</svg>"
    return "33-validation-runs.svg", svg


# ---------------------------------------------------------------------------------------------
# 34  Automated against manual peak area
# ---------------------------------------------------------------------------------------------
def fig_34(ds: Dataset) -> tuple[str, str]:
    detected = detected_peaks(ds)
    manual_vals = np.array([r["manual_f"] for r in detected])
    auto_vals = np.array([r["auto_f"] for r in detected])
    r_pearson = np.corrcoef(np.log10(manual_vals), np.log10(auto_vals))[0, 1]
    median_are = np.median(np.abs(auto_vals - manual_vals) / manual_vals)

    PW, PH = 480, 480
    ML, MR, MT, MB = 62, 20, 50, 52
    pw = PW - ML - MR
    ph = PH - MT - MB

    # The figure was first drawn on 10^2 to 10^8, which held nearly every detected peak while
    # undetected ones were left out by their absent flag. Where a large share of the population
    # sits below that window, open it to 1 count*s rather than leave those peaks off the plot
    # while the statistics above it still count them.
    below = np.sum((manual_vals < 100) | (auto_vals < 100)) / len(detected)
    lo, hi = (0, 8) if below > 0.15 else (2, 8)
    ticks = list(range(lo, hi + 1))
    tick_labels = ["1" if t == 0 else "10" if t == 1 else power_of_ten(t) for t in ticks]

    def sx(v):
        return ML + (v - lo) / (hi - lo) * pw

    def sy(v):
        return MT + ph - (v - lo) / (hi - lo) * ph

    svg = header(PW, PH)
    svg += txt(ML, 22, "Auto vs manual peak area", size=15, weight="700")
    svg += txt(PW - MR, 22, f"r = {r_pearson:.4f}   median |error| = {median_are*100:.1f}%",
               size=11, color=INK3, anchor="end")
    svg += txt(PW - MR, 38, f"{len(detected):,} detected peaks", size=10, color=INK3, anchor="end")

    for i_t, t in enumerate(ticks):
        x = sx(t)
        y = sy(t)
        svg += f'<line x1="{ML}" y1="{y}" x2="{ML+pw}" y2="{y}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += f'<line x1="{x}" y1="{MT}" x2="{x}" y2="{MT+ph}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(ML - 6, y + 4, tick_labels[i_t], size=10, color=INK3, anchor="end")
        svg += txt(x, MT + ph + 14, tick_labels[i_t], size=10, color=INK3, anchor="middle")

    svg += txt(ML + pw // 2, PH - 6, "Manual area", size=12, color=INK2, anchor="middle", weight="600")
    svg += txt(14, MT + ph // 2, "Auto area", size=12, color=INK2, anchor="middle", weight="600", rotate=-90)

    svg += f'<line x1="{ML}" y1="{MT}" x2="{ML}" y2="{MT+ph}" stroke="{INK3}" stroke-width="1"/>\n'
    svg += f'<line x1="{ML}" y1="{MT+ph}" x2="{ML+pw}" y2="{MT+ph}" stroke="{INK3}" stroke-width="1"/>\n'
    svg += (f'<line x1="{sx(lo)}" y1="{sy(lo)}" x2="{sx(hi)}" y2="{sy(hi)}" stroke="{ACCENT}"'
            f' stroke-width="1.5" stroke-dasharray="6,3" opacity="0.6"/>\n')

    # A fixed-seed subsample keeps the file small and the picture the same from run to run.
    rng = np.random.default_rng(42)
    MAX_PTS = 6000
    if len(detected) > MAX_PTS:
        idx = rng.choice(len(detected), MAX_PTS, replace=False)
        sample_pts = [detected[i] for i in sorted(idx)]
    else:
        sample_pts = detected

    for r in sample_pts:
        lm = math.log10(r["manual_f"])
        la = math.log10(r["auto_f"])
        if lo <= lm <= hi and lo <= la <= hi:
            c = METHOD_COLORS.get(r["method"], INK3)
            svg += f'<circle cx="{sx(lm):.1f}" cy="{sy(la):.1f}" r="1.8" fill="{c}" opacity="0.35"/>\n'

    lx, ly = ML + 12, MT + 12
    for i, (method, color) in enumerate(METHOD_COLORS.items()):
        yy = ly + i * 16
        svg += f'<circle cx="{lx}" cy="{yy}" r="4" fill="{color}" opacity="0.7"/>\n'
        svg += txt(lx + 10, yy + 4, method, size=10, color=INK2, fam=MONO)

    svg += "</svg>"
    return "34-area-scatter.svg", svg


# ---------------------------------------------------------------------------------------------
# 35  Distribution of the relative error
# ---------------------------------------------------------------------------------------------
def fig_35(ds: Dataset) -> tuple[str, str]:
    detected = detected_peaks(ds)
    manual_vals = np.array([r["manual_f"] for r in detected])
    auto_vals = np.array([r["auto_f"] for r in detected])
    pct = (auto_vals - manual_vals) / manual_vals * 100

    clip_lo_h, clip_hi_h = -50, 50
    bin_width = 2
    bins = np.arange(clip_lo_h, clip_hi_h + bin_width, bin_width)
    counts, edges = np.histogram(pct, bins=bins)
    max_count = int(counts.max())

    within_5 = np.sum(np.abs(pct) <= 5) / len(pct) * 100
    within_10 = np.sum(np.abs(pct) <= 10) / len(pct) * 100
    within_20 = np.sum(np.abs(pct) <= 20) / len(pct) * 100
    median_pct = np.median(pct)

    PW, PH = 560, 340
    ML, MR, MT, MB = 62, 20, 50, 52
    pw = PW - ML - MR
    ph = PH - MT - MB
    y_max_h = max_count * 1.08

    def hx(v):
        return ML + (v - clip_lo_h) / (clip_hi_h - clip_lo_h) * pw

    def hy(v):
        return MT + ph - v / y_max_h * ph

    svg = header(PW, PH)
    svg += txt(ML, 22, "Relative error distribution", size=15, weight="700")
    svg += txt(PW - MR, 22, f"median {median_pct:+.1f}%   within ±10%: {within_10:.0f}%",
               size=11, color=INK3, anchor="end")
    svg += txt(PW - MR, 38, f"within ±5%: {within_5:.0f}%   within ±20%: {within_20:.0f}%",
               size=10, color=INK3, anchor="end")

    y_step = 2000 if max_count > 8000 else 1000 if max_count > 4000 else 500
    for yt in range(0, int(y_max_h) + 1, y_step):
        y = hy(yt)
        svg += f'<line x1="{ML}" y1="{y}" x2="{ML+pw}" y2="{y}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(ML - 6, y + 4, f"{yt:,}", size=9, color=INK3, anchor="end")

    for xt in range(-50, 51, 10):
        svg += txt(hx(xt), MT + ph + 14, f"{xt}%", size=10, color=INK3, anchor="middle")

    bands = [(5, ACCENT, "0.75"), (10, ACCENT, "0.45"), (20, BLUE, "0.35"), (math.inf, INK3, "0.3")]
    for i, c in enumerate(counts):
        x1 = hx(edges[i])
        x2 = hx(edges[i + 1])
        bw = x2 - x1 - 1
        bh = c / y_max_h * ph
        top = hy(c)
        mid = (edges[i] + edges[i + 1]) / 2
        col, op = next((col, op) for limit, col, op in bands if abs(mid) <= limit)
        svg += f'<rect x="{x1+0.5}" y="{top}" width="{bw}" height="{bh}" fill="{col}" opacity="{op}" rx="1"/>\n'

    zx = hx(0)
    svg += (f'<line x1="{zx}" y1="{MT}" x2="{zx}" y2="{MT+ph}" stroke="{RED}" stroke-width="1"'
            f' stroke-dasharray="4,3" opacity="0.5"/>\n')
    svg += f'<line x1="{ML}" y1="{MT+ph}" x2="{ML+pw}" y2="{MT+ph}" stroke="{INK3}" stroke-width="1"/>\n'
    svg += f'<line x1="{ML}" y1="{MT}" x2="{ML}" y2="{MT+ph}" stroke="{INK3}" stroke-width="1"/>\n'

    svg += txt(ML + pw // 2, PH - 6, "(auto − manual) / manual", size=12, color=INK2, anchor="middle", weight="600")
    svg += txt(14, MT + ph // 2, "Count", size=12, color=INK2, anchor="middle", weight="600", rotate=-90)

    if ds.extras:
        labels = ["within ±5%", "±5–10%", "±10–20%", "beyond ±20%"]
        for i, (label, (_, col, op)) in enumerate(zip(labels, bands)):
            yy = MT + 10 + i * 15
            svg += swatch(ML + 12, yy, col, opacity=op)
            svg += txt(ML + 27, yy + 9, label, size=9, color=INK2)
        beyond = np.sum(np.abs(pct) > clip_hi_h) / len(pct) * 100
        svg += txt(ML + pw - 8, MT + 16, f"{beyond:.0f}% of peaks lie beyond ±{clip_hi_h}% (not drawn)",
                   size=9, color=INK3, anchor="end")

    svg += "</svg>"
    return "35-error-histogram.svg", svg


# ---------------------------------------------------------------------------------------------
# 36  Boundary method by homolog
# ---------------------------------------------------------------------------------------------
def fig_36(ds: Dataset) -> tuple[str, str]:
    method_by_hom = defaultdict(Counter)
    for r in ds.areas:
        method_by_hom[r["homolog"]][r["method"]] += 1

    PW, PH = 560, 340
    ML, MR, MT, MB = 50, 20, 44, 64
    pw = PW - ML - MR
    ph = PH - MT - MB
    n_bars = len(HOMOLOGS)
    bar_w = pw / n_bars * 0.7

    svg = header(PW, PH)
    svg += txt(ML, 22, "Integration method by homolog group", size=15, weight="700")

    for pct_val in [0, 25, 50, 75, 100]:
        y = MT + ph - pct_val / 100 * ph
        svg += f'<line x1="{ML}" y1="{y}" x2="{ML+pw}" y2="{y}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(ML - 6, y + 4, f"{pct_val}%", size=9, color=INK3, anchor="end")

    method_order = ["valley_drop", "deconvolved", "tangent_skim", "rt_window"]
    for i, hom in enumerate(HOMOLOGS):
        cx = ML + (i + 0.5) * pw / n_bars
        bx = cx - bar_w / 2
        total = sum(method_by_hom[hom].values())
        if total == 0:
            continue
        y_bot = MT + ph
        for method in method_order:
            frac = method_by_hom[hom][method] / total
            bh = frac * ph
            y_top = y_bot - bh
            svg += f'<rect x="{bx}" y="{y_top}" width="{bar_w}" height="{bh}" fill="{METHOD_COLORS[method]}" rx="2"/>\n'
            if frac > 0.06:
                svg += txt(cx, (y_top + y_bot) / 2 + 4, f"{frac*100:.0f}%",
                           size=9, color="#FFF", anchor="middle", weight="600")
            y_bot = y_top
        svg += txt(cx, MT + ph + 14, f"Cl{HOM_SHORT[hom]}", size=11, color=INK2, anchor="middle", weight="600")

    legend_y = PH - 14
    legend_methods = list(METHOD_COLORS.items())
    total_w = len(legend_methods) * 120
    lx_start = ML + (pw - total_w) // 2
    for i, (method, color) in enumerate(legend_methods):
        lx = lx_start + i * 125
        svg += f'<rect x="{lx}" y="{legend_y - 8}" width="10" height="10" rx="2" fill="{color}"/>\n'
        svg += txt(lx + 14, legend_y + 1, method, size=9, color=INK2, fam=MONO)

    svg += f'<line x1="{ML}" y1="{MT+ph}" x2="{ML+pw}" y2="{MT+ph}" stroke="{INK3}" stroke-width="1"/>\n'
    svg += f'<line x1="{ML}" y1="{MT}" x2="{ML}" y2="{MT+ph}" stroke="{INK3}" stroke-width="1"/>\n'
    if ds.extras:
        svg += txt(14, MT + ph // 2, "Share of peaks", size=11, color=INK2, anchor="middle", weight="600", rotate=-90)

    svg += "</svg>"
    return "36-method-by-homolog.svg", svg


# ---------------------------------------------------------------------------------------------
# 37  Relative error by homolog
# ---------------------------------------------------------------------------------------------
def fig_37(ds: Dataset) -> tuple[str, str]:
    err_by_hom = defaultdict(list)
    for r in detected_peaks(ds):
        err_by_hom[r["homolog"]].append((r["auto_f"] - r["manual_f"]) / r["manual_f"] * 100)

    PW, PH = 560, 340
    ML, MR, MT, MB = 50, 20, 44, 52
    pw = PW - ML - MR
    ph = PH - MT - MB
    y_lo_b, y_hi_b = -40, 40

    def by(v):
        return MT + ph - (v - y_lo_b) / (y_hi_b - y_lo_b) * ph

    def clamp(v):
        return max(y_lo_b, min(y_hi_b, v))

    svg = header(PW, PH)
    svg += txt(ML, 22, "Relative error by homolog group", size=15, weight="700")
    svg += txt(PW - MR, 22, "box = IQR, whisker = 5th–95th percentile", size=10, color=INK3, anchor="end")

    for yt in range(-40, 41, 10):
        y = by(yt)
        svg += f'<line x1="{ML}" y1="{y}" x2="{ML+pw}" y2="{y}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(ML - 6, y + 4, f"{yt}%", size=9, color=INK3, anchor="end")
    svg += (f'<line x1="{ML}" y1="{by(0)}" x2="{ML+pw}" y2="{by(0)}" stroke="{RED}" stroke-width="0.8"'
            f' stroke-dasharray="4,3" opacity="0.5"/>\n')

    n_bars = len(HOMOLOGS)
    bw = pw / n_bars * 0.5
    cut_off = False
    for i, hom in enumerate(HOMOLOGS):
        vals = np.array(err_by_hom[hom])
        if len(vals) == 0:
            continue
        cx = ML + (i + 0.5) * pw / n_bars
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        p5, p95 = np.percentile(vals, [5, 95])
        cut_off = cut_off or p5 < y_lo_b or p95 > y_hi_b

        svg += f'<line x1="{cx}" y1="{by(clamp(p5))}" x2="{cx}" y2="{by(clamp(p95))}" stroke="{INK3}" stroke-width="1.2"/>\n'
        svg += f'<line x1="{cx-bw*0.3}" y1="{by(clamp(p5))}" x2="{cx+bw*0.3}" y2="{by(clamp(p5))}" stroke="{INK3}" stroke-width="1.2"/>\n'
        svg += f'<line x1="{cx-bw*0.3}" y1="{by(clamp(p95))}" x2="{cx+bw*0.3}" y2="{by(clamp(p95))}" stroke="{INK3}" stroke-width="1.2"/>\n'
        box_top = by(clamp(q3))
        box_bot = by(clamp(q1))
        svg += (f'<rect x="{cx-bw/2}" y="{box_top}" width="{bw}" height="{box_bot-box_top}" fill="{ACCENT}"'
                f' fill-opacity="0.15" stroke="{ACCENT}" stroke-width="1.5" rx="3"/>\n')
        svg += f'<line x1="{cx-bw/2}" y1="{by(clamp(med))}" x2="{cx+bw/2}" y2="{by(clamp(med))}" stroke="{ACCENT}" stroke-width="2.5"/>\n'

        svg += txt(cx, MT + ph + 14, f"Cl{HOM_SHORT[hom]}", size=11, color=INK2, anchor="middle", weight="600")
        svg += txt(cx, MT + ph + 27, f"n={len(vals):,}", size=8, color=INK3, anchor="middle")

    svg += f'<line x1="{ML}" y1="{MT+ph}" x2="{ML+pw}" y2="{MT+ph}" stroke="{INK3}" stroke-width="1"/>\n'
    svg += f'<line x1="{ML}" y1="{MT}" x2="{ML}" y2="{MT+ph}" stroke="{INK3}" stroke-width="1"/>\n'
    svg += txt(14, MT + ph // 2, "(auto − manual) / manual", size=11, color=INK2, anchor="middle", weight="600", rotate=-90)
    if ds.extras and cut_off:
        svg += txt(PW - MR, 36, f"boxes and whiskers are cut off at ±{y_hi_b}%", size=10, color=INK3, anchor="end")

    svg += "</svg>"
    return "37-error-by-homolog.svg", svg


# ---------------------------------------------------------------------------------------------
# 38  Whether the two parties agree a peak is there
# ---------------------------------------------------------------------------------------------
def detection_counts(ds: Dataset) -> tuple[int, int, int, int]:
    """True positives, false positives, false negatives and true negatives, analyst as reference."""
    tp = fp = fn = tn = 0
    for r in ds.areas:
        if r["role"] != "native" or r["file_type"] != "sample":
            continue
        analyst_present = r["manual_f"] > 0
        auto_present = "absent" not in r["flags"]
        if analyst_present and auto_present:
            tp += 1
        elif not analyst_present and not auto_present:
            tn += 1
        elif not analyst_present and auto_present:
            fp += 1
        else:
            fn += 1
    return tp, fp, fn, tn


def fig_38(ds: Dataset) -> tuple[str, str]:
    tp, fp, fn, tn = detection_counts(ds)
    total_det = tp + fp + fn + tn
    sens = tp / (tp + fn) * 100 if (tp + fn) else 0
    spec = tn / (tn + fp) * 100 if (tn + fp) else 0
    acc = (tp + tn) / total_det * 100 if total_det else 0

    PW, PH = 400, 280
    svg = header(PW, PH)
    svg += txt(20, 28, "Detection agreement", size=15, weight="700")
    svg += txt(PW - 20, 28, f"{total_det:,} native sample peaks", size=11, color=INK3, anchor="end")

    mx, my = 140, 60
    cw, ch = 110, 70
    svg += txt(mx + cw // 2, my - 18, "Analyst reported", size=11, color=INK2, anchor="middle", weight="600")
    svg += txt(mx + cw // 2, my - 6, "present", size=10, color=INK3, anchor="middle")
    svg += txt(mx + cw + cw // 2, my - 6, "absent", size=10, color=INK3, anchor="middle")
    svg += txt(mx - 10, my + ch // 2 + 4, "present", size=10, color=INK3, anchor="end")
    svg += txt(mx - 10, my + ch + ch // 2 + 4, "absent", size=10, color=INK3, anchor="end")
    svg += txt(mx - 28, my + ch, "Auto", size=11, color=INK2, anchor="end", weight="600", rotate=-90)

    cells = [(0, 0, tp, ACCENT, "TP"), (1, 0, fp, WARN, "FP"), (0, 1, fn, RED, "FN"), (1, 1, tn, BLUE, "TN")]
    for cx_i, cy_i, val, color, label in cells:
        rx_c = mx + cx_i * cw
        ry_c = my + cy_i * ch
        opacity = "0.15" if label in ("TP", "TN") else "0.08"
        svg += (f'<rect x="{rx_c}" y="{ry_c}" width="{cw}" height="{ch}" fill="{color}" fill-opacity="{opacity}"'
                f' stroke="{RULE}" stroke-width="1" rx="4"/>\n')
        svg += txt(rx_c + cw // 2, ry_c + ch // 2 - 2, f"{val:,}", size=18, weight="700", color=color, anchor="middle")
        share = val / total_det * 100 if total_det else 0
        svg += txt(rx_c + cw // 2, ry_c + ch // 2 + 16, f"({share:.1f}%)", size=10, color=INK3, anchor="middle")

    sy_stat = my + 2 * ch + 22
    svg += txt(mx - 10, sy_stat, f"Sensitivity {sens:.1f}%", size=12, color=ACCENT, weight="600")
    svg += txt(mx + cw + 10, sy_stat, f"Specificity {spec:.1f}%", size=12, color=BLUE, weight="600")
    svg += txt(mx - 10, sy_stat + 18, f"Overall agreement {acc:.1f}%", size=11, color=INK2, weight="600")

    svg += "</svg>"
    return "38-detection-agreement.svg", svg


# ---------------------------------------------------------------------------------------------
# 39  How well each congener's areas correlate
# ---------------------------------------------------------------------------------------------
def fig_39(ds: Dataset) -> tuple[str, str]:
    by_label = defaultdict(lambda: {"m": [], "a": []})
    for r in detected_peaks(ds):
        by_label[r["label"]]["m"].append(r["manual_f"])
        by_label[r["label"]]["a"].append(r["auto_f"])

    r2_vals = []
    for data in by_label.values():
        if len(data["m"]) < 5:
            continue
        r_c = np.corrcoef(np.log10(np.array(data["m"])), np.log10(np.array(data["a"])))[0, 1]
        r2_vals.append(r_c ** 2)
    r2_arr = np.array(r2_vals)
    below_80 = np.sum(r2_arr < 0.80)

    PW, PH = 480, 300
    ML, MR, MT, MB = 50, 20, 44, 52
    pw = PW - ML - MR
    ph = PH - MT - MB

    x_lo_r, x_hi_r = 0.80, 1.005
    bin_w_r2 = 0.005
    r2_bins = np.arange(x_lo_r, x_hi_r + bin_w_r2, bin_w_r2)
    r2_counts, r2_edges = np.histogram(r2_arr[r2_arr >= 0.80], bins=r2_bins)
    r2_max_count = max(int(r2_counts.max()), 1)

    def rx(v):
        return ML + (v - x_lo_r) / (x_hi_r - x_lo_r) * pw

    def ry_r(v):
        return MT + ph - v / (r2_max_count * 1.1) * ph

    svg = header(PW, PH)
    above_95 = np.sum(r2_arr > 0.95)
    svg += txt(ML, 22, "Per-congener r² distribution", size=15, weight="700")
    svg += txt(PW - MR, 22, f"median = {np.median(r2_arr):.4f}   {above_95}/{len(r2_arr)} above 0.95",
               size=11, color=INK3, anchor="end")
    if below_80 > 0:
        svg += txt(PW - MR, 38, f"{below_80} congener{'s' if below_80 > 1 else ''} below 0.80 (not shown)",
                   size=10, color=INK3, anchor="end")

    # Tens suit the usual spread; a flatter histogram needs finer ticks to have any at all.
    if r2_max_count >= 10:
        y_step_r2 = 10 if r2_max_count < 60 else 20
        y_ticks = range(0, r2_max_count + y_step_r2, y_step_r2)
    else:
        y_ticks = range(0, r2_max_count + 1, 1 if r2_max_count <= 4 else 2)
    for yt in y_ticks:
        y = ry_r(yt)
        svg += f'<line x1="{ML}" y1="{y}" x2="{ML+pw}" y2="{y}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(ML - 6, y + 4, str(yt), size=9, color=INK3, anchor="end")

    for xt_v in [0.80, 0.85, 0.90, 0.95, 1.00]:
        svg += txt(rx(xt_v), MT + ph + 14, f"{xt_v:.2f}", size=10, color=INK3, anchor="middle")

    bands = [(0.99, ACCENT, "0.7"), (0.95, BLUE, "0.5"), (-math.inf, INK3, "0.3")]
    for i, c in enumerate(r2_counts):
        x1 = rx(r2_edges[i])
        x2 = rx(r2_edges[i + 1])
        bw = x2 - x1 - 1
        bh = c / (r2_max_count * 1.1) * ph
        top = ry_r(c)
        mid_r2 = (r2_edges[i] + r2_edges[i + 1]) / 2
        col, op = next((col, op) for limit, col, op in bands if mid_r2 >= limit)
        svg += f'<rect x="{x1+0.5}" y="{top}" width="{bw}" height="{bh}" fill="{col}" opacity="{op}" rx="1"/>\n'

    svg += f'<line x1="{ML}" y1="{MT+ph}" x2="{ML+pw}" y2="{MT+ph}" stroke="{INK3}" stroke-width="1"/>\n'
    svg += f'<line x1="{ML}" y1="{MT}" x2="{ML}" y2="{MT+ph}" stroke="{INK3}" stroke-width="1"/>\n'
    svg += txt(ML + pw // 2, PH - 6, "r² (log-scale areas)", size=12, color=INK2, anchor="middle", weight="600")
    svg += txt(14, MT + ph // 2, "Congeners", size=12, color=INK2, anchor="middle", weight="600", rotate=-90)

    if ds.extras:
        labels = ["r² ≥ 0.99", "0.95 – 0.99", "below 0.95"]
        for i, (label, (_, col, op)) in enumerate(zip(labels, bands)):
            yy = MT + 8 + i * 15
            svg += swatch(ML + 12, yy, col, opacity=op)
            svg += txt(ML + 27, yy + 9, label, size=9, color=INK2)

    svg += "</svg>"
    return "39-congener-r2.svg", svg


# ---------------------------------------------------------------------------------------------
# 40  Surrogate recovery
# ---------------------------------------------------------------------------------------------
def fig_40(ds: Dataset) -> tuple[str, str]:
    rec_pairs = []
    for r in ds.recoveries:
        m = float(r["manual"]) if r["manual"] else None
        a = float(r["auto"]) if r["auto"] else None
        if m is not None and a is not None and 0 < m < 2 and 0 < a < 2:
            rec_pairs.append((m * 100, a * 100, r["homolog"]))

    rec_m = np.array([p[0] for p in rec_pairs])
    rec_a = np.array([p[1] for p in rec_pairs])
    rec_median_err = np.median(np.abs(rec_a - rec_m))
    rec_corr = np.corrcoef(rec_m, rec_a)[0, 1]

    PW, PH = 420, 420
    ML, MR, MT, MB = 52, 20, 50, 52
    pw = PW - ML - MR
    ph = PH - MT - MB
    ax_lo_r, ax_hi_r = 20, 140

    def scx(v):
        return ML + (v - ax_lo_r) / (ax_hi_r - ax_lo_r) * pw

    def scy(v):
        return MT + ph - (v - ax_lo_r) / (ax_hi_r - ax_lo_r) * ph

    svg = header(PW, PH)
    svg += txt(ML, 22, "Surrogate recovery", size=15, weight="700")
    svg += txt(PW - MR, 22, f"r = {rec_corr:.3f}   median |Δ| = {rec_median_err:.1f} pp",
               size=11, color=INK3, anchor="end")
    svg += txt(PW - MR, 38, f"{len(rec_pairs):,} recovery values", size=10, color=INK3, anchor="end")

    for t_val in range(20, 141, 20):
        x = scx(t_val)
        y = scy(t_val)
        svg += f'<line x1="{ML}" y1="{y}" x2="{ML+pw}" y2="{y}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += f'<line x1="{x}" y1="{MT}" x2="{x}" y2="{MT+ph}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(ML - 6, y + 4, f"{t_val}%", size=9, color=INK3, anchor="end")
        svg += txt(x, MT + ph + 14, f"{t_val}%", size=9, color=INK3, anchor="middle")

    svg += (f'<line x1="{scx(ax_lo_r)}" y1="{scy(ax_lo_r)}" x2="{scx(ax_hi_r)}" y2="{scy(ax_hi_r)}"'
            f' stroke="{ACCENT}" stroke-width="1.5" stroke-dasharray="6,3" opacity="0.6"/>\n')

    rng2 = np.random.default_rng(99)
    MAX_REC = 4000
    if len(rec_pairs) > MAX_REC:
        idx = rng2.choice(len(rec_pairs), MAX_REC, replace=False)
        plot_recs = [rec_pairs[i] for i in sorted(idx)]
    else:
        plot_recs = rec_pairs

    for m_val, a_val, _hom in plot_recs:
        if ax_lo_r <= m_val <= ax_hi_r and ax_lo_r <= a_val <= ax_hi_r:
            svg += f'<circle cx="{scx(m_val):.1f}" cy="{scy(a_val):.1f}" r="2" fill="{ACCENT}" opacity="0.3"/>\n'

    svg += f'<line x1="{ML}" y1="{MT}" x2="{ML}" y2="{MT+ph}" stroke="{INK3}" stroke-width="1"/>\n'
    svg += f'<line x1="{ML}" y1="{MT+ph}" x2="{ML+pw}" y2="{MT+ph}" stroke="{INK3}" stroke-width="1"/>\n'
    svg += txt(ML + pw // 2, PH - 6, "Manual recovery (%)", size=12, color=INK2, anchor="middle", weight="600")
    svg += txt(14, MT + ph // 2, "Auto recovery (%)", size=12, color=INK2, anchor="middle", weight="600", rotate=-90)

    svg += "</svg>"
    return "40-recovery-scatter.svg", svg


# ---------------------------------------------------------------------------------------------
# 41  Total PCB mass per sample
# ---------------------------------------------------------------------------------------------
def sample_mass_totals(ds: Dataset) -> list[tuple[float, float]]:
    """Manual and automated total PCB mass of every sample where both are positive.

    Native congeners only. The labelled surrogates and internal standards are spiked into every
    file, so they are not PCBs from the sample, and adding them puts about 120 ng under every
    point. Worse, one bad surrogate mass then moves a whole run: in one archive run the program
    put two late surrogates at 8,000 to 20,000 ng each, and every sample in that run appeared
    ten to a hundred times too high while its native total agreed within half.
    """
    totals = defaultdict(lambda: {"manual": 0.0, "auto": 0.0})
    for r in ds.sample_natives:
        key = (r["run"], r["sample"])
        totals[key]["manual"] += r["_m"]
        totals[key]["auto"] += r["_a"]
    return [(v["manual"], v["auto"]) for v in totals.values() if v["manual"] > 0 and v["auto"] > 0]


def fig_41(ds: Dataset) -> tuple[str, str]:
    pairs_mass = sample_mass_totals(ds)
    tm = np.array([p[0] for p in pairs_mass])
    ta = np.array([p[1] for p in pairs_mass])
    log_tm = np.log10(tm)
    log_ta = np.log10(ta)
    r_mass = np.corrcoef(log_tm, log_ta)[0, 1]

    PW, PH = 480, 480
    ML, MR, MT, MB = 62, 20, 50, 52
    pw = PW - ML - MR
    ph = PH - MT - MB

    lo_m = max(0, math.floor(np.min(log_tm) * 2) / 2)
    hi_m = math.ceil(np.max(log_tm) * 2) / 2
    lo_m = min(lo_m, math.floor(np.min(log_ta) * 2) / 2)
    hi_m = max(hi_m, math.ceil(np.max(log_ta) * 2) / 2)
    lo_m = max(0, lo_m - 0.25)
    hi_m = hi_m + 0.25

    def smx(v):
        return ML + (v - lo_m) / (hi_m - lo_m) * pw

    def smy(v):
        return MT + ph - (v - lo_m) / (hi_m - lo_m) * ph

    svg = header(PW, PH)
    svg += txt(ML, 22, "Total PCB mass per sample", size=15, weight="700")
    rpd = np.abs(ta - tm) / ((ta + tm) / 2) * 100
    n_lt10 = int(np.sum(rpd < 10))
    n_10_25 = int(np.sum((rpd >= 10) & (rpd < 25)))
    n_ge25 = int(np.sum(rpd >= 25))
    bands = [(lambda v: v < 10, GOOD, "0.55"), (lambda v: 10 <= v < 25, MARGINAL, "0.65"), (lambda v: v >= 25, BAD, "0.70")]

    svg += txt(PW - MR, 22, f"r = {r_mass:.4f}   median RPD = {np.median(rpd):.1f}%", size=11, color=INK3, anchor="end")
    svg += txt(PW - MR, 38, f"{len(pairs_mass)} samples, sum of native congeners (ng)", size=10, color=INK3, anchor="end")

    for exp in range(int(math.ceil(lo_m)), int(math.floor(hi_m)) + 1):
        x = smx(exp)
        y = smy(exp)
        svg += f'<line x1="{ML}" y1="{y}" x2="{ML+pw}" y2="{y}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += f'<line x1="{x}" y1="{MT}" x2="{x}" y2="{MT+ph}" stroke="{RULE}" stroke-width="0.5"/>\n'
        label = power_of_ten(exp)
        svg += txt(ML - 6, y + 4, label, size=10, color=INK3, anchor="end")
        svg += txt(x, MT + ph + 14, label, size=10, color=INK3, anchor="middle")

    svg += txt(ML + pw // 2, PH - 6, "Manual total mass (ng)", size=12, color=INK2, anchor="middle", weight="600")
    svg += txt(14, MT + ph // 2, "Auto total mass (ng)", size=12, color=INK2, anchor="middle", weight="600", rotate=-90)
    svg += f'<line x1="{ML}" y1="{MT}" x2="{ML}" y2="{MT+ph}" stroke="{INK3}" stroke-width="1"/>\n'
    svg += f'<line x1="{ML}" y1="{MT+ph}" x2="{ML+pw}" y2="{MT+ph}" stroke="{INK3}" stroke-width="1"/>\n'
    svg += (f'<line x1="{smx(lo_m)}" y1="{smy(lo_m)}" x2="{smx(hi_m)}" y2="{smy(hi_m)}" stroke="{ACCENT}"'
            f' stroke-width="1.5" stroke-dasharray="6,3" opacity="0.6"/>\n')

    # Worst agreement first, so the well-behaved majority is drawn over the outliers, not under.
    point_data = list(zip(pairs_mass, rpd))
    point_data.sort(key=lambda p: -p[1])
    for (m_v, a_v), rpd_v in point_data:
        lmv = math.log10(m_v)
        lav = math.log10(a_v)
        col, op = next((col, op) for test, col, op in bands if test(rpd_v))
        svg += (f'<circle cx="{smx(lmv):.1f}" cy="{smy(lav):.1f}" r="3.5" fill="{col}" opacity="{op}"'
                f' stroke="{col}" stroke-width="0.5" stroke-opacity="0.6"/>\n')

    band_labels = [f"RPD < 10%  ({n_lt10})", f"RPD 10–25%  ({n_10_25})", f"RPD ≥ 25%  ({n_ge25})"]
    lx_l = ML + pw - 10
    for i in range(3):
        ly_l = MT + ph - 10 - (2 - i) * 18
        svg += f'<circle cx="{lx_l - 130}" cy="{ly_l - 3}" r="5" fill="{bands[i][1]}" opacity="0.7"/>\n'
        svg += txt(lx_l - 120, ly_l + 1, band_labels[i], size=10, color=INK2, fam=MONO)

    svg += "</svg>"
    return "41-total-mass-scatter.svg", svg


FIGURES = [fig_33, fig_34, fig_35, fig_36, fig_37, fig_38, fig_39, fig_40, fig_41]
