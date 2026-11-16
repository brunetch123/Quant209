"""Figures 51 to 58: which boundary method a peak gets, and whether the pipeline detects it at all.

Every figure here is drawn on native analytes in samples. A peak is *detected* when the
pipeline reports it so (``detected`` in the pooled masses); a *missed* peak is one the analyst
gave a mass that the pipeline called absent.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict

from .data import Dataset
from .svg import (BAD, GOOD, HOMOLOG_COLORS, HOMOLOGS, INK, INK2, INK3, MARGINAL, METHOD_COLORS, METHOD_LABELS,
                  MONO, RULE, header, swatch, txt)

METH_ORDER = ["valley_drop", "rt_window", "deconvolved", "tangent_skim"]
BINS = [
    (0, 0.01, "0–0.01"), (0.01, 0.1, "0.01–0.1"), (0.1, 1, "0.1–1"),
    (1, 10, "1–10"), (10, 100, "10–100"), (100, 1000, "100–1k"),
    (1000, 1e6, ">1k"),
]


def _y_axis(ml, mt, pw, ph, title):
    out = ""
    for pct in [0, 25, 50, 75, 100]:
        yy = mt + ph - pct / 100 * ph
        out += f'<line x1="{ml}" y1="{yy}" x2="{ml + pw}" y2="{yy}" stroke="{RULE}" stroke-width="0.5"/>\n'
        out += txt(ml - 6, yy + 4, f"{pct}%", size=9, color=INK3, anchor="end", fam=MONO)
    out += txt(ml - 40, mt + ph / 2, title, size=11, color=INK2, anchor="middle", weight="500", rotate=-90)
    return out


def _method_legend(ml, h):
    out = ""
    leg_y = h - 16
    for i, meth in enumerate(METH_ORDER):
        lx = ml + i * 130
        out += f'<rect x="{lx}" y="{leg_y - 8}" width="12" height="12" fill="{METHOD_COLORS[meth]}" opacity="0.8" rx="2"/>\n'
        out += txt(lx + 16, leg_y + 3, METHOD_LABELS[meth], size=10, color=INK2)
    return out


def _method_stack(bx, bar_w, mt, ph, counts, total, label_size):
    out = ""
    cy = mt + ph
    for meth in METH_ORDER:
        frac = counts.get(meth, 0) / total
        seg_h = frac * ph
        cy -= seg_h
        out += f'<rect x="{bx}" y="{cy}" width="{bar_w}" height="{seg_h}" fill="{METHOD_COLORS[meth]}" opacity="0.8" rx="2"/>\n'
        if frac >= 0.08:
            out += txt(bx + bar_w/2, cy + seg_h/2 + 4, f"{100*frac:.0f}%", size=label_size, weight="600",
                       color="white", anchor="middle")
    return out


# ---------------------------------------------------------------------------------------------
# 51  Boundary method by peak mass
# ---------------------------------------------------------------------------------------------
def fig_51(ds: Dataset) -> tuple[str, str]:
    rows = ds.sample_natives
    binned = [(label, [r for r in rows if r["_m"] > 0 and lo < r["_m"] <= hi]) for lo, hi, label in BINS]
    shares = [Counter(r["method"] for r in in_bin) for _, in_bin in binned if in_bin]
    sizes = [len(in_bin) for _, in_bin in binned if in_bin]

    def share(i, meth):
        return shares[i].get(meth, 0) / sizes[i] * 100

    shifts = bool(shares) and share(0, "rt_window") - share(-1, "rt_window") >= 20 \
        and share(-1, "valley_drop") - share(0, "valley_drop") >= 20
    subtitle = ("method selection shifts from RT window to valley drop as mass increases" if shifts
                else "share of each boundary method in each mass band")

    W, H = 560, 420
    ML, MR, MT, MB = 56, 20, 56, 80
    PW = W - ML - MR
    PH = H - MT - MB
    bar_w = PW / len(BINS) - 8

    svg = header(W, H)
    svg += txt(ML, 24, "Integration method by peak mass", size=14, weight="700")
    svg += txt(ML, 40, subtitle, size=11, color=INK3)

    for bi, (label, in_bin) in enumerate(binned):
        total = len(in_bin)
        if total == 0:
            continue
        bx = ML + bi * (PW / len(BINS)) + 4
        svg += _method_stack(bx, bar_w, MT, PH, Counter(r["method"] for r in in_bin), total, 9)
        svg += txt(bx + bar_w/2, MT + PH + 14, label, size=9, color=INK2, anchor="middle", fam=MONO)
        svg += txt(bx + bar_w/2, MT + PH + 26, f"n={total:,}", size=8, color=INK3, anchor="middle", fam=MONO)

    for pct in [0, 25, 50, 75, 100]:
        yy = MT + PH - pct / 100 * PH
        svg += f'<line x1="{ML}" y1="{yy}" x2="{ML + PW}" y2="{yy}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(ML - 6, yy + 4, f"{pct}%", size=9, color=INK3, anchor="end", fam=MONO)
    svg += txt(ML + PW/2, MT + PH + 42, "manual peak mass (ng)", size=11, color=INK2, anchor="middle", weight="500")
    svg += txt(ML - 40, MT + PH/2, "fraction of peaks", size=11, color=INK2, anchor="middle", weight="500", rotate=-90)
    svg += _method_legend(ML, H)

    svg += "</svg>"
    return "51-method-by-mass.svg", svg


# ---------------------------------------------------------------------------------------------
# 52  Boundary method by homolog
# ---------------------------------------------------------------------------------------------
def fig_52(ds: Dataset) -> tuple[str, str]:
    rows = ds.sample_natives
    by_hom = {h: [r for r in rows if r["homolog"] == h] for h in HOMOLOGS}
    heavy = [by_hom[h] for h in ("nona", "deca") if by_hom[h]]
    window_only = bool(heavy) and all(
        sum(1 for r in hrows if r["method"] == "rt_window") / len(hrows) >= 0.85 for hrows in heavy)
    subtitle = ("heavier homologs rely almost entirely on RT window integration" if window_only
                else "share of each boundary method in each homolog group")

    W, H = 580, 420
    ML, MR, MT, MB = 56, 20, 56, 80
    PW = W - ML - MR
    PH = H - MT - MB
    bar_w = PW / len(HOMOLOGS) - 6

    svg = header(W, H)
    svg += txt(ML, 24, "Integration method by homolog group", size=14, weight="700")
    svg += txt(ML, 40, subtitle, size=11, color=INK3)

    for hi, h in enumerate(HOMOLOGS):
        total = len(by_hom[h])
        if total == 0:
            continue
        bx = ML + hi * (PW / len(HOMOLOGS)) + 3
        svg += _method_stack(bx, bar_w, MT, PH, Counter(r["method"] for r in by_hom[h]), total, 8)
        svg += txt(bx + bar_w/2, MT + PH + 14, h, size=9, color=INK2, anchor="middle", weight="500")

    svg += _y_axis(ML, MT, PW, PH, "fraction of peaks")
    svg += _method_legend(ML, H)

    svg += "</svg>"
    return "52-method-by-homolog.svg", svg


# ---------------------------------------------------------------------------------------------
# 53  Detection rate by peak mass
# ---------------------------------------------------------------------------------------------
RATE_BANDS = [(60, GOOD, "≥ 60%"), (30, MARGINAL, "30–60%"), (-math.inf, BAD, "< 30%")]


def _rate_label(cx, by, bar_h, ph, label_y, rate, color):
    """The percentage on a detection bar.

    It sits above the bar at ``label_y``. A bar that fills the plot has no room above it, so
    its label goes inside. A rate just short of 100 % is not rounded up to it.
    """
    if rate < 99.5 or rate == 100:
        text = f"{rate:.0f}%"
    else:
        text = f"{rate:.1f}%" if rate < 99.95 else ">99.9%"
    if bar_h > ph - 4:
        return txt(cx, by + 14, text, size=10, weight="600", color="white", anchor="middle", fam=MONO)
    return txt(cx, label_y, text, size=10, weight="600", color=color, anchor="middle", fam=MONO)


def fig_53(ds: Dataset) -> tuple[str, str]:
    rows = ds.sample_natives
    det_data = []
    for lo, hi, label in BINS:
        in_bin = [r for r in rows if r["_m"] > 0 and lo < r["_m"] <= hi]
        if not in_bin:
            continue
        det = sum(1 for r in in_bin if r["_det"])
        det_data.append((label, det / len(in_bin) * 100, len(in_bin), det))

    rates = [d[1] for d in det_data]
    falls_at_low_mass = bool(rates) and max(rates) - min(rates) >= 20 and rates[0] == min(rates)
    if falls_at_low_mass:
        subtitle = "low-mass peaks are flagged absent far more often"
    elif rates:
        subtitle = f"at least {math.floor(min(rates))}% of peaks are detected in every mass band"
    else:
        subtitle = "no peaks to score"

    W, H = 520, 400
    ML, MR, MT, MB = 56, 20, 56, 56
    PW = W - ML - MR
    PH = H - MT - MB

    svg = header(W, H)
    svg += txt(ML, 24, "Peak detection rate by mass", size=14, weight="700")
    svg += txt(ML, 40, subtitle, size=11, color=INK3)
    svg += f'<rect x="{ML}" y="{MT}" width="{PW}" height="{PH}" fill="white" rx="3"/>\n'

    bar_w = PW / len(det_data) - 8 if det_data else 0
    for i, (label, rate, total, _det) in enumerate(det_data):
        bx = ML + i * (PW / len(det_data)) + 4
        bar_h = rate / 100 * PH
        by = MT + PH - bar_h
        color = next(col for limit, col, _ in RATE_BANDS if rate >= limit)
        svg += f'<rect x="{bx}" y="{by}" width="{bar_w}" height="{bar_h}" fill="{color}" opacity="0.75" rx="2"/>\n'
        svg += _rate_label(bx + bar_w/2, by, bar_h, PH, by - 6, rate, color)
        svg += txt(bx + bar_w/2, MT + PH + 14, label, size=9, color=INK2, anchor="middle", fam=MONO)
        svg += txt(bx + bar_w/2, MT + PH + 26, f"n={total:,}", size=8, color=INK3, anchor="middle", fam=MONO)

    for pct in [0, 25, 50, 75, 100]:
        yy = MT + PH - pct / 100 * PH
        svg += f'<line x1="{ML}" y1="{yy}" x2="{ML + PW}" y2="{yy}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(ML - 6, yy + 4, f"{pct}%", size=9, color=INK3, anchor="end", fam=MONO)
    svg += txt(ML + PW/2, MT + PH + 42, "manual peak mass (ng)", size=11, color=INK2, anchor="middle", weight="500")
    svg += txt(ML - 40, MT + PH/2, "detected (%)", size=11, color=INK2, anchor="middle", weight="500", rotate=-90)

    if ds.extras:
        lx = W - MR - 196
        svg += txt(lx - 8, 24, "detected", size=9, color=INK3, anchor="end")
        for i, (_, color, label) in enumerate(RATE_BANDS):
            svg += swatch(lx + i * 68, 15, color, opacity="0.75")
            svg += txt(lx + i * 68 + 14, 24, label, size=9, color=INK2)

    svg += "</svg>"
    return "53-detection-by-mass.svg", svg


# ---------------------------------------------------------------------------------------------
# 54  Detection rate by homolog
# ---------------------------------------------------------------------------------------------
def fig_54(ds: Dataset) -> tuple[str, str]:
    rows = ds.sample_natives
    rate_of = {}
    for h in HOMOLOGS:
        hrows = [r for r in rows if r["homolog"] == h]
        if hrows:
            rate_of[h] = sum(1 for r in hrows if r["_det"]) / len(hrows) * 100

    present = [h for h in HOMOLOGS if h in rate_of]
    if not present:
        subtitle = "no peaks to score"
    elif rate_of[present[0]] - rate_of[present[-1]] >= 10:
        first, last = present[0], present[-1]
        subtitle = f"detection drops from {rate_of[first]:.0f}% ({first}) to {rate_of[last]:.0f}% ({last})"
    else:
        subtitle = f"at least {math.floor(min(rate_of.values()))}% of peaks are detected in every homolog group"

    W, H = 560, 400
    ML, MR, MT, MB = 56, 20, 56, 56
    PW = W - ML - MR
    PH = H - MT - MB

    svg = header(W, H)
    svg += txt(ML, 24, "Peak detection rate by homolog", size=14, weight="700")
    svg += txt(ML, 40, subtitle, size=11, color=INK3)
    svg += f'<rect x="{ML}" y="{MT}" width="{PW}" height="{PH}" fill="white" rx="3"/>\n'

    bar_w = PW / len(HOMOLOGS) - 6
    for hi_idx, h in enumerate(HOMOLOGS):
        rate = rate_of.get(h, 0)
        bx = ML + hi_idx * (PW / len(HOMOLOGS)) + 3
        bar_h = rate / 100 * PH
        by = MT + PH - bar_h
        color = HOMOLOG_COLORS[h]
        svg += f'<rect x="{bx}" y="{by}" width="{bar_w}" height="{bar_h}" fill="{color}" opacity="0.75" rx="2"/>\n'
        svg += _rate_label(bx + bar_w/2, by, bar_h, PH, by - 6 if bar_h > 20 else MT + PH - 20, rate, color)
        svg += txt(bx + bar_w/2, MT + PH + 14, h, size=9, color=INK2, anchor="middle", weight="500")

    svg += _y_axis(ML, MT, PW, PH, "detected (%)")

    svg += "</svg>"
    return "54-detection-by-homolog.svg", svg


# ---------------------------------------------------------------------------------------------
# 55  How large are the peaks that are missed?
# ---------------------------------------------------------------------------------------------
def fig_55(ds: Dataset) -> tuple[str, str]:
    rows = ds.sample_natives
    det_true_masses = [r["_m"] for r in rows if r["_det"] and r["_m"] > 0]
    det_false_masses = [r["_m"] for r in rows if not r["_det"] and r["_m"] > 0]

    log_bins = [(10**exp, 10**(exp+1)) for exp in range(-4, 6)]
    hist_det = [sum(1 for m in det_true_masses if lo <= m < hi) for lo, hi in log_bins]
    hist_miss = [sum(1 for m in det_false_masses if lo <= m < hi) for lo, hi in log_bins]
    max_count = max(max(hist_det), max(hist_miss))

    n_missed = len(det_false_masses)
    n_scored = n_missed + len(det_true_masses)
    small = sum(1 for m in det_false_masses if m < 0.1)
    if n_scored and n_missed / n_scored >= 0.05 and small / n_missed >= 0.6:
        subtitle = "missed peaks cluster below 0.1 ng; detected peaks spread across range"
    elif n_scored:
        subtitle = f"{n_missed:,} of {n_scored:,} peaks are missed ({n_missed / n_scored * 100:.1f}%)"
    else:
        subtitle = "no peaks to score"

    W, H = 540, 400
    ML, MR, MT, MB = 56, 20, 56, 64
    PW = W - ML - MR
    PH = H - MT - MB

    svg = header(W, H)
    svg += txt(ML, 24, "Mass distribution: detected vs missed peaks", size=14, weight="700")
    svg += txt(ML, 40, subtitle, size=11, color=INK3)
    svg += f'<rect x="{ML}" y="{MT}" width="{PW}" height="{PH}" fill="white" rx="3"/>\n'

    group_w = PW / len(log_bins)
    half_bar = (group_w - 6) / 2
    for i, (lo, _hi) in enumerate(log_bins):
        gx = ML + i * group_w
        bh_d = hist_det[i] / max_count * PH if max_count else 0
        svg += f'<rect x="{gx + 2}" y="{MT + PH - bh_d}" width="{half_bar}" height="{bh_d}" fill="#0E6B68" opacity="0.7" rx="1"/>\n'
        bh_m = hist_miss[i] / max_count * PH if max_count else 0
        svg += f'<rect x="{gx + 2 + half_bar + 1}" y="{MT + PH - bh_m}" width="{half_bar}" height="{bh_m}" fill="#C43D3D" opacity="0.7" rx="1"/>\n'
        lbl = f"{int(lo)}" if lo >= 1 else f"{lo:.0e}".replace("e-0", "e-")
        svg += txt(gx + group_w/2, MT + PH + 14, lbl, size=8, color=INK3, anchor="middle", fam=MONO)

    for frac in [0, 0.25, 0.5, 0.75, 1.0]:
        yy = MT + PH - frac * PH
        svg += f'<line x1="{ML}" y1="{yy}" x2="{ML + PW}" y2="{yy}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(ML - 6, yy + 4, f"{int(frac * max_count):,}", size=8, color=INK3, anchor="end", fam=MONO)

    svg += txt(ML + PW/2, MT + PH + 30, "manual peak mass (ng, log scale)", size=11, color=INK2, anchor="middle", weight="500")
    svg += txt(ML - 42, MT + PH/2, "peak count", size=11, color=INK2, anchor="middle", weight="500", rotate=-90)

    leg_y = H - 14
    svg += f'<rect x="{ML}" y="{leg_y - 8}" width="12" height="12" fill="#0E6B68" opacity="0.7" rx="2"/>\n'
    svg += txt(ML + 16, leg_y + 3, f"detected  (n={len(det_true_masses):,})", size=10, color=INK2)
    svg += f'<rect x="{ML + 200}" y="{leg_y - 8}" width="12" height="12" fill="#C43D3D" opacity="0.7" rx="2"/>\n'
    svg += txt(ML + 216, leg_y + 3, f"missed  (n={n_missed:,})", size=10, color=INK2)

    svg += "</svg>"
    return "55-missed-mass-distribution.svg", svg


# ---------------------------------------------------------------------------------------------
# 56  The congeners most often missed
# ---------------------------------------------------------------------------------------------
def fig_56(ds: Dataset) -> tuple[str, str]:
    cong_stats = defaultdict(lambda: {"total": 0, "missed": 0, "homolog": ""})
    for r in ds.sample_natives:
        if r["_m"] > 0:
            cs = cong_stats[r["label"]]
            cs["total"] += 1
            cs["homolog"] = r["homolog"]
            if not r["_det"]:
                cs["missed"] += 1

    ranked = sorted(cong_stats.items(), key=lambda x: x[1]["missed"], reverse=True)
    by_miss = [(label, stats) for label, stats in ranked[:20] if stats["missed"] > 0]
    max_miss = by_miss[0][1]["missed"] if by_miss else 1

    ROWS_FULL = 20
    W = 560
    ML, MR, MT, MB = 90, 20, 50, 36
    PW = W - ML - MR
    row_h = (590 - MT - MB) / ROWS_FULL
    H = 590 if len(by_miss) == ROWS_FULL else round(MT + max(len(by_miss), 1) * row_h + MB)

    if len(by_miss) == ROWS_FULL:
        title = "Top 20 congeners by missed peak count"
    elif by_miss:
        title = f"The {len(by_miss)} congeners with missed peaks"
    else:
        title = "No congener has a missed peak"

    svg = header(W, H)
    svg += txt(ML, 24, title, size=14, weight="700")
    svg += txt(ML, 40, "detected = False despite manual mass > 0", size=11, color=INK3)

    for i, (label, stats) in enumerate(by_miss):
        yy = MT + i * row_h
        color = HOMOLOG_COLORS.get(stats["homolog"], INK3)
        bar_w_px = stats["missed"] / max_miss * (PW - 120)
        svg += f'<rect x="{ML}" y="{yy + 2}" width="{bar_w_px}" height="{row_h - 4}" fill="{color}" opacity="0.7" rx="2"/>\n'
        svg += txt(ML - 6, yy + row_h/2 + 4, f"PCB-{label}", size=9, weight="500", color=INK, anchor="end", fam=MONO)
        svg += txt(ML + bar_w_px + 6, yy + row_h/2 + 4,
                   f"{stats['missed']}/{stats['total']} ({100*stats['missed']/stats['total']:.0f}%)",
                   size=9, color=INK2, fam=MONO)

    shown = {stats["homolog"] for _, stats in by_miss}
    leg_y = H - 12
    for i, h in enumerate(h for h in HOMOLOGS if h in shown):
        lx = ML + i * 56
        svg += f'<rect x="{lx}" y="{leg_y - 8}" width="10" height="10" fill="{HOMOLOG_COLORS[h]}" opacity="0.7" rx="1"/>\n'
        svg += txt(lx + 14, leg_y + 2, h, size=9, color=INK2)

    svg += "</svg>"
    return "56-congener-miss-rate.svg", svg


# ---------------------------------------------------------------------------------------------
# 57  Which congeners need each of the harder boundary methods
# ---------------------------------------------------------------------------------------------
def fig_57(ds: Dataset) -> tuple[str, str]:
    rows = ds.sample_natives
    total_by_cong = Counter(r["label"] for r in rows)
    hom_by_cong = {r["label"]: r["homolog"] for r in rows}

    W, H = 740, 520
    ML, MR, MT, MB = 80, 20, 50, 20
    CELL_W = (W - ML - MR) / 3
    CELL_H = H - MT - MB
    N_SHOW = 12

    svg = header(W, H)
    svg += txt(ML, 24, "Top congeners by method usage", size=14, weight="700")
    svg += txt(ML, 40, "which congeners most frequently require each integration method", size=11, color=INK3)

    shown = set()
    for ci, meth in enumerate(["deconvolved", "tangent_skim", "rt_window"]):
        cx = ML + ci * CELL_W
        top = Counter(r["label"] for r in rows if r["method"] == meth).most_common(N_SHOW)
        if not top:
            continue
        max_count = top[0][1]
        svg += txt(cx + CELL_W/2, MT + 6, METHOD_LABELS[meth], size=12, weight="700", color=METHOD_COLORS[meth], anchor="middle")

        rh = (CELL_H - 20) / N_SHOW
        for ri, (label, count) in enumerate(top):
            ry = MT + 16 + ri * rh
            bw = count / max_count * (CELL_W - 90)
            hom = hom_by_cong.get(label, "")
            shown.add(hom)
            svg += f'<rect x="{cx + 46}" y="{ry + 2}" width="{bw}" height="{rh - 4}" fill="{HOMOLOG_COLORS.get(hom, INK3)}" opacity="0.65" rx="2"/>\n'
            svg += txt(cx + 42, ry + rh/2 + 3, label, size=8, color=INK2, anchor="end", fam=MONO)
            pct = 100 * count / total_by_cong[label]
            svg += txt(cx + 48 + bw, ry + rh/2 + 3, f"{count} ({pct:.0f}%)", size=8, color=INK3, fam=MONO)

    if ds.extras:
        leg_y = H - 8
        for i, h in enumerate(h for h in HOMOLOGS if h in shown):
            lx = ML + i * 62
            svg += f'<rect x="{lx}" y="{leg_y - 8}" width="10" height="10" fill="{HOMOLOG_COLORS[h]}" opacity="0.65" rx="1"/>\n'
            svg += txt(lx + 14, leg_y + 1, h, size=9, color=INK2)

    svg += "</svg>"
    return "57-method-by-congener.svg", svg


# ---------------------------------------------------------------------------------------------
# 58  Per sample: peaks the analyst reported against peaks the pipeline detected
# ---------------------------------------------------------------------------------------------
GAP_BANDS = [(30, BAD, 0.7, "more than 30"), (10, MARGINAL, 0.6, "11 to 30"), (-math.inf, GOOD, 0.35, "10 or fewer")]


def fig_58(ds: Dataset) -> tuple[str, str]:
    sample_man = Counter()
    sample_det = Counter()
    for r in ds.sample_natives:
        key = (r["run"], r["sample"])
        if r["_m"] > 0:
            sample_man[key] += 1
        if r["_det"]:
            sample_det[key] += 1
    ax_max = max(Counter((r["run"], r["sample"]) for r in ds.sample_natives).values(), default=1)

    W, H = 500, 520
    ML, MR, MT, MB = 56, 20, 50, 60
    PLOT = min(W - ML - MR, H - MT - MB)

    svg = header(W, H)
    svg += txt(ML, 24, "Per-sample detection count", size=14, weight="700")
    svg += txt(ML, 40, "manual (mass > 0) vs automated (detected = True)", size=11, color=INK3)
    svg += f'<rect x="{ML}" y="{MT}" width="{PLOT}" height="{PLOT}" fill="white" rx="3"/>\n'

    def spx(v):
        return ML + v / ax_max * PLOT

    def spy(v):
        return MT + PLOT - v / ax_max * PLOT

    for g in [0, 50, 100, 150, ax_max]:
        svg += f'<line x1="{spx(g)}" y1="{spy(0)}" x2="{spx(g)}" y2="{spy(ax_max)}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += f'<line x1="{spx(0)}" y1="{spy(g)}" x2="{spx(ax_max)}" y2="{spy(g)}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(spx(g), MT + PLOT + 14, str(g), size=9, color=INK3, anchor="middle", fam=MONO)
        svg += txt(ML - 6, spy(g) + 4, str(g), size=9, color=INK3, anchor="end", fam=MONO)
    svg += (f'<line x1="{spx(0)}" y1="{spy(0)}" x2="{spx(ax_max)}" y2="{spy(ax_max)}" stroke="{INK3}"'
            f' stroke-width="1" stroke-dasharray="4,3" opacity="0.4"/>\n')

    for key in sorted(set(sample_man) | set(sample_det)):
        mx = sample_man.get(key, 0)
        dx = sample_det.get(key, 0)
        color, op = next((col, op) for limit, col, op, _ in GAP_BANDS if mx - dx > limit)
        svg += f'<circle cx="{spx(mx)}" cy="{spy(dx)}" r="3" fill="{color}" opacity="{op}"/>\n'

    svg += txt(ML + PLOT/2, MT + PLOT + 32, "congeners with manual mass > 0", size=11, color=INK2, anchor="middle", weight="500")
    svg += txt(ML - 40, MT + PLOT/2, "congeners with detected = True", size=11, color=INK2, anchor="middle",
               weight="500", rotate=-90)

    if ds.extras:
        leg_y = MT + PLOT + 50
        svg += txt(ML, leg_y + 3, "reported but not detected:", size=9, color=INK3)
        for i, (_, color, op, label) in enumerate(GAP_BANDS):
            lx = ML + 136 + i * 92
            svg += f'<circle cx="{lx}" cy="{leg_y}" r="3" fill="{color}" opacity="{op}"/>\n'
            svg += txt(lx + 8, leg_y + 3, label, size=9, color=INK2)

    svg += "</svg>"
    return "58-detection-agreement.svg", svg


FIGURES = [fig_51, fig_52, fig_53, fig_54, fig_55, fig_56, fig_57, fig_58]
