"""Figures 59 to 76: the congeners behind the mass disagreement, and behind the fingerprint one.

Part one ranks congeners by how much they contribute to total-mass error in the bad samples
(here, total-mass RPD above 25 %) and asks what those driver congeners have in common. Part
two does the same for the samples whose best-matching Aroclor differs between the manual and
the automated masses, called bad cos theta samples.
"""

from __future__ import annotations

import math
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass

from .data import Dataset
from .fingerprint import matches
from .svg import (HOMOLOG_COLORS, HOMOLOGS, INK2, INK3, METHOD_COLORS, METHOD_LABELS, MONO, PURPLE, RULE,
                  WARN, header, txt)

METH_ORDER = ["valley_drop", "rt_window", "deconvolved", "tangent_skim"]
HOM_LABELS = {h: h.capitalize() for h in HOMOLOGS}
GOOD_COLOR = "#0E6B68"
BAD_COLOR = "#C43D3D"
BAD_THRESH = 25
#: Aroclors from lightest to heaviest, for asking which way a group of samples leans.
AROCLOR_WEIGHT = ["1221", "1232", "1016", "1242", "1248", "1254", "1260", "1262", "1268"]


def rpd(m: float, a: float) -> float:
    avg = (m + a) / 2
    if avg == 0:
        return 0
    return abs(m - a) / avg * 100


@dataclass
class Context:
    """Everything the driver figures share, worked out once per dataset."""

    rows: list[dict]
    sample_total: dict
    bad_keys: set
    good_keys: set
    cong_homolog: dict
    top_drivers: list
    costheta: dict
    bad_ct_keys: set
    good_ct_keys: set

    @property
    def driver_labels(self) -> set:
        return {label for label, _ in self.top_drivers}

    @property
    def n_mismatch(self) -> int:
        return len(self.bad_ct_keys)

    def driver_homologs(self) -> list[str]:
        present = {self.cong_homolog[label] for label, _ in self.top_drivers}
        return [h for h in HOMOLOGS if h in present]


def _context(ds: Dataset) -> Context:
    rows = ds.sample_natives
    by_sample = defaultdict(list)
    for r in rows:
        by_sample[(r["run"], r["sample"])].append(r)
    sample_total = {}
    for key, rs in by_sample.items():
        tm = sum(r["_m"] for r in rs)
        ta = sum(r["_a"] for r in rs)
        sample_total[key] = {"manual": tm, "auto": ta, "rpd": rpd(tm, ta), "rows": rs}
    bad_keys = {k for k, v in sample_total.items() if v["rpd"] > BAD_THRESH}
    good_keys = {k for k, v in sample_total.items() if v["rpd"] <= BAD_THRESH}

    cong_contribution = defaultdict(float)
    for key in sorted(bad_keys):
        for r in sample_total[key]["rows"]:
            cong_contribution[r["label"]] += abs(r["_m"] - r["_a"])
    cong_homolog = {r["label"]: r["homolog"] for r in rows}
    top_drivers = sorted(cong_contribution.items(), key=lambda x: (-x[1], x[0]))[:20]

    # A sample with no mass on one side has no fingerprint, so no best match to compare. The
    # figures as first drawn counted such a sample as a mismatch; that is kept only for
    # redrawing them as they were.
    costheta = {}
    if ds.aroclors:
        for key, m in matches(ds).items():
            if not m.defined and ds.extras:
                continue
            costheta[key] = {
                "best_manual": m.best_manual, "best_auto": m.best_auto,
                "match_agrees": m.agrees,
                "ct_rpd": rpd(m.cos_manual, m.cos_auto) if m.cos_manual > 0 and m.cos_auto > 0 else 0,
            }
    bad_ct_keys = {k for k, v in costheta.items() if not v["match_agrees"]}
    return Context(rows=rows, sample_total=sample_total, bad_keys=bad_keys, good_keys=good_keys,
                   cong_homolog=cong_homolog, top_drivers=top_drivers, costheta=costheta,
                   bad_ct_keys=bad_ct_keys, good_ct_keys=set(costheta) - bad_ct_keys)


def context(ds: Dataset) -> Context:
    return ds.memo("driver_context", lambda: _context(ds))


def _nothing(name: str, title: str, message: str) -> tuple[str, str]:
    """A placeholder for a figure whose subject does not occur in this comparison."""
    svg = header(520, 120)
    svg += txt(24, 28, title, size=14, weight="700")
    svg += txt(24, 64, message, size=11, color=INK3)
    svg += "</svg>"
    return name, svg


def _label_margin(ds: Dataset) -> int:
    """Room left of a bar for a congener label. The 80 px first used cut off the longest names."""
    return 100 if ds.extras else 80


def _homolog_legend(homologs, ml, leg_y):
    out = ""
    for hi, h in enumerate(homologs):
        lx = ml + hi * 72
        out += f'<rect x="{lx}" y="{leg_y}" width="10" height="10" rx="2" fill="{HOMOLOG_COLORS[h]}"/>\n'
        out += txt(lx + 14, leg_y + 9, HOM_LABELS[h], size=9, color=INK3)
    return out


def _x_axis(ml, mr, mt, ph, w, label):
    return (f'<line x1="{ml}" y1="{mt + ph}" x2="{w - mr}" y2="{mt + ph}" stroke="{RULE}" stroke-width="1"/>\n'
            + txt((ml + w - mr) / 2, mt + ph + 16, label, size=10, anchor="middle", color=INK3))


def _y_title(ml, mt, ph, label, offset=8):
    return txt(ml - offset, (mt + mt + ph) / 2, label, size=10, color=INK3, anchor="middle", rotate=-90)


def _pair_legend(x1, x2, leg_y, bad_label, good_label, size=10, opacity="0.7"):
    return (f'<rect x="{x1}" y="{leg_y}" width="10" height="10" rx="2" fill="{BAD_COLOR}" opacity="{opacity}"/>\n'
            + txt(x1 + 14, leg_y + 9, bad_label, size=size, color=INK2)
            + f'<rect x="{x2}" y="{leg_y}" width="10" height="10" rx="2" fill="{GOOD_COLOR}" opacity="{opacity}"/>\n'
            + txt(x2 + 14, leg_y + 9, good_label, size=size, color=INK2))


NO_BAD_SAMPLES = f"No sample's total mass differs by more than {BAD_THRESH}%."
NO_BAD_COS = "Every sample's best-matching Aroclor is the same from the manual and the automated masses."


# ---------------------------------------------------------------------------------------------
# 59  The congeners that carry the mass error
# ---------------------------------------------------------------------------------------------
def fig_59(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    name, title = "59-congener-drivers.svg", "Top 20 congeners driving total mass error"
    if not ctx.top_drivers:
        return _nothing(name, title, NO_BAD_SAMPLES)

    W, H = 600, 560
    ML, MR, MT, MB = _label_margin(ds), 80, 56, 40
    PH = H - MT - MB
    n_bars = 20
    bar_h = PH / n_bars - 4
    max_contrib = ctx.top_drivers[0][1]

    svg = header(W, H)
    svg += txt(ML, 24, title, size=14, weight="700")
    svg += txt(ML, 40, f"cumulative |manual − auto| across {len(ctx.bad_keys)} samples with RPD > {BAD_THRESH}%",
               size=11, color=INK3)

    for i, (label, contrib) in enumerate(ctx.top_drivers):
        hom = ctx.cong_homolog[label]
        by = MT + i * (PH / n_bars) + 2
        bw = contrib / max_contrib * (W - ML - MR)
        svg += f'<rect x="{ML}" y="{by}" width="{bw}" height="{bar_h}" fill="{HOMOLOG_COLORS[hom]}" opacity="0.8" rx="2"/>\n'
        svg += txt(ML - 4, by + bar_h/2 + 4, f"PCB-{label}", size=9, anchor="end", fam=MONO, weight="500")
        svg += txt(ML + bw + 4, by + bar_h/2 + 4, f"{contrib:.0f} ng", size=9, color=INK2, fam=MONO)

    svg += _homolog_legend(ctx.driver_homologs(), ML, H - 20)
    svg += _x_axis(ML, MR, MT, PH, W, "cumulative |manual − auto| (ng)")
    svg += "</svg>"
    return name, svg


# ---------------------------------------------------------------------------------------------
# 60  Do the driver congeners disagree everywhere, or only in the bad samples?
# ---------------------------------------------------------------------------------------------
def fig_60(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    name, title = "60-driver-rpd-all-samples.svg", "Median |RPD| of driver congeners across all samples"
    cong_rpds = defaultdict(list)
    for r in ctx.rows:
        if r["_m"] > 0 and r["_a"] > 0 and r["label"] in ctx.driver_labels:
            cong_rpds[r["label"]].append(rpd(r["_m"], r["_a"]))
    driver_sorted = sorted(((label, statistics.median(vals)) for label, vals in cong_rpds.items()),
                           key=lambda x: (-x[1], x[0]))
    if not driver_sorted:
        return _nothing(name, title, NO_BAD_SAMPLES)

    extreme = [label for label, med in driver_sorted
               if med >= 50 and ctx.cong_homolog[label] in ("penta", "hexa")]
    subtitle = ("a few penta/hexa congeners carry extremely high individual RPDs" if extreme
                else "no driver congener's own median RPD reaches 50%")

    W, H = 600, 540
    ML, MR, MT, MB = _label_margin(ds), 60, 56, 40
    PH = H - MT - MB
    n_bars = len(driver_sorted)
    bar_h = min(PH / n_bars - 3, 20)
    max_rpd = max(v for _, v in driver_sorted) * 1.1

    svg = header(W, H)
    svg += txt(ML, 24, title, size=14, weight="700")
    svg += txt(ML, 40, subtitle, size=11, color=INK3)

    if max_rpd > 25:
        ref_x = ML + 25 / max_rpd * (W - ML - MR)
        svg += f'<line x1="{ref_x}" y1="{MT}" x2="{ref_x}" y2="{MT+PH}" stroke="{RULE}" stroke-width="1" stroke-dasharray="3,2"/>\n'
        svg += txt(ref_x, MT - 4, "25%", size=8, color=INK3, anchor="middle", fam=MONO)

    for i, (label, med_rpd) in enumerate(driver_sorted):
        hom = ctx.cong_homolog[label]
        by = MT + i * (PH / n_bars) + 1
        bw = med_rpd / max_rpd * (W - ML - MR)
        svg += f'<rect x="{ML}" y="{by}" width="{bw}" height="{bar_h}" fill="{HOMOLOG_COLORS[hom]}" opacity="0.8" rx="2"/>\n'
        svg += txt(ML - 4, by + bar_h/2 + 4, f"PCB-{label}", size=9, anchor="end", fam=MONO, weight="500")
        svg += txt(ML + bw + 4, by + bar_h/2 + 4, f"{med_rpd:.0f}%", size=9, color=INK2, fam=MONO)

    svg += _homolog_legend(ctx.driver_homologs(), ML, H - 20)
    svg += _x_axis(ML, MR, MT, PH, W, "median RPD (%)")
    svg += "</svg>"
    return name, svg


# ---------------------------------------------------------------------------------------------
# 61 and 66  Boundary methods in two groups of peaks, side by side
# ---------------------------------------------------------------------------------------------
def _method_pair(title, subtitle, groups):
    W, H = 520, 380
    ML, MR, MT, MB = 40, 40, 56, 70
    PW = W - ML - MR
    PH = H - MT - MB
    bar_w = PW / 2 - 40

    svg = header(W, H)
    svg += txt(W/2, 24, title, size=14, weight="700", anchor="middle")
    svg += txt(W/2, 40, subtitle, size=11, color=INK3, anchor="middle")

    for gi, (gname, mc, tot) in enumerate(groups):
        bx = ML + gi * (PW / 2) + 20
        svg += txt(bx + bar_w/2, MT - 4, gname, size=11, weight="600", anchor="middle")
        svg += txt(bx + bar_w/2, MT + 6, f"n = {tot:,}", size=9, color=INK3, anchor="middle", fam=MONO)
        cy = MT + PH + 14
        for m in METH_ORDER:
            frac = mc.get(m, 0) / tot if tot else 0
            seg_h = frac * PH
            cy -= seg_h
            svg += f'<rect x="{bx}" y="{cy}" width="{bar_w}" height="{seg_h}" fill="{METHOD_COLORS[m]}" opacity="0.8" rx="2"/>\n'
            if frac >= 0.05:
                svg += txt(bx + bar_w/2, cy + seg_h/2 + 4, f"{frac*100:.0f}%", size=10, color="white",
                           anchor="middle", weight="600", fam=MONO)

    for pct in [0, 25, 50, 75, 100]:
        yy = MT + PH - pct/100 * PH + 14
        svg += f'<line x1="{ML}" y1="{yy}" x2="{ML+PW}" y2="{yy}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(ML - 4, yy + 3, f"{pct}%", size=8, color=INK3, anchor="end", fam=MONO)

    leg_y = H - 30
    for mi, m in enumerate(METH_ORDER):
        lx = ML + mi * 120
        svg += f'<rect x="{lx}" y="{leg_y}" width="10" height="10" rx="2" fill="{METHOD_COLORS[m]}"/>\n'
        svg += txt(lx + 14, leg_y + 9, METHOD_LABELS[m], size=9, color=INK3)

    svg += _y_title(ML, MT, PH, "% of measurements", offset=28)
    svg += "</svg>"
    return svg


def _penta_hexa(ctx: Context):
    rows = [r for r in ctx.rows if r["homolog"] in ("penta", "hexa")]
    drivers = [r for r in rows if r["label"] in ctx.driver_labels]
    others = [r for r in rows if r["label"] not in ctx.driver_labels]
    return drivers, others


def fig_61(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    name, title = "61-driver-method-comparison.svg", "Integration method: driver congeners vs others"
    if not ctx.top_drivers:
        return _nothing(name, title, NO_BAD_SAMPLES)
    drivers, others = _penta_hexa(ctx)
    return name, _method_pair(title, "within penta and hexa homologs only", [
        ("Driver congeners", Counter(r["method"] for r in drivers), len(drivers)),
        ("Other penta/hexa", Counter(r["method"] for r in others), len(others)),
    ])


# ---------------------------------------------------------------------------------------------
# 62  Are the driver congeners simply the big ones?
# ---------------------------------------------------------------------------------------------
def fig_62(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    name, title = "62-driver-detection-mass.svg", "Detection rate and median mass: drivers vs others"
    if not ctx.top_drivers:
        return _nothing(name, title, NO_BAD_SAMPLES)
    drivers, others = _penta_hexa(ctx)
    driver_rate = sum(1 for r in drivers if r["_det"]) / len(drivers) * 100 if drivers else 0
    other_rate = sum(1 for r in others if r["_det"]) / len(others) * 100 if others else 0
    driver_masses = [r["_m"] for r in drivers if r["_m"] > 0]
    other_masses = [r["_m"] for r in others if r["_m"] > 0]
    driver_mass = statistics.median(driver_masses) if driver_masses else 0
    other_mass = statistics.median(other_masses) if other_masses else 0

    size = "larger" if driver_mass > other_mass else "no larger"
    if driver_rate - other_rate >= 1:
        found = "and better detected"
    elif other_rate - driver_rate >= 1:
        found = "but detected less often"
    else:
        found = "and detected as often"
    subtitle = f"within penta/hexa; driver congeners are {size} {found}"

    W, H = 520, 380
    ML, MR, MT, MB = 40, 40, 56, 60
    PW = W - ML - MR
    PH = H - MT - MB

    svg = header(W, H)
    svg += txt(W/2, 24, title, size=14, weight="700", anchor="middle")
    svg += txt(W/2, 40, subtitle, size=11, color=INK3, anchor="middle")

    bar_w = 60
    for mi, (mname, dv, ov, unit) in enumerate([
        ("Detection rate", driver_rate, other_rate, "%"),
        ("Median mass (ng)", driver_mass, other_mass, " ng"),
    ]):
        group_x = ML + 30 + mi * (PW / 2)
        max_val = max(dv, ov) * 1.15
        if max_val == 0:
            max_val = 1
        svg += txt(group_x + bar_w + 5, MT + 10, mname, size=11, weight="600", anchor="middle")

        bh_d = dv / max_val * PH * 0.85
        by_d = MT + PH - bh_d + 14
        svg += f'<rect x="{group_x}" y="{by_d}" width="{bar_w}" height="{bh_d}" fill="{BAD_COLOR}" opacity="0.7" rx="2"/>\n'
        svg += txt(group_x + bar_w/2, by_d - 4, f"{dv:.1f}{unit}", size=10, anchor="middle", fam=MONO, weight="500")

        bh_o = ov / max_val * PH * 0.85
        by_o = MT + PH - bh_o + 14
        svg += f'<rect x="{group_x + bar_w + 10}" y="{by_o}" width="{bar_w}" height="{bh_o}" fill="{GOOD_COLOR}" opacity="0.7" rx="2"/>\n'
        svg += txt(group_x + bar_w + 10 + bar_w/2, by_o - 4, f"{ov:.1f}{unit}", size=10, anchor="middle", fam=MONO, weight="500")

    svg += _pair_legend(ML + 80, ML + 230, H - 26, "Driver congeners", "Other penta/hexa")
    svg += "</svg>"
    return name, svg


# ---------------------------------------------------------------------------------------------
# 63  Driver congeners in the bad samples and in the good ones
# ---------------------------------------------------------------------------------------------
def fig_63(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    name, title = "63-driver-rpd-bad-vs-good.svg", "Driver congener RPD: bad samples vs good samples"
    if not ctx.top_drivers:
        return _nothing(name, title, NO_BAD_SAMPLES)

    cong_rpd_bad = defaultdict(list)
    cong_rpd_good = defaultdict(list)
    for r in ctx.rows:
        if r["_m"] > 0 and r["_a"] > 0 and r["label"] in ctx.driver_labels:
            key = (r["run"], r["sample"])
            if key in ctx.bad_keys:
                cong_rpd_bad[r["label"]].append(rpd(r["_m"], r["_a"]))
            elif key in ctx.good_keys:
                cong_rpd_good[r["label"]].append(rpd(r["_m"], r["_a"]))

    driver_bad_good = []
    for label in sorted(ctx.driver_labels):
        med_bad = statistics.median(cong_rpd_bad[label]) if cong_rpd_bad[label] else 0
        med_good = statistics.median(cong_rpd_good[label]) if cong_rpd_good[label] else 0
        driver_bad_good.append((label, med_bad, med_good))
    driver_bad_good.sort(key=lambda x: -x[1])

    W, H = 600, 560
    ML, MR, MT, MB = _label_margin(ds), 60, 56, 40
    PH = H - MT - MB
    PW = W - ML - MR
    n_bars = len(driver_bad_good)
    bar_h = min(PH / n_bars - 4, 18)
    max_rpd_bg = max(max(b, g) for _, b, g in driver_bad_good) * 1.05 or 1

    svg = header(W, H)
    svg += txt(ML, 24, title, size=14, weight="700")
    svg += txt(ML, 40, "median |RPD| per congener, split by overall sample quality", size=11, color=INK3)

    for i, (label, med_bad, med_good) in enumerate(driver_bad_good):
        by = MT + i * (PH / n_bars) + 1
        bw_bad = med_bad / max_rpd_bg * PW
        svg += f'<rect x="{ML}" y="{by}" width="{bw_bad}" height="{bar_h/2}" fill="{BAD_COLOR}" opacity="0.8" rx="1"/>\n'
        bw_good = med_good / max_rpd_bg * PW
        svg += f'<rect x="{ML}" y="{by + bar_h/2}" width="{bw_good}" height="{bar_h/2}" fill="{GOOD_COLOR}" opacity="0.8" rx="1"/>\n'
        svg += txt(ML - 4, by + bar_h/2 + 3, f"PCB-{label}", size=9, anchor="end", fam=MONO, weight="500")
        svg += txt(ML + max(bw_bad, bw_good) + 4, by + bar_h/2 + 3, f"{med_bad:.0f}% / {med_good:.0f}%",
                   size=8, color=INK3, fam=MONO)

    svg += _pair_legend(ML, ML + 200, H - 22, f"Bad samples (RPD > {BAD_THRESH}%)",
                        f"Good samples (RPD ≤ {BAD_THRESH}%)", size=9, opacity="0.8")
    svg += _x_axis(ML, MR, MT, PH, W, "median RPD (%)")
    svg += "</svg>"
    return name, svg


# ---------------------------------------------------------------------------------------------
# 64  How far apart are the two best-match scores?
# ---------------------------------------------------------------------------------------------
def fig_64(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    CT_BINS = [(0, 1), (1, 2), (2, 3), (3, 5), (5, 10), (10, 20), (20, 50), (50, 200)]
    CT_BIN_LABELS = ["0–1%", "1–2%", "2–3%", "3–5%", "5–10%", "10–20%", "20–50%", ">50%"]
    bin_match = [sum(1 for v in ctx.costheta.values() if lo <= v["ct_rpd"] < hi and v["match_agrees"]) for lo, hi in CT_BINS]
    bin_mismatch = [sum(1 for v in ctx.costheta.values() if lo <= v["ct_rpd"] < hi and not v["match_agrees"]) for lo, hi in CT_BINS]

    W, H = 520, 400
    ML, MR, MT, MB = 56, 30, 56, 74
    PW = W - ML - MR
    PH = H - MT - MB
    bar_w = PW / len(CT_BINS) - 6
    max_count = max(m + mm for m, mm in zip(bin_match, bin_mismatch))

    svg = header(W, H)
    svg += txt(ML, 24, "Best-match cos θ RPD distribution", size=14, weight="700")
    svg += txt(ML, 40, f"{ctx.n_mismatch} of {len(ctx.costheta)} samples have a different best aroclor match",
               size=11, color=INK3)

    for bi, label in enumerate(CT_BIN_LABELS):
        bx = ML + bi * (PW / len(CT_BINS)) + 3
        nm, nmm = bin_match[bi], bin_mismatch[bi]
        bh_m = nm / max_count * PH if max_count > 0 else 0
        bh_mm = nmm / max_count * PH if max_count > 0 else 0
        by_m = MT + PH - bh_m
        by_mm = by_m - bh_mm
        svg += f'<rect x="{bx}" y="{by_m}" width="{bar_w}" height="{bh_m}" fill="{GOOD_COLOR}" opacity="0.7" rx="2"/>\n'
        if bh_mm > 0:
            svg += f'<rect x="{bx}" y="{by_mm}" width="{bar_w}" height="{bh_mm}" fill="{BAD_COLOR}" opacity="0.8" rx="2"/>\n'
        svg += txt(bx + bar_w/2, by_mm - 4, str(nm + nmm), size=9, anchor="middle", fam=MONO, weight="500")
        svg += txt(bx + bar_w/2, MT + PH + 14, label, size=9, anchor="middle", color=INK2, fam=MONO)

    svg += txt(W/2, MT + PH + 32, "RPD of best-match cos θ score", size=10, anchor="middle", color=INK3)
    svg += f'<line x1="{ML}" y1="{MT+PH}" x2="{ML+PW}" y2="{MT+PH}" stroke="{RULE}" stroke-width="1"/>\n'

    leg_y = H - 22
    svg += f'<rect x="{ML + 60}" y="{leg_y}" width="10" height="10" rx="2" fill="{GOOD_COLOR}" opacity="0.7"/>\n'
    svg += txt(ML + 74, leg_y + 9, "Same best aroclor", size=9, color=INK2)
    svg += f'<rect x="{ML + 230}" y="{leg_y}" width="10" height="10" rx="2" fill="{BAD_COLOR}" opacity="0.8"/>\n'
    svg += txt(ML + 244, leg_y + 9, "Different best aroclor", size=9, color=INK2)

    svg += _y_title(ML, MT, PH, "samples")
    svg += "</svg>"
    return "64-costheta-rpd-distribution.svg", svg


# ---------------------------------------------------------------------------------------------
# 65 and 73  A homolog-by-homolog comparison of the bad and the good cos theta samples
# ---------------------------------------------------------------------------------------------
def _homolog_pairs(ctx: Context, title, subtitle, bad_of, good_of, label_from, y_label):
    W, H = 560, 400
    ML, MR, MT, MB = 56, 30, 56, 64
    PW = W - ML - MR
    PH = H - MT - MB
    bar_w = PW / len(HOMOLOGS) - 6
    half_w = (bar_w - 4) / 2
    max_val = max([bad_of[h] for h in HOMOLOGS] + [good_of[h] for h in HOMOLOGS]) * 1.15

    svg = header(W, H)
    svg += txt(ML, 24, title, size=14, weight="700")
    svg += txt(ML, 40, subtitle, size=11, color=INK3)

    for hi, h in enumerate(HOMOLOGS):
        bx = ML + hi * (PW / len(HOMOLOGS)) + 3
        mb, mg = bad_of[h], good_of[h]
        bh_b = mb / max_val * PH if max_val > 0 else 0
        svg += f'<rect x="{bx}" y="{MT + PH - bh_b}" width="{half_w}" height="{bh_b}" fill="{BAD_COLOR}" opacity="0.7" rx="2"/>\n'
        if mb >= label_from:
            svg += txt(bx + half_w/2, MT + PH - bh_b - 3, f"{mb:.0f}%", size=8, anchor="middle", fam=MONO, color=BAD_COLOR)
        bh_g = mg / max_val * PH if max_val > 0 else 0
        svg += f'<rect x="{bx + half_w + 4}" y="{MT + PH - bh_g}" width="{half_w}" height="{bh_g}" fill="{GOOD_COLOR}" opacity="0.7" rx="2"/>\n'
        if mg >= label_from:
            svg += txt(bx + half_w + 4 + half_w/2, MT + PH - bh_g - 3, f"{mg:.0f}%", size=8, anchor="middle", fam=MONO, color=GOOD_COLOR)
        svg += txt(bx + bar_w/2, MT + PH + 14, h[:4].capitalize(), size=9, anchor="middle", color=INK2, fam=MONO)

    svg += f'<line x1="{ML}" y1="{MT+PH}" x2="{ML+PW}" y2="{MT+PH}" stroke="{RULE}" stroke-width="1"/>\n'
    svg += _pair_legend(ML + 60, ML + 260, H - 26, f"Bad cos θ (n={len(ctx.bad_ct_keys)})",
                        f"Good cos θ (n={len(ctx.good_ct_keys)})")
    svg += _y_title(ML, MT, PH, y_label)
    svg += "</svg>"
    return svg


def fig_65(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    name, title = "65-costheta-homolog-profile.svg", "Homolog composition: bad vs good cos θ samples"
    if not ctx.bad_ct_keys:
        return _nothing(name, title, NO_BAD_COS)

    fracs_bad = defaultdict(list)
    fracs_good = defaultdict(list)
    for key, info in ctx.sample_total.items():
        if info["manual"] <= 0:
            continue
        hom_mass = defaultdict(float)
        for r in info["rows"]:
            hom_mass[r["homolog"]] += r["_m"]
        # As first drawn, every sample outside the bad group counted as good here.
        group = fracs_bad if key in ctx.bad_ct_keys else fracs_good if (key in ctx.good_ct_keys or not ds.extras) else None
        if group is None:
            continue
        for h in HOMOLOGS:
            group[h].append(hom_mass[h] / info["manual"] * 100)

    bad_of = {h: statistics.median(fracs_bad[h]) if fracs_bad[h] else 0 for h in HOMOLOGS}
    good_of = {h: statistics.median(fracs_good[h]) if fracs_good[h] else 0 for h in HOMOLOGS}
    return name, _homolog_pairs(ctx, title, "median mass fraction by homolog group (manual values)",
                                bad_of, good_of, 1, "median mass fraction (%)")


def fig_66(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    name, title = "66-costheta-method-comparison.svg", "Integration method: bad vs good cos θ samples"
    if not ctx.bad_ct_keys:
        return _nothing(name, title, NO_BAD_COS)
    meth_bad = Counter()
    meth_good = Counter()
    for r in ctx.rows:
        key = (r["run"], r["sample"])
        if key in ctx.bad_ct_keys:
            meth_bad[r["method"]] += 1
        elif key in ctx.good_ct_keys:
            meth_good[r["method"]] += 1
    return name, _method_pair(title, "samples where auto picks a different best aroclor", [
        ("Bad cos θ", meth_bad, sum(meth_bad.values())),
        ("Good cos θ", meth_good, sum(meth_good.values())),
    ])


# ---------------------------------------------------------------------------------------------
# 67 and 68  A per-sample quantity, as a share of the bad and of the good cos theta samples
# ---------------------------------------------------------------------------------------------
def _group_histogram(w, h, mb, title, subtitle, bins, labels, bad_values, good_values, x_label, ticks=False,
                     x_label_dy=30, title_offset=8):
    ML, MR, MT = 56, 30, 56
    PW = w - ML - MR
    PH = h - MT - mb
    bad_pct = [sum(1 for x in bad_values if lo <= x < hi) / len(bad_values) * 100 if bad_values else 0 for lo, hi in bins]
    good_pct = [sum(1 for x in good_values if lo <= x < hi) / len(good_values) * 100 if good_values else 0 for lo, hi in bins]
    max_pct = max(max(bad_pct), max(good_pct)) * 1.15
    bar_w = PW / len(bins) - 4
    half_bar = (bar_w - 2) / 2

    svg = header(w, h)
    svg += txt(ML, 24, title, size=14, weight="700")
    svg += txt(ML, 40, subtitle, size=11, color=INK3)

    for bi, label in enumerate(labels):
        bx = ML + bi * (PW / len(bins)) + 2
        bh_b = bad_pct[bi] / max_pct * PH if max_pct > 0 else 0
        svg += f'<rect x="{bx}" y="{MT + PH - bh_b}" width="{half_bar}" height="{bh_b}" fill="{BAD_COLOR}" opacity="0.7" rx="1"/>\n'
        bh_g = good_pct[bi] / max_pct * PH if max_pct > 0 else 0
        svg += f'<rect x="{bx + half_bar + 2}" y="{MT + PH - bh_g}" width="{half_bar}" height="{bh_g}" fill="{GOOD_COLOR}" opacity="0.7" rx="1"/>\n'
        svg += txt(bx + bar_w/2, MT + PH + 14, label, size=8, anchor="middle", color=INK2, fam=MONO)

    svg += f'<line x1="{ML}" y1="{MT+PH}" x2="{ML+PW}" y2="{MT+PH}" stroke="{RULE}" stroke-width="1"/>\n'
    svg += txt(w/2, MT + PH + x_label_dy, x_label, size=10, anchor="middle", color=INK3)

    legend = (f'<rect x="{ML + PW - 160}" y="{MT + 10}" width="10" height="10" rx="2" fill="{BAD_COLOR}" opacity="0.7"/>\n'
              + txt(ML + PW - 146, MT + 19, "Bad cos θ", size=9, color=INK2)
              + f'<rect x="{ML + PW - 160}" y="{MT + 26}" width="10" height="10" rx="2" fill="{GOOD_COLOR}" opacity="0.7"/>\n'
              + txt(ML + PW - 146, MT + 35, "Good cos θ", size=9, color=INK2))
    if ticks:
        svg += _y_title(ML, MT, PH, "% of group", offset=title_offset)
        for pct_tick in range(10, int(max_pct) + 1, 10):
            ty = MT + PH - (pct_tick / max_pct * PH)
            svg += f'<line x1="{ML - 3}" y1="{ty}" x2="{ML}" y2="{ty}" stroke="{RULE}" stroke-width="1"/>\n'
            svg += txt(ML - 5, ty + 3, f"{pct_tick}", size=8, anchor="end", color=INK3, fam=MONO)
        svg += legend
    else:
        svg += legend
        svg += _y_title(ML, MT, PH, "% of group")
    svg += "</svg>"
    return svg


def fig_67(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    name, title = "67-costheta-miss-rate.svg", "Missed peak rate: bad vs good cos θ"
    if not ctx.bad_ct_keys:
        return _nothing(name, title, NO_BAD_COS)

    miss_bad, miss_good = [], []
    for key, info in ctx.sample_total.items():
        total_peaks = len(info["rows"])
        miss_pct = sum(1 for r in info["rows"] if not r["_det"]) / total_peaks * 100 if total_peaks else 0
        if key in ctx.bad_ct_keys:
            miss_bad.append(miss_pct)
        elif key in ctx.good_ct_keys:
            miss_good.append(miss_pct)

    bins = [(0, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60), (60, 70), (70, 80), (80, 100)]
    labels = ["0–10", "10–20", "20–30", "30–40", "40–50", "50–60", "60–70", "70–80", ">80"]
    subtitle = (f"bad cos θ: median {statistics.median(miss_bad):.0f}%  |  "
                f"good cos θ: median {statistics.median(miss_good):.0f}%")
    return name, _group_histogram(480, 380, 50, title, subtitle, bins, labels, miss_bad, miss_good,
                                  "missed peak rate (%)")


def fig_68(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    name, title = "68-costheta-total-mass.svg", "Total sample mass: bad vs good cos θ"
    if not ctx.bad_ct_keys:
        return _nothing(name, title, NO_BAD_COS)

    bad_totals = [ctx.sample_total[k]["manual"] for k in sorted(ctx.bad_ct_keys)]
    good_totals = [ctx.sample_total[k]["manual"] for k in sorted(ctx.good_ct_keys)]
    bins = [(0, 50), (50, 100), (100, 250), (250, 500), (500, 1000), (1000, 5000), (5000, 50000)]
    labels = ["0–50", "50–100", "100–250", "250–500", "500–1k", "1k–5k", ">5k"]
    subtitle = (f"bad median = {statistics.median(bad_totals):.0f} ng  |  "
                f"good median = {statistics.median(good_totals):.0f} ng")
    return name, _group_histogram(520, 380, 60, title, subtitle, bins, labels, bad_totals, good_totals,
                                  "total sample mass (ng)", ticks=True, x_label_dy=32,
                                  title_offset=30 if ds.extras else 8)


# ---------------------------------------------------------------------------------------------
# 69  The congeners whose share of the profile moves most
# ---------------------------------------------------------------------------------------------
def _profile(sample_rows, field):
    total = sum(r[field] for r in sample_rows)
    if total <= 0:
        return {}
    return {r["label"]: r[field] / total for r in sample_rows}


def cos_drivers(ctx: Context) -> list[tuple[str, float]]:
    """The 20 congeners whose mass fraction differs most, on average, in the bad cos theta samples."""
    cong_frac_diff = defaultdict(list)
    for key in sorted(ctx.bad_ct_keys):
        rows = ctx.sample_total[key]["rows"]
        prof_m = _profile(rows, "_m")
        prof_a = _profile(rows, "_a")
        for lbl in prof_m:
            cong_frac_diff[lbl].append(abs(prof_m.get(lbl, 0) - prof_a.get(lbl, 0)))
    mean_diff = {lbl: statistics.mean(vals) for lbl, vals in cong_frac_diff.items() if vals}
    return sorted(mean_diff.items(), key=lambda x: (-x[1], x[0]))[:20]


def fig_69(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    name, title = "69-costheta-congener-drivers.svg", "Top congeners driving cos θ disagreement"
    top = cos_drivers(ctx)
    if not top:
        return _nothing(name, title, NO_BAD_COS)

    W, H = 600, 560
    ML, MR, MT, MB = _label_margin(ds), 80, 56, 40
    PH = H - MT - MB
    n_bars = len(top)
    bar_h = min(PH / n_bars - 3, 20)
    max_diff = top[0][1] * 1.05 or 1

    svg = header(W, H)
    svg += txt(ML, 24, title, size=14, weight="700")
    svg += txt(ML, 40, f"mean |manual fraction − auto fraction| across {len(ctx.bad_ct_keys)} bad cos θ samples",
               size=11, color=INK3)

    homs = set()
    for i, (label, diff) in enumerate(top):
        hom = ctx.cong_homolog.get(label, "?")
        homs.add(hom)
        by = MT + i * (PH / n_bars) + 1
        bw = diff / max_diff * (W - ML - MR)
        svg += f'<rect x="{ML}" y="{by}" width="{bw}" height="{bar_h}" fill="{HOMOLOG_COLORS.get(hom, INK3)}" opacity="0.8" rx="2"/>\n'
        svg += txt(ML - 4, by + bar_h/2 + 4, f"PCB-{label}", size=9, anchor="end", fam=MONO, weight="500")
        svg += txt(ML + bw + 4, by + bar_h/2 + 4, f"{diff*100:.2f}%", size=9, color=INK2, fam=MONO)

    svg += _homolog_legend([h for h in HOMOLOGS if h in homs], ML, H - 20)
    svg += _x_axis(ML, MR, MT, PH, W, "mean |fraction difference| (%)")
    svg += "</svg>"
    return name, svg


# ---------------------------------------------------------------------------------------------
# 70  Are the bad-mass samples and the bad-fingerprint samples the same ones?
# ---------------------------------------------------------------------------------------------
def fig_70(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    n_both = len(ctx.bad_keys & ctx.bad_ct_keys)
    n_only_mass = len(ctx.bad_keys - ctx.bad_ct_keys)
    n_only_ct = len(ctx.bad_ct_keys - ctx.bad_keys)
    total = len(ctx.sample_total)
    n_neither = total - n_both - n_only_mass - n_only_ct

    W, H = 480, 340
    svg = header(W, H)
    svg += txt(W/2, 24, "Overlap: bad total mass vs bad cos θ", size=14, weight="700", anchor="middle")
    svg += txt(W/2, 40, "bad cos θ = best aroclor match disagrees", size=11, color=INK3, anchor="middle")

    cx_l, cx_r = W/2 - 56, W/2 + 56
    cy_e = H/2 + 6
    erx, ery = 120, 90
    svg += f'<ellipse cx="{cx_l}" cy="{cy_e}" rx="{erx}" ry="{ery}" fill="{BAD_COLOR}" opacity="0.25" stroke="{BAD_COLOR}" stroke-width="1.5"/>\n'
    svg += f'<ellipse cx="{cx_r}" cy="{cy_e}" rx="{erx}" ry="{ery}" fill="{WARN}" opacity="0.25" stroke="{WARN}" stroke-width="1.5"/>\n'

    svg += txt(cx_l - 40, cy_e - ery - 10, "Mass RPD > 25%", size=11, weight="600", color=BAD_COLOR)
    svg += txt(cx_r + 40, cy_e - ery - 10, "Bad cos θ", size=11, weight="600", color=WARN, anchor="end")
    svg += txt(cx_l - 40, cy_e - ery - 24, f"n = {n_only_mass + n_both}", size=10, fam=MONO, color=INK2)
    svg += txt(cx_r + 40, cy_e - ery - 24, f"n = {n_only_ct + n_both}", size=10, fam=MONO, color=INK2, anchor="end")

    svg += txt(cx_l - 44, cy_e + 5, str(n_only_mass), size=22, weight="700", anchor="middle", color=BAD_COLOR)
    svg += txt(W/2, cy_e + 5, str(n_both), size=22, weight="700", anchor="middle", color=PURPLE)
    svg += txt(cx_r + 44, cy_e + 5, str(n_only_ct), size=22, weight="700", anchor="middle", color=WARN)
    svg += txt(cx_l - 44, cy_e + 22, "only", size=9, anchor="middle", color=INK3)
    svg += txt(W/2, cy_e + 22, "both", size=9, anchor="middle", color=INK3)
    svg += txt(cx_r + 44, cy_e + 22, "only", size=9, anchor="middle", color=INK3)

    share = n_neither / total * 100 if total else 0
    svg += txt(W/2, H - 18, f"Neither: {n_neither} of {total} samples ({share:.0f}%)", size=11, fam=MONO,
               color=INK3, anchor="middle")
    svg += "</svg>"
    return "70-overlap-bad-mass-costheta.svg", svg


# ---------------------------------------------------------------------------------------------
# 71  Per sample: mass disagreement against fingerprint disagreement
# ---------------------------------------------------------------------------------------------
def fig_71(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    W, H = 500, 500
    ML, MR, MT, MB = 56, 30, 56, 50
    PLOT = min(W - ML - MR, H - MT - MB)

    svg = header(W, H)
    svg += txt(ML, 24, "Total mass RPD vs cos θ RPD", size=14, weight="700")
    svg += txt(ML, 40, "each dot = one sample; quadrants show co-occurrence", size=11, color=INK3)
    svg += f'<rect x="{ML}" y="{MT}" width="{PLOT}" height="{PLOT}" fill="white" rx="3"/>\n'

    max_mass_rpd = min(max(v["rpd"] for v in ctx.sample_total.values()), 200) or 1
    max_ct_rpd = min(max((v["ct_rpd"] for v in ctx.costheta.values()), default=0), 100) or 1

    for pct in range(0, int(max_mass_rpd) + 1, 25):
        fx = ML + pct / max_mass_rpd * PLOT
        svg += f'<line x1="{fx}" y1="{MT}" x2="{fx}" y2="{MT+PLOT}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(fx, MT + PLOT + 14, f"{pct}%", size=8, anchor="middle", color=INK3, fam=MONO)
    for pct in range(0, int(max_ct_rpd) + 1, 10):
        fy = MT + PLOT - pct / max_ct_rpd * PLOT
        svg += f'<line x1="{ML}" y1="{fy}" x2="{ML+PLOT}" y2="{fy}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(ML - 4, fy + 3, f"{pct}%", size=8, anchor="end", color=INK3, fam=MONO)

    thresh_x = ML + BAD_THRESH / max_mass_rpd * PLOT
    svg += (f'<line x1="{thresh_x}" y1="{MT}" x2="{thresh_x}" y2="{MT+PLOT}" stroke="{BAD_COLOR}" stroke-width="1"'
            f' stroke-dasharray="4,3" opacity="0.5"/>\n')
    ct_thresh_y = MT + PLOT - 5 / max_ct_rpd * PLOT
    svg += (f'<line x1="{ML}" y1="{ct_thresh_y}" x2="{ML+PLOT}" y2="{ct_thresh_y}" stroke="{BAD_COLOR}" stroke-width="1"'
            f' stroke-dasharray="4,3" opacity="0.5"/>\n')

    for key, info in ctx.sample_total.items():
        mass_rpd_val = min(info["rpd"], max_mass_rpd)
        ct_rpd_val = min(ctx.costheta[key]["ct_rpd"], max_ct_rpd) if key in ctx.costheta else 0
        px = ML + mass_rpd_val / max_mass_rpd * PLOT
        py = MT + PLOT - ct_rpd_val / max_ct_rpd * PLOT
        bad_mass, bad_ct = key in ctx.bad_keys, key in ctx.bad_ct_keys
        color = PURPLE if (bad_mass and bad_ct) else BAD_COLOR if bad_mass else WARN if bad_ct else GOOD_COLOR
        is_bad = bad_mass or bad_ct
        svg += f'<circle cx="{px}" cy="{py}" r="{"3" if is_bad else "2"}" fill="{color}" opacity="{"0.6" if is_bad else "0.3"}"/>\n'

    svg += txt(ML + PLOT/2, MT + PLOT + 32, "total mass |RPD|", size=10, anchor="middle", color=INK3)
    svg += txt(ML - 36, MT + PLOT/2, "cos θ |RPD|", size=10, color=INK3, anchor="middle", rotate=-90)

    for li, (lc, ll) in enumerate([(PURPLE, "Both bad"), (BAD_COLOR, "Only mass RPD bad"),
                                   (WARN, "Only cos θ bad"), (GOOD_COLOR, "Neither (good)")]):
        lx = ML + PLOT - 140
        ly = MT + 10 + li * 16
        svg += f'<circle cx="{lx}" cy="{ly}" r="3" fill="{lc}" opacity="0.6"/>\n'
        svg += txt(lx + 8, ly + 4, ll, size=9, color=INK2)

    svg += "</svg>"
    return "71-mass-rpd-vs-costheta-rpd.svg", svg


# ---------------------------------------------------------------------------------------------
# 72  Which Aroclors do the bad cos theta samples look like?
# ---------------------------------------------------------------------------------------------
def fig_72(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    name, title = "72-costheta-aroclor-distribution.svg", "Best aroclor match: bad vs good cos θ samples"
    if not ctx.bad_ct_keys:
        return _nothing(name, title, NO_BAD_COS)

    aroclor_bad = Counter(ctx.costheta[k]["best_manual"] for k in ctx.bad_ct_keys)
    aroclor_good = Counter(ctx.costheta[k]["best_manual"] for k in ctx.good_ct_keys)
    active = sorted(set(aroclor_bad) | set(aroclor_good), key=int)
    bad_total = len(ctx.bad_ct_keys)
    good_total = len(ctx.good_ct_keys)

    def lean(counts, total):
        ranks = [AROCLOR_WEIGHT.index(a) * n for a, n in counts.items() if a in AROCLOR_WEIGHT]
        return sum(ranks) / total if total else 0

    difference = lean(aroclor_bad, bad_total) - lean(aroclor_good, good_total)
    leaning = ("bad samples skew toward lighter aroclors" if difference <= -0.3
               else "bad samples skew toward heavier aroclors" if difference >= 0.3
               else "bad and good samples match the same aroclors")

    W, H = 560, 420
    ML, MR, MT, MB = 56, 30, 56, 60
    PW = W - ML - MR
    PH = H - MT - MB
    bar_w = PW / len(active) - 6
    half_bar = (bar_w - 4) / 2

    svg = header(W, H)
    svg += txt(ML, 24, title, size=14, weight="700")
    svg += txt(ML, 40, f"by manual identification; {leaning}", size=11, color=INK3)

    def pct_of(counts, total, a):
        return counts.get(a, 0) / total * 100 if total else 0

    max_pct = max(max(pct_of(aroclor_bad, bad_total, a), pct_of(aroclor_good, good_total, a)) for a in active) * 1.15

    for ai, a in enumerate(active):
        bx = ML + ai * (PW / len(active)) + 3
        bp = pct_of(aroclor_bad, bad_total, a)
        gp = pct_of(aroclor_good, good_total, a)
        bh_b = bp / max_pct * PH if max_pct > 0 else 0
        svg += f'<rect x="{bx}" y="{MT + PH - bh_b}" width="{half_bar}" height="{bh_b}" fill="{BAD_COLOR}" opacity="0.7" rx="2"/>\n'
        if bp >= 3:
            svg += txt(bx + half_bar/2, MT + PH - bh_b - 3, f"{bp:.0f}%", size=8, anchor="middle", fam=MONO, color=BAD_COLOR)
        bh_g = gp / max_pct * PH if max_pct > 0 else 0
        svg += f'<rect x="{bx + half_bar + 4}" y="{MT + PH - bh_g}" width="{half_bar}" height="{bh_g}" fill="{GOOD_COLOR}" opacity="0.7" rx="2"/>\n'
        if gp >= 3:
            svg += txt(bx + half_bar + 4 + half_bar/2, MT + PH - bh_g - 3, f"{gp:.0f}%", size=8, anchor="middle", fam=MONO, color=GOOD_COLOR)
        svg += txt(bx + bar_w/2, MT + PH + 14, a, size=9, anchor="middle", color=INK2, fam=MONO)

    svg += f'<line x1="{ML}" y1="{MT+PH}" x2="{ML+PW}" y2="{MT+PH}" stroke="{RULE}" stroke-width="1"/>\n'
    svg += txt(W/2, MT + PH + 32, "Aroclor", size=10, anchor="middle", color=INK3)
    svg += _pair_legend(ML + 80, ML + 260, H - 26, f"Bad cos θ (n={bad_total})", f"Good cos θ (n={good_total})")
    svg += _y_title(ML, MT, PH, "% of group")
    svg += "</svg>"
    return name, svg


def fig_73(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    name, title = "73-costheta-rpd-by-homolog.svg", "Median |RPD| by homolog: bad vs good cos θ"
    if not ctx.bad_ct_keys:
        return _nothing(name, title, NO_BAD_COS)

    rpd_bad = defaultdict(list)
    rpd_good = defaultdict(list)
    for r in ctx.rows:
        if r["_m"] > 0 and r["_a"] > 0:
            key = (r["run"], r["sample"])
            if key in ctx.bad_ct_keys:
                rpd_bad[r["homolog"]].append(rpd(r["_m"], r["_a"]))
            elif key in ctx.good_ct_keys:
                rpd_good[r["homolog"]].append(rpd(r["_m"], r["_a"]))
    bad_of = {h: statistics.median(rpd_bad[h]) if rpd_bad[h] else 0 for h in HOMOLOGS}
    good_of = {h: statistics.median(rpd_good[h]) if rpd_good[h] else 0 for h in HOMOLOGS}

    def mean_gap(homologs):
        return statistics.mean(bad_of[h] - good_of[h] for h in homologs)

    light, heavy = HOMOLOGS[:4], HOMOLOGS[4:]
    if mean_gap(heavy) - mean_gap(light) >= 5:
        subtitle = "heavier homologs show greater divergence in bad cos θ samples"
    elif mean_gap(light) - mean_gap(heavy) >= 5:
        subtitle = "lighter homologs show greater divergence in bad cos θ samples"
    else:
        subtitle = "bad cos θ samples diverge about equally across homologs"
    return name, _homolog_pairs(ctx, title, subtitle, bad_of, good_of, 2, "median RPD (%)")


# ---------------------------------------------------------------------------------------------
# 74  Which Aroclor is mistaken for which
# ---------------------------------------------------------------------------------------------
def fig_74(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    name, title = "74-costheta-mismatch-table.svg", "Aroclor misidentification transitions"
    if not ctx.bad_ct_keys:
        return _nothing(name, title, NO_BAD_COS)

    mismatch_from_to = Counter()
    for key, ct in ctx.costheta.items():
        if not ct["match_agrees"]:
            mismatch_from_to[(ct["best_manual"], ct["best_auto"])] += 1
    transitions = sorted(mismatch_from_to.items(), key=lambda x: -x[1])
    n_trans = min(len(transitions), 12)

    W, H = 480, 56 + n_trans * 28 + 20
    ML, MR, MT = 20, 20, 56

    svg = header(W, H)
    svg += txt(ML, 24, title, size=14, weight="700")
    svg += txt(ML, 40, f"{ctx.n_mismatch} samples where auto picks a different best aroclor", size=11, color=INK3)

    hx = ML + 20
    svg += txt(hx, MT + 6, "Manual", size=10, weight="600", fam=MONO)
    svg += txt(hx + 80, MT + 6, "→", size=10, anchor="middle")
    svg += txt(hx + 100, MT + 6, "Auto", size=10, weight="600", fam=MONO)
    svg += txt(hx + 200, MT + 6, "Count", size=10, weight="600", fam=MONO)
    svg += txt(hx + 260, MT + 6, "Bar", size=10, weight="600", fam=MONO, color=INK3)
    svg += f'<line x1="{ML}" y1="{MT + 12}" x2="{W - MR}" y2="{MT + 12}" stroke="{RULE}" stroke-width="1"/>\n'

    max_trans = transitions[0][1]
    for ti, ((frm, to), n) in enumerate(transitions[:n_trans]):
        ry = MT + 18 + ti * 28
        svg += txt(hx, ry + 14, f"Aroclor {frm}", size=10, fam=MONO)
        svg += txt(hx + 80, ry + 14, "→", size=10, anchor="middle", color=INK3)
        svg += txt(hx + 100, ry + 14, f"Aroclor {to}", size=10, fam=MONO)
        svg += txt(hx + 200, ry + 14, str(n), size=10, fam=MONO, weight="600")
        svg += f'<rect x="{hx + 260}" y="{ry + 4}" width="{n / max_trans * 150}" height="14" fill="{PURPLE}" opacity="0.6" rx="2"/>\n'

    svg += "</svg>"
    return name, svg


# ---------------------------------------------------------------------------------------------
# 75  Do the same congeners drive both kinds of disagreement?
# ---------------------------------------------------------------------------------------------
def fig_75(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    name, title = "75-driver-overlap.svg", "Driver congeners: mass error vs cos θ overlap"
    ct_driver_labels = {label for label, _ in cos_drivers(ctx)}
    if not ctx.top_drivers or not ct_driver_labels:
        return _nothing(name, title, NO_BAD_SAMPLES if not ctx.top_drivers else NO_BAD_COS)

    both = ctx.driver_labels & ct_driver_labels
    only_mass = ctx.driver_labels - ct_driver_labels
    only_ct = ct_driver_labels - ctx.driver_labels

    W, H = 520, 380
    ML, MR, MT = 30, 30, 56

    svg = header(W, H)
    svg += txt(ML, 24, title, size=14, weight="700")
    svg += txt(ML, 40, f"{len(both)} of the top {len(ctx.top_drivers)} mass-error drivers also drive cos θ disagreement",
               size=11, color=INK3)

    col_w = (W - ML - MR) / 3
    columns = [
        (f"Both ({len(both)})", PURPLE, sorted(both)),
        (f"Only mass error ({len(only_mass)})", BAD_COLOR, sorted(only_mass)),
        (f"Only cos θ ({len(only_ct)})", WARN, sorted(only_ct)),
    ]
    for ci, (hdr, color, items) in enumerate(columns):
        cx = ML + ci * col_w + col_w / 2
        svg += txt(cx, MT + 10, hdr, size=11, weight="600", anchor="middle", color=color)
        svg += f'<line x1="{ML + ci * col_w + 10}" y1="{MT + 16}" x2="{ML + (ci+1) * col_w - 10}" y2="{MT + 16}" stroke="{color}" stroke-width="1.5"/>\n'
        for ii, label in enumerate(items):
            hom = ctx.cong_homolog.get(label, "?")
            svg += txt(cx, MT + 34 + ii * 16, f"PCB-{label}", size=9, anchor="middle", fam=MONO, weight="500")
            svg += txt(cx + 50, MT + 34 + ii * 16, f"({hom})", size=8, anchor="start", color=HOMOLOG_COLORS.get(hom, INK3))

    svg += "</svg>"
    return name, svg


# ---------------------------------------------------------------------------------------------
# 76  Is the fingerprint less stable in small samples?
# ---------------------------------------------------------------------------------------------
def fig_76(ds: Dataset) -> tuple[str, str]:
    ctx = context(ds)
    W, H = 620, 420
    ML, MR, MT, MB = 64, 180, 56, 56
    PW = W - ML - MR
    PH = H - MT - MB

    svg = header(W, H)
    svg += txt(ML, 24, "Total sample mass vs cos θ RPD", size=14, weight="700")
    svg += txt(ML, 40, "each dot = one sample; color = best aroclor match agrees / disagrees", size=11, color=INK3)
    svg += f'<rect x="{ML}" y="{MT}" width="{PW}" height="{PH}" fill="white" rx="3"/>\n'

    scatter_data = [(info["manual"], ctx.costheta[k]["ct_rpd"], not ctx.costheta[k]["match_agrees"])
                    for k, info in ctx.sample_total.items() if k in ctx.costheta]
    masses_sc = [d[0] for d in scatter_data if d[0] > 0] or [1.0]
    rpds_sc = [d[1] for d in scatter_data] or [0.0]
    x_min_log = math.floor(math.log10(max(min(masses_sc), 1)))
    x_max_log = max(math.ceil(math.log10(max(masses_sc))), x_min_log + 1)
    y_max = min(max(rpds_sc) * 1.05, 50) or 1

    def xpos(mass):
        if mass <= 0:
            return ML
        lm = math.log10(max(mass, 10**x_min_log))
        return ML + (lm - x_min_log) / (x_max_log - x_min_log) * PW

    def ypos(value):
        return MT + PH - min(value, y_max) / y_max * PH

    for exp in range(x_min_log, x_max_log + 1):
        gx = xpos(10**exp)
        svg += f'<line x1="{gx}" y1="{MT}" x2="{gx}" y2="{MT + PH}" stroke="{RULE}" stroke-width="0.5" opacity="0.5"/>\n'
        lbl = f"{10**exp:,.0f}" if 10**exp < 10000 else f"{10**exp/1000:.0f}k"
        svg += txt(gx, MT + PH + 14, lbl, size=8, anchor="middle", color=INK3, fam=MONO)
    for yt in range(0, int(y_max) + 1, 5):
        gy = ypos(yt)
        svg += f'<line x1="{ML}" y1="{gy}" x2="{ML + PW}" y2="{gy}" stroke="{RULE}" stroke-width="0.5" opacity="0.5"/>\n'
        svg += txt(ML - 5, gy + 3, f"{yt}%", size=8, anchor="end", color=INK3, fam=MONO)

    # Agreeing samples first, so the mismatches are drawn on top of them.
    for mass, value, is_bad in sorted(scatter_data, key=lambda d: d[2]):
        color = BAD_COLOR if is_bad else GOOD_COLOR
        svg += (f'<circle cx="{xpos(mass)}" cy="{ypos(value)}" r="{"4.5" if is_bad else "3"}" fill="{color}"'
                f' opacity="{"0.7" if is_bad else "0.35"}"/>\n')

    svg += f'<line x1="{ML}" y1="{MT + PH}" x2="{ML + PW}" y2="{MT + PH}" stroke="{RULE}" stroke-width="1"/>\n'
    svg += f'<line x1="{ML}" y1="{MT}" x2="{ML}" y2="{MT + PH}" stroke="{RULE}" stroke-width="1"/>\n'
    svg += txt(ML + PW/2, MT + PH + 38, "total sample mass (ng, log scale)", size=10, anchor="middle", color=INK3)
    svg += _y_title(ML, MT, PH, "cos θ RPD (%)", offset=40)

    leg_x = ML + PW + 20
    leg_y2 = MT + 20
    svg += f'<circle cx="{leg_x}" cy="{leg_y2 + 5}" r="4.5" fill="{BAD_COLOR}" opacity="0.7"/>\n'
    svg += txt(leg_x + 10, leg_y2 + 9, "Mismatch", size=9, color=INK2)
    svg += txt(leg_x + 10, leg_y2 + 21, "(bad cos θ)", size=9, color=INK3)
    svg += f'<circle cx="{leg_x}" cy="{leg_y2 + 39}" r="3" fill="{GOOD_COLOR}" opacity="0.35"/>\n'
    svg += txt(leg_x + 10, leg_y2 + 43, "Match agrees", size=9, color=INK2)
    svg += txt(leg_x + 10, leg_y2 + 55, "(good)", size=9, color=INK3)
    svg += txt(leg_x, leg_y2 + 80, f"n = {sum(1 for d in scatter_data if d[2])} bad", size=9, fam=MONO, color=INK3)
    svg += txt(leg_x, leg_y2 + 93, f"n = {sum(1 for d in scatter_data if not d[2])} good", size=9, fam=MONO, color=INK3)

    svg += "</svg>"
    return "76-costheta-rpd-vs-mass.svg", svg


FIGURES = [fig_59, fig_60, fig_61, fig_62, fig_63, fig_64, fig_65, fig_66, fig_67, fig_68, fig_69,
           fig_70, fig_71, fig_72, fig_73, fig_74, fig_75, fig_76]
for _figure in (fig_64, fig_65, fig_66, fig_67, fig_68, fig_69, fig_70, fig_71, fig_72, fig_73, fig_74, fig_75, fig_76):
    _figure.needs_aroclors = True
