"""Figures 44 to 50b: which samples disagree on total mass, and what they have in common.

A sample is scored on its total native mass. Its RPD is ``|manual - auto| / mean``; a sample
at 25 % or more is called a bad sample throughout.
"""

from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from .data import Dataset
from .fingerprint import matches
from .svg import (ACCENT, BAD, BLUE, GOOD, HOMOLOG_COLORS, HOMOLOGS, INK, INK2, INK3, MARGINAL, MONO,
                  PURPLE, RED, RPD_BANDS, RULE, WARN, band_color, header, txt, upper_median)

BAD_RPD = 25.0
DATE_TOKEN = re.compile(r"^\d{6,8}$")


@dataclass
class Sample:
    """One sample's native congener rows and the totals drawn from them."""

    rows: list[dict] = field(default_factory=list)
    mass_m: float = 0.0
    mass_a: float = 0.0
    detected_m: int = 0
    detected_a: int = 0
    methods: Counter = field(default_factory=Counter)
    homolog_m: defaultdict = field(default_factory=lambda: defaultdict(float))

    @property
    def rpd(self) -> float:
        avg = (self.mass_m + self.mass_a) / 2
        return abs(self.mass_m - self.mass_a) / avg * 100 if avg > 0 else 200

    @property
    def bad(self) -> bool:
        return self.rpd >= BAD_RPD


def _samples(ds: Dataset) -> dict[tuple[str, str], Sample]:
    table: dict[tuple[str, str], Sample] = defaultdict(Sample)
    for r in ds.sample_natives:
        s = table[(r["run"], r["sample"])]
        s.rows.append(r)
        s.mass_m += r["_m"]
        s.mass_a += r["_a"]
        s.methods[r["method"]] += 1
        s.homolog_m[r["homolog"]] += r["_m"]
        if r["_m"] > 0:
            s.detected_m += 1
        if r["_det"]:
            s.detected_a += 1
    return dict(table)


def samples(ds: Dataset) -> dict[tuple[str, str], Sample]:
    """Every sample, keyed by (run, sample), in the order the comparison lists them."""
    return ds.memo("error_samples", lambda: _samples(ds))


def run_family(run: str) -> str:
    """The letters that open a run name's first two fields, which the batches of one study share."""
    fields = [f for f in run.split("_") if not DATE_TOKEN.match(f) and f.lower() != "rerun"]
    stems = [m.group(0) for m in (re.match(r"[A-Za-z]+", f) for f in fields[:2]) if m]
    return "_".join(stems) or run


def bad_sample_groups(ds: Dataset) -> list[tuple[str, str, str]]:
    """What the bad samples have in common, as (group, how many, detail), the largest group first.

    Three things can be read off the comparison itself: whole runs going wrong together, one
    side reporting no mass at all, and the internal standard being integrated differently.
    The groups overlap; a sample can sit in more than one.
    """
    table = samples(ds)
    bad = {key for key, s in table.items() if s.bad}
    if not bad:
        return []

    run_total = Counter(key[0] for key in table)
    run_bad = Counter(key[0] for key in bad)
    clustered = [run for run in run_bad if run_bad[run] >= 3 and run_bad[run] * 2 >= run_total[run]]
    family_bad, family_total = Counter(), Counter()
    for run, n in run_total.items():
        family_total[run_family(run)] += n
        family_bad[run_family(run)] += run_bad[run]
    top_family = max(family_bad, key=lambda f: (family_bad[f], f))

    auto_zero = sum(1 for key in bad if table[key].mass_a == 0)
    manual_zero = sum(1 for key in bad if table[key].mass_m == 0)

    is_error = defaultdict(list)
    for r in ds.areas:
        if r["file_type"] == "sample" and r["role"].startswith("IS_") and r["manual_f"] > 0:
            is_error[(r["run"], r["sample"])].append(abs((r["auto_f"] - r["manual_f"]) / r["manual_f"] * 100))
    is_off = {key for key, errs in is_error.items() if sum(errs) / len(errs) >= 5}

    n_bad = len(bad)
    groups = [
        (sum(run_bad[run] for run in clustered), "Run clustering",
         f"{sum(run_bad[run] for run in clustered)} of {n_bad} bad samples, {len(clustered)} runs",
         f"{top_family}: {family_bad[top_family]}/{family_total[top_family]} samples bad"),
        (auto_zero + manual_zero, "One side reports zero",
         f"{auto_zero + manual_zero} of {n_bad} bad samples",
         f"auto = 0 in {auto_zero}, manual = 0 in {manual_zero}"),
        (len(is_off & bad), "IS area errors",
         f"{len(is_off & bad)} of {n_bad} bad samples",
         f"{len(is_off & bad)} of {len(is_off)} with IS error ≥ 5% are bad"),
    ]
    groups.sort(key=lambda g: -g[0])
    return [(name, count, detail) for size, name, count, detail in groups if size > 0]


# ---------------------------------------------------------------------------------------------
# 44  Do the mass failures and the fingerprint failures overlap?
# ---------------------------------------------------------------------------------------------
def fig_44(ds: Dataset) -> tuple[str, str]:
    table = samples(ds)
    high_rpd = {key for key, s in table.items() if s.bad}
    cos_disagree = {key for key, m in matches(ds).items() if m.defined and not m.agrees}
    both = high_rpd & cos_disagree
    only_rpd = high_rpd - cos_disagree
    only_cos = cos_disagree - high_rpd
    all_good = len(table) - len(high_rpd | cos_disagree)
    cos_color, both_color = (WARN, PURPLE) if ds.extras else (BLUE, INK)

    W, H = 620, 340
    svg = header(W, H)
    svg += txt(24, 28, "Error mode overlap", size=16, weight="700")
    svg += txt(W - 16, 28, f"{len(table)} samples total", size=11, color=INK3, anchor="end")

    cx1, cx2 = 180, 300
    cy = 170
    r = 100
    svg += f'<circle cx="{cx1}" cy="{cy}" r="{r}" fill="{RED}" opacity="0.15" stroke="{RED}" stroke-width="1.5"/>\n'
    svg += f'<circle cx="{cx2}" cy="{cy}" r="{r}" fill="{cos_color}" opacity="0.15" stroke="{cos_color}" stroke-width="1.5"/>\n'

    svg += txt(cx1 - 30, cy - r - 12, "High RPD", size=12, weight="600", color=RED, anchor="middle")
    svg += txt(cx1 - 30, cy - r + 2, "(≥ 25%)", size=10, color=RED, anchor="middle")
    svg += txt(cx2 + 30, cy - r - 12, "Cos θ disagree", size=12, weight="600", color=cos_color, anchor="middle")

    svg += txt(cx1 - 40, cy + 5, str(len(only_rpd)), size=28, weight="700", color=RED, anchor="middle")
    mid_x = (cx1 + cx2) / 2
    svg += txt(mid_x, cy + 5, str(len(both)), size=28, weight="700", color=both_color, anchor="middle")
    svg += txt(cx2 + 40, cy + 5, str(len(only_cos)), size=28, weight="700", color=cos_color, anchor="middle")
    svg += txt(mid_x, cy + 24, "both", size=10, color=INK3, anchor="middle")

    share = 100 * all_good / len(table) if table else 0
    svg += txt((cx1 + cx2) / 2, cy + r + 24, f"{all_good} samples ({share:.0f}%) pass both tests",
               size=11, color=ACCENT, anchor="middle", weight="600")

    sx = 430
    sy = 60
    for name, count, detail in bad_sample_groups(ds):
        svg += txt(sx, sy, name, size=10, weight="600", color=INK2)
        svg += txt(sx, sy + 14, count, size=10, color=INK3)
        svg += txt(sx + 8, sy + 32, detail, size=10, color=INK3)
        sy += 50

    svg += "</svg>"
    return "44-error-overlap.svg", svg


# ---------------------------------------------------------------------------------------------
# 45  Bad samples by run
# ---------------------------------------------------------------------------------------------
def fig_45(ds: Dataset) -> tuple[str, str]:
    table = samples(ds)
    run_total = Counter()
    run_bad = Counter()
    for key, s in table.items():
        run_total[key[0]] += 1
        if s.bad:
            run_bad[key[0]] += 1

    runs = sorted(run_total, key=lambda r: (-run_bad[r] / run_total[r] if run_bad[r] > 0 else 0, r))
    show_runs = [r for r in runs if run_bad[r] > 0]
    show_runs += [r for r in runs if run_bad[r] == 0][:5]

    n = len(show_runs)
    BAR_H = 18
    GAP = 4
    MARGIN_L = 300
    MARGIN_R = 80
    MARGIN_T = 56
    MARGIN_B = 30
    BAR_MAX_W = 250
    W = MARGIN_L + BAR_MAX_W + MARGIN_R
    H = MARGIN_T + n * (BAR_H + GAP) + MARGIN_B
    largest = max(run_total.values(), default=1)

    svg = header(W, H)
    svg += txt(MARGIN_L, 28, "Samples with RPD ≥ 25% by run", size=16, weight="700")
    svg += txt(MARGIN_L, 44, "red = bad, gray = total in run", size=11, color=INK3)

    for i, run in enumerate(show_runs):
        y = MARGIN_T + i * (BAR_H + GAP)
        total = run_total[run]
        bad = run_bad[run]
        rate = bad / total if total > 0 else 0
        display_name = run if len(run) <= 38 else run[:35] + "..."

        svg += txt(MARGIN_L - 6, y + BAR_H - 4, display_name, size=9, color=INK2 if bad > 0 else INK3,
                   anchor="end", fam=MONO)
        tw = BAR_MAX_W * (total / largest)
        svg += f'<rect x="{MARGIN_L}" y="{y}" width="{tw}" height="{BAR_H}" fill="{RULE}" rx="2"/>\n'
        if bad > 0:
            bw = BAR_MAX_W * (bad / largest)
            svg += f'<rect x="{MARGIN_L}" y="{y}" width="{bw}" height="{BAR_H}" fill="{RED}" opacity="0.7" rx="2"/>\n'
            svg += txt(MARGIN_L + tw + 6, y + BAR_H - 4, f"{bad}/{total} ({100*rate:.0f}%)",
                       size=9, color=RED, fam=MONO, weight="600")
        else:
            svg += txt(MARGIN_L + tw + 6, y + BAR_H - 4, f"0/{total}", size=9, color=INK3, fam=MONO)

    svg += "</svg>"
    return "45-run-clustering.svg", svg


# ---------------------------------------------------------------------------------------------
# 46  Peaks detected, manual against automated
# ---------------------------------------------------------------------------------------------
def _rpd_legend(x, y, row_gap=16, dy=0):
    out = ""
    for i, (label, color) in enumerate(RPD_BANDS):
        out += f'<circle cx="{x + 6}" cy="{y + i * row_gap + dy}" r="4" fill="{color}"/>\n'
        out += txt(x + 16, y + i * row_gap + dy + 4, label, size=10, color=INK2)
    return out


def fig_46(ds: Dataset) -> tuple[str, str]:
    table = samples(ds)
    n_congeners = max((len(s.rows) for s in table.values()), default=1)

    W, H = 560, 500
    PL, PR, PT, PB = 64, 24, 50, 54
    PW = W - PL - PR
    PH = H - PT - PB

    svg = header(W, H)
    svg += txt(PL, 28, "Peak detection: manual vs automated", size=16, weight="700")
    svg += txt(W - PR, 28, f"of {n_congeners} native congeners per sample", size=11, color=INK3, anchor="end")
    svg += f'<rect x="{PL}" y="{PT}" width="{PW}" height="{PH}" fill="white" rx="3"/>\n'

    def px(v):
        return PL + (v / n_congeners) * PW

    def py(v):
        return PT + PH - (v / n_congeners) * PH

    for g in [0, 50, 100, 150, n_congeners]:
        svg += f'<line x1="{px(g)}" y1="{py(0)}" x2="{px(g)}" y2="{py(n_congeners)}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += f'<line x1="{px(0)}" y1="{py(g)}" x2="{px(n_congeners)}" y2="{py(g)}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(PL - 6, py(g) + 4, str(g), size=9, color=INK3, anchor="end", fam=MONO)
        svg += txt(px(g), PT + PH + 14, str(g), size=9, color=INK3, anchor="middle", fam=MONO)
    svg += (f'<line x1="{px(0)}" y1="{py(0)}" x2="{px(n_congeners)}" y2="{py(n_congeners)}" stroke="{INK3}"'
            f' stroke-width="1" stroke-dasharray="4,3"/>\n')

    # Worst band first, so the bad samples are not hidden under the good ones.
    order = {BAD: 0, MARGINAL: 1, GOOD: 2}
    points = sorted(((order[band_color(s.rpd)], px(s.detected_m), py(s.detected_a), band_color(s.rpd))
                     for s in table.values()), key=lambda p: p[0])
    for _, cx, cy, color in points:
        op = 0.6 if color == GOOD else 0.7
        r = 3 if color == GOOD else 4
        svg += f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{color}" opacity="{op}"/>\n'

    lx, ly = PL + 12, PT + PH - 60
    svg += f'<rect x="{lx - 4}" y="{ly - 14}" width="160" height="56" fill="white" opacity="0.9" rx="3"/>\n'
    svg += _rpd_legend(lx, ly)

    svg += txt(PL + PW / 2, H - 8, "peaks detected (manual)", size=12, color=INK2, anchor="middle", weight="500")
    svg += txt(PL - 40, PT + PH / 2, "peaks detected (auto)", size=12, color=INK2, anchor="middle", weight="500", rotate=-90)

    svg += "</svg>"
    return "46-detection-scatter.svg", svg


# ---------------------------------------------------------------------------------------------
# 47  The congeners that carry the mass error
# ---------------------------------------------------------------------------------------------
def fig_47(ds: Dataset) -> tuple[str, str]:
    congener_abs_err = defaultdict(float)
    congener_homolog = {}
    for s in samples(ds).values():
        if not s.bad:
            continue
        for c in s.rows:
            congener_abs_err[c["label"]] += abs(c["_a"] - c["_m"])
            congener_homolog[c["label"]] = c["homolog"]
    top = sorted(congener_abs_err.items(), key=lambda x: -x[1])[:15]
    max_err = top[0][1] if top else 1.0

    n = len(top)
    BAR_H = 24
    GAP = 6
    MARGIN_L = 130
    MARGIN_R = 100
    MARGIN_T = 56
    MARGIN_B = 40
    BAR_MAX_W = 300
    W = MARGIN_L + BAR_MAX_W + MARGIN_R
    H = MARGIN_T + n * (BAR_H + GAP) + MARGIN_B

    svg = header(W, H)
    svg += txt(MARGIN_L, 28, "Congeners driving mass error", size=16, weight="700")
    svg += txt(MARGIN_L, 44, "total |auto − manual| for RPD ≥ 25% samples", size=11, color=INK3)

    for i, (label, total_err) in enumerate(top):
        y = MARGIN_T + i * (BAR_H + GAP)
        hom = congener_homolog.get(label, "")
        color = HOMOLOG_COLORS.get(hom, INK2)
        svg += txt(MARGIN_L - 6, y + BAR_H - 6, f"PCB {label}", size=10, weight="600", color=color, anchor="end", fam=MONO)
        bw = BAR_MAX_W * (total_err / max_err) if max_err else 0
        svg += f'<rect x="{MARGIN_L}" y="{y}" width="{bw}" height="{BAR_H}" fill="{color}" opacity="0.6" rx="2"/>\n'
        svg += txt(MARGIN_L + bw + 6, y + BAR_H - 6, f"{total_err/1000:.1f}k ng  ({hom})", size=9, color=INK2, fam=MONO)

    svg += "</svg>"
    return "47-congener-error-drivers.svg", svg


# ---------------------------------------------------------------------------------------------
# 48  Does agreement depend on how much PCB the sample holds?
# ---------------------------------------------------------------------------------------------
def fig_48(ds: Dataset) -> tuple[str, str]:
    W, H = 600, 420
    PL, PR, PT, PB = 64, 24, 50, 54
    PW = W - PL - PR
    PH = H - PT - PB

    svg = header(W, H)
    svg += txt(PL, 28, "Sample mass vs quantification agreement", size=16, weight="700")
    svg += txt(W - PR, 28, "RPD = |manual − auto| / mean", size=11, color=INK3, anchor="end")
    svg += f'<rect x="{PL}" y="{PT}" width="{PW}" height="{PH}" fill="white" rx="3"/>\n'

    lo_m, hi_m = -1, 5.5
    rpd_max = 210

    def px(log_mass):
        return PL + (log_mass - lo_m) / (hi_m - lo_m) * PW

    def py(rpd):
        return PT + PH - (rpd / rpd_max) * PH

    sups = {2: "²", 3: "³", 4: "⁴", 5: "⁵"}
    for exp in range(0, 6):
        x = px(exp)
        svg += f'<line x1="{x}" y1="{py(0)}" x2="{x}" y2="{py(rpd_max)}" stroke="{RULE}" stroke-width="0.5"/>\n'
        label = "1" if exp == 0 else "10" if exp == 1 else f"10{sups[exp]}"
        svg += txt(x, PT + PH + 14, label, size=9, color=INK3, anchor="middle", fam=MONO)
    for rpd_line in [0, 25, 50, 100, 150, 200]:
        y = py(rpd_line)
        svg += f'<line x1="{px(lo_m)}" y1="{y}" x2="{px(hi_m)}" y2="{y}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(PL - 6, y + 4, f"{rpd_line}%", size=9, color=INK3, anchor="end", fam=MONO)

    svg += (f'<line x1="{px(lo_m)}" y1="{py(25)}" x2="{px(hi_m)}" y2="{py(25)}" stroke="{RED}" stroke-width="1"'
            f' stroke-dasharray="6,3" opacity="0.5"/>\n')
    svg += txt(px(hi_m) - 4, py(25) - 5, "25% RPD threshold", size=9, color=RED, anchor="end", weight="500")

    order = {BAD: 0, MARGINAL: 1, GOOD: 2}
    points = []
    for s in samples(ds).values():
        total = max(s.mass_m, s.mass_a)
        if total <= 0:
            continue
        color = band_color(s.rpd)
        points.append((order[color], px(math.log10(max(total, 0.1))), py(min(s.rpd, rpd_max)), color))
    points.sort(key=lambda p: p[0])
    for _, cx, cy, color in points:
        op = 0.5 if color == GOOD else 0.7
        r = 3 if color == GOOD else 3.5
        svg += f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{color}" opacity="{op}"/>\n'

    svg += txt(PL + PW / 2, H - 6, "total PCB mass (ng, log scale)", size=12, color=INK2, anchor="middle", weight="500")
    svg += txt(PL - 46, PT + PH / 2, "RPD (%)", size=12, color=INK2, anchor="middle", weight="500", rotate=-90)

    # The legend first sat top right, over the samples that disagree by 200 % at high mass.
    lx, ly = (PL + 12 if ds.extras else PL + PW - 160), PT + 14
    svg += f'<rect x="{lx - 4}" y="{ly - 8}" width="160" height="56" fill="white" opacity="0.9" rx="3"/>\n'
    svg += _rpd_legend(lx, ly, dy=4)

    svg += "</svg>"
    return "48-mass-vs-rpd.svg", svg


# ---------------------------------------------------------------------------------------------
# 49  Boundary methods in good and in bad samples
# ---------------------------------------------------------------------------------------------
def fig_49(ds: Dataset) -> tuple[str, str]:
    methods = ["valley_drop", "deconvolved", "rt_window", "tangent_skim"]
    method_labels = ["valley drop", "deconvolved", "RT window", "tangent skim"]
    colors = {"valley_drop": ACCENT, "deconvolved": WARN, "rt_window": BLUE, "tangent_skim": RED}

    good_counts = Counter()
    bad_counts = Counter()
    for s in samples(ds).values():
        (bad_counts if s.bad else good_counts).update(s.methods)
    good_total = sum(good_counts.values())
    bad_total = sum(bad_counts.values())

    W, H = 540, 320
    MARGIN_L = 120
    MARGIN_T = 50
    MARGIN_B = 40
    BAR_W = 140
    GAP = 60

    svg = header(W, H)
    svg += txt(24, 28, "Integration method mix: good vs poor samples", size=14, weight="700")

    bar_x1 = MARGIN_L + 20
    bar_x2 = bar_x1 + BAR_W + GAP
    bar_top = MARGIN_T + 20
    bar_h = H - MARGIN_T - MARGIN_B - 20

    for bx, title, counts, total in [
        (bar_x1, "RPD < 25%", good_counts, good_total),
        (bar_x2, "RPD ≥ 25%", bad_counts, bad_total),
    ]:
        svg += txt(bx + BAR_W / 2, bar_top - 16, title, size=12, weight="600", anchor="middle", color=INK)
        svg += txt(bx + BAR_W / 2, bar_top - 2, f"({total:,} peaks)", size=10, anchor="middle", color=INK3)
        cy = bar_top
        for m in methods:
            frac = counts[m] / total if total > 0 else 0
            seg_h = frac * bar_h
            svg += f'<rect x="{bx}" y="{cy}" width="{BAR_W}" height="{seg_h}" fill="{colors[m]}" opacity="0.75"/>\n'
            if frac >= 0.04:
                svg += txt(bx + BAR_W / 2, cy + seg_h / 2 + 4, f"{100*frac:.1f}%", size=10, weight="600",
                           color="white", anchor="middle")
            cy += seg_h

    ly = H - 20
    lx = 24
    for i, (m, label) in enumerate(zip(methods, method_labels)):
        x = lx + i * 128
        svg += f'<rect x="{x}" y="{ly - 8}" width="12" height="12" fill="{colors[m]}" opacity="0.75" rx="2"/>\n'
        svg += txt(x + 18, ly + 3, label, size=10, color=INK2)

    svg += "</svg>"
    return "49-method-mix.svg", svg


# ---------------------------------------------------------------------------------------------
# 50  Are the bad samples richer in penta and hexa congeners?
# ---------------------------------------------------------------------------------------------
@dataclass
class Composition:
    """Penta-plus-hexa share and homolog make-up of the good and the bad samples."""

    points: list[tuple[float, float]]            # (penta+hexa fraction, RPD) per sample
    good_fraction: list[float]
    bad_fraction: list[float]
    good_medians: dict[str, float]
    bad_medians: dict[str, float]
    n_good: int
    n_bad: int

    @property
    def z(self) -> float:
        """Mann-Whitney z for the bad samples' penta+hexa share against the good samples'."""
        nb, ng = len(self.bad_fraction), len(self.good_fraction)
        if not nb or not ng:
            return 0.0
        ranked = sorted([(f, 1) for f in self.bad_fraction] + [(f, 0) for f in self.good_fraction])
        rank_sum = sum(rank for rank, (_, is_bad) in enumerate(ranked, 1) if is_bad)
        u = rank_sum - nb * (nb + 1) / 2
        return (u - nb * ng / 2) / math.sqrt(nb * ng * (nb + ng + 1) / 12)

    def statistic_text(self) -> str:
        p = math.erfc(abs(self.z) / math.sqrt(2))
        return f"Mann-Whitney U:  z = {self.z:.2f},  " + ("p < 0.0001" if p < 0.0001 else f"p = {p:.4f}")

    def visible_homologs(self) -> list[str]:
        """Homologs that make up at least 1 % of either bar."""
        shown = []
        for hom in HOMOLOGS:
            shares = [med[hom] / sum(med.values()) for med in (self.good_medians, self.bad_medians) if sum(med.values()) > 0]
            if any(share >= 0.01 for share in shares):
                shown.append(hom)
        return shown


def composition(ds: Dataset) -> Composition:
    points, good_frac, bad_frac = [], [], []
    good_hom, bad_hom = defaultdict(list), defaultdict(list)
    for s in samples(ds).values():
        if s.mass_m <= 0:
            continue
        ph_frac = (s.homolog_m.get("penta", 0) + s.homolog_m.get("hexa", 0)) / s.mass_m
        points.append((ph_frac, s.rpd))
        (bad_frac if s.bad else good_frac).append(ph_frac)
        for h in HOMOLOGS:
            (bad_hom if s.bad else good_hom)[h].append(s.homolog_m.get(h, 0) / s.mass_m)
    return Composition(
        points=points, good_fraction=good_frac, bad_fraction=bad_frac,
        good_medians={h: upper_median(good_hom[h]) if good_hom[h] else 0.0 for h in HOMOLOGS},
        bad_medians={h: upper_median(bad_hom[h]) if bad_hom[h] else 0.0 for h in HOMOLOGS},
        n_good=len(good_frac), n_bad=len(bad_frac),
    )


def _composition_points(comp: Composition, px, py, y_max) -> str:
    order = {BAD: 0, MARGINAL: 1, GOOD: 2}
    points = sorted(((order[band_color(rpd)], px(frac), py(min(rpd, y_max)), band_color(rpd))
                     for frac, rpd in comp.points), key=lambda p: p[0])
    out = ""
    for _, cx, cy, color in points:
        op = 0.45 if color == GOOD else 0.7
        r = 3 if color == GOOD else 3.5
        out += f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{color}" opacity="{op}"/>\n'
    return out


def _median_text(values: list[float]) -> str:
    return f"{upper_median(values) * 100:.1f}%" if values else "n/a"


def fig_50(ds: Dataset) -> tuple[str, str]:
    comp = composition(ds)

    W, H = 780, 636
    SCATTER_L = 70
    SCATTER_W = 340
    SCATTER_T = 56
    SCATTER_B = 74
    BAR_L = 500
    BAR_W = 120
    BAR_GAP = 40
    PH = H - SCATTER_T - SCATTER_B

    svg = header(W, H)
    svg += txt(SCATTER_L, 28, "Penta + hexa mass fraction vs quantification error", size=16, weight="700")
    svg += txt(SCATTER_L, 44, comp.statistic_text(), size=11, color=RED, weight="600", fam=MONO)
    svg += txt(SCATTER_L + 268, 44, f"good {_median_text(comp.good_fraction)}", size=11, color=ACCENT, fam=MONO)
    svg += txt(SCATTER_L + 340, 44, f"bad {_median_text(comp.bad_fraction)}", size=11, color=RED, fam=MONO)

    svg += f'<rect x="{SCATTER_L}" y="{SCATTER_T}" width="{SCATTER_W}" height="{PH}" fill="white" rx="3"/>\n'
    x_lo, x_hi = 0.0, 1.0
    y_max = 210

    def px(v):
        return SCATTER_L + (v - x_lo) / (x_hi - x_lo) * SCATTER_W

    def py(v):
        return SCATTER_T + PH - (v / y_max) * PH

    for xg in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
        svg += f'<line x1="{px(xg)}" y1="{py(0)}" x2="{px(xg)}" y2="{py(y_max)}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(px(xg), SCATTER_T + PH + 14, f"{xg:.1f}", size=9, color=INK3, anchor="middle", fam=MONO)
    for yg in [0, 25, 50, 100, 150, 200]:
        svg += f'<line x1="{px(x_lo)}" y1="{py(yg)}" x2="{px(x_hi)}" y2="{py(yg)}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(SCATTER_L - 6, py(yg) + 4, f"{yg}%", size=9, color=INK3, anchor="end", fam=MONO)
    svg += (f'<line x1="{px(x_lo)}" y1="{py(25)}" x2="{px(x_hi)}" y2="{py(25)}" stroke="{RED}" stroke-width="1"'
            f' stroke-dasharray="6,3" opacity="0.4"/>\n')
    svg += _composition_points(comp, px, py, y_max)

    svg += txt(SCATTER_L + SCATTER_W / 2, SCATTER_T + PH + 30, "penta + hexa fraction of total mass (manual)",
               size=11, color=INK2, anchor="middle", weight="500")
    svg += txt(SCATTER_L - 46, SCATTER_T + PH / 2, "RPD (%)", size=11, color=INK2, anchor="middle", weight="500", rotate=-90)

    svg += txt(BAR_L + BAR_W + BAR_GAP / 2, SCATTER_T - 10, "Homolog composition", size=13, weight="700", anchor="middle")
    for bx, medians, label, sublabel in [
        (BAR_L, comp.good_medians, "RPD < 25%", f"n = {comp.n_good}"),
        (BAR_L + BAR_W + BAR_GAP, comp.bad_medians, "RPD ≥ 25%", f"n = {comp.n_bad}"),
    ]:
        svg += txt(bx + BAR_W / 2, SCATTER_T + PH + 16, label, size=11, weight="600", anchor="middle", color=INK)
        svg += txt(bx + BAR_W / 2, SCATTER_T + PH + 30, sublabel, size=9, anchor="middle", color=INK3)
        total_med = sum(medians.values())
        cy = SCATTER_T
        for h in HOMOLOGS:
            frac = medians[h] / total_med if total_med > 0 else 0
            seg_h = frac * PH
            svg += f'<rect x="{bx}" y="{cy}" width="{BAR_W}" height="{seg_h}" fill="{HOMOLOG_COLORS[h]}" opacity="0.75"/>\n'
            if frac >= 0.04:
                svg += txt(bx + BAR_W / 2, cy + seg_h / 2 + 4, f"{100*frac:.0f}%", size=9, weight="600",
                           color="white", anchor="middle")
            cy += seg_h

    plot_bot = SCATTER_T + PH
    visible_hom = comp.visible_homologs()
    leg_y = plot_bot + 48
    spacing = 75
    for i, h in enumerate(visible_hom):
        lx = SCATTER_L + i * spacing
        svg += f'<rect x="{lx}" y="{leg_y - 7}" width="10" height="10" fill="{HOMOLOG_COLORS[h]}" opacity="0.75" rx="1"/>\n'
        svg += txt(lx + 14, leg_y + 3, h, size=9, color=INK2)
    if len(visible_hom) < len(HOMOLOGS):
        svg += txt(SCATTER_L + len(visible_hom) * spacing, leg_y + 3, "(others < 1%)", size=9, color=INK3)

    rpd_leg_y = plot_bot + 64
    for i, (label, color) in enumerate(RPD_BANDS):
        rx = SCATTER_L + i * 130
        svg += f'<circle cx="{rx + 6}" cy="{rpd_leg_y - 2}" r="4" fill="{color}"/>\n'
        svg += txt(rx + 16, rpd_leg_y + 2, label, size=9, color=INK2)

    svg += "</svg>"
    return "50-composition-vs-rpd.svg", svg


def fig_50a(ds: Dataset) -> tuple[str, str]:
    comp = composition(ds)

    PLOT = 400
    ML, MR = 60, 20
    MT = 56
    MB = 64
    W1 = ML + PLOT + MR
    H1 = MT + PLOT + MB

    svg = header(W1, H1)
    svg += txt(ML, 24, "Penta + hexa fraction vs quantification error", size=14, weight="700")
    svg += txt(ML, 40, comp.statistic_text(), size=10, color=RED, weight="600", fam=MONO)
    svg += txt(ML + 252, 40, f"good {_median_text(comp.good_fraction)}", size=10, color=ACCENT, fam=MONO)
    svg += txt(ML + 330, 40, f"bad {_median_text(comp.bad_fraction)}", size=10, color=RED, fam=MONO)
    svg += f'<rect x="{ML}" y="{MT}" width="{PLOT}" height="{PLOT}" fill="white" rx="3"/>\n'

    x_lo, x_hi = 0.0, 1.0
    y_max = 210

    def px(v):
        return ML + (v - x_lo) / (x_hi - x_lo) * PLOT

    def py(v):
        return MT + PLOT - (v / y_max) * PLOT

    for xg in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
        svg += f'<line x1="{px(xg)}" y1="{py(0)}" x2="{px(xg)}" y2="{py(y_max)}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(px(xg), MT + PLOT + 14, f"{xg:.1f}", size=9, color=INK3, anchor="middle", fam=MONO)
    for yg in [0, 25, 50, 100, 150, 200]:
        svg += f'<line x1="{px(x_lo)}" y1="{py(yg)}" x2="{px(x_hi)}" y2="{py(yg)}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(ML - 6, py(yg) + 4, f"{yg}%", size=9, color=INK3, anchor="end", fam=MONO)
    svg += (f'<line x1="{px(x_lo)}" y1="{py(25)}" x2="{px(x_hi)}" y2="{py(25)}" stroke="{RED}" stroke-width="1"'
            f' stroke-dasharray="6,3" opacity="0.4"/>\n')
    svg += _composition_points(comp, px, py, y_max)

    svg += txt(ML + PLOT / 2, MT + PLOT + 30, "penta + hexa fraction of total mass (manual)",
               size=11, color=INK2, anchor="middle", weight="500")
    svg += txt(ML - 40, MT + PLOT / 2, "RPD (%)", size=11, color=INK2, anchor="middle", weight="500", rotate=-90)

    leg_y = MT + PLOT + 50
    for i, (label, color) in enumerate(RPD_BANDS):
        rx = ML + i * 140
        svg += f'<circle cx="{rx + 6}" cy="{leg_y - 2}" r="4" fill="{color}"/>\n'
        svg += txt(rx + 16, leg_y + 2, label, size=9, color=INK2)

    svg += "</svg>"
    return "50a-composition-scatter.svg", svg


def fig_50b(ds: Dataset) -> tuple[str, str]:
    comp = composition(ds)

    BAR_W = 120
    BAR_GAP = 50
    ML2 = 30
    MR2 = 30
    MT2 = 48
    MB2 = 62
    BAR_H = 380
    W2 = ML2 + BAR_W + BAR_GAP + BAR_W + MR2
    H2 = MT2 + BAR_H + MB2

    svg = header(W2, H2)
    svg += txt(W2 / 2, 24, "Median homolog composition", size=14, weight="700", anchor="middle")
    svg += txt(W2 / 2, 40, "good (RPD < 25%) vs bad (RPD ≥ 25%) samples", size=11, color=INK3, anchor="middle")

    for bx, medians, label, sublabel in [
        (ML2, comp.good_medians, "RPD < 25%", f"n = {comp.n_good}"),
        (ML2 + BAR_W + BAR_GAP, comp.bad_medians, "RPD ≥ 25%", f"n = {comp.n_bad}"),
    ]:
        svg += txt(bx + BAR_W / 2, MT2 + BAR_H + 16, label, size=11, weight="600", anchor="middle", color=INK)
        svg += txt(bx + BAR_W / 2, MT2 + BAR_H + 30, sublabel, size=9, anchor="middle", color=INK3)
        total_med = sum(medians.values())
        cy = MT2
        for h in HOMOLOGS:
            frac = medians[h] / total_med if total_med > 0 else 0
            seg_h = frac * BAR_H
            svg += f'<rect x="{bx}" y="{cy}" width="{BAR_W}" height="{seg_h}" fill="{HOMOLOG_COLORS[h]}" opacity="0.75"/>\n'
            if frac >= 0.04:
                svg += txt(bx + BAR_W / 2, cy + seg_h / 2 + 4, f"{100*frac:.0f}%", size=10, weight="600",
                           color="white", anchor="middle")
            cy += seg_h

    visible_hom = comp.visible_homologs()
    leg_y = MT2 + BAR_H + 48
    spacing = (W2 - 2 * ML2) / (len(visible_hom) + 0.5) if visible_hom else 0
    for i, h in enumerate(visible_hom):
        lx = ML2 + i * spacing
        svg += f'<rect x="{lx}" y="{leg_y - 7}" width="10" height="10" fill="{HOMOLOG_COLORS[h]}" opacity="0.75" rx="1"/>\n'
        svg += txt(lx + 14, leg_y + 3, h, size=9, color=INK2)

    svg += "</svg>"
    return "50b-composition-bars.svg", svg


fig_44.needs_aroclors = True
FIGURES = [fig_44, fig_45, fig_46, fig_47, fig_48, fig_49, fig_50, fig_50a, fig_50b]
