"""Figures 42 and 43: does automated quantification point at the same Aroclor as the analyst's?

A sample's congener masses are compared with each Aroclor's congener profile by cosine
similarity (cos theta). The Aroclor with the highest cos theta is the sample's best match.
The comparison needs a fingerprint on both sides, so a sample whose manual or automated
total is zero has no best match to compare.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass

from .data import Dataset
from .svg import ACCENT, INK, INK2, INK3, MONO, RED, RULE, header, txt

AROCLOR_COLORS = {
    "1221": "#8B5CF6", "1232": "#0E6B68", "1016": "#2B5D8C",
    "1242": "#A85723", "1248": "#C43D3D", "1254": "#059669",
    "1260": "#7C3AED", "1262": "#B45309", "1268": "#6366F1",
}


@dataclass(frozen=True)
class Match:
    """One sample's cos theta against every Aroclor, from the manual and the automated masses."""

    manual: dict[str, float]
    auto: dict[str, float]
    best_manual: str
    best_auto: str
    #: False when either total is zero: there is then no fingerprint, and no best match.
    defined: bool

    @property
    def agrees(self) -> bool:
        return self.best_manual == self.best_auto

    @property
    def cos_manual(self) -> float:
        return self.manual[self.best_manual]

    @property
    def cos_auto(self) -> float:
        return self.auto[self.best_auto]


def cos_theta(vec_a: dict[str, float], vec_b: dict[str, float], labels: list[str]) -> float:
    """Cosine similarity of two vectors indexed by congener label."""
    dot = 0.0
    mag_a = 0.0
    mag_b = 0.0
    for lbl in labels:
        a = vec_a.get(lbl, 0.0)
        b = vec_b.get(lbl, 0.0)
        dot += a * b
        mag_a += a * a
        mag_b += b * b
    if mag_a == 0 or mag_b == 0:
        return 0.0
    return dot / (math.sqrt(mag_a) * math.sqrt(mag_b))


def _matches(ds: Dataset) -> dict[tuple[str, str], Match]:
    manual = defaultdict(dict)
    auto = defaultdict(dict)
    for r in ds.sample_natives:
        key = (r["run"], r["sample"])
        manual[key][r["label"]] = r["_m"]
        auto[key][r["label"]] = r["_a"]

    reference = ds.aroclors[ds.aroclor_names[0]]
    labels = sorted(set(reference) & {lbl for vec in manual.values() for lbl in vec}, key=lambda x: (len(x), x))
    out = {}
    for key in sorted(manual):
        cos_m = {name: cos_theta(manual[key], ds.aroclors[name], labels) for name in ds.aroclor_names}
        cos_a = {name: cos_theta(auto[key], ds.aroclors[name], labels) for name in ds.aroclor_names}
        total_m = sum(manual[key].get(lbl, 0) for lbl in labels)
        total_a = sum(auto[key].get(lbl, 0) for lbl in labels)
        out[key] = Match(manual=cos_m, auto=cos_a,
                         best_manual=max(cos_m, key=cos_m.get), best_auto=max(cos_a, key=cos_a.get),
                         defined=total_m != 0 and total_a != 0)
    return out


def matches(ds: Dataset) -> dict[tuple[str, str], Match]:
    """Every sample's fingerprint match, keyed by (run, sample) and computed once per dataset."""
    return ds.memo("fingerprint_matches", lambda: _matches(ds))


def _circle(cx, cy, r, fill, opacity):
    return f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{fill}" opacity="{opacity}"/>\n'


def fig_42(ds: Dataset) -> tuple[str, str]:
    found = {key: m for key, m in matches(ds).items() if m.defined}

    MARGIN_L = 58
    MARGIN_R = 16
    MARGIN_T = 48
    MARGIN_B = 50
    CELL_W = 200
    CELL_H = 200
    GAP_X = 56
    GAP_Y = 56
    COLS = 3
    ROWS = 3
    W = MARGIN_L + COLS * CELL_W + (COLS - 1) * GAP_X + MARGIN_R
    H = MARGIN_T + ROWS * CELL_H + (ROWS - 1) * GAP_Y + MARGIN_B

    svg = header(W, H)
    svg += txt(MARGIN_L, 28, "Aroclor cos θ agreement: manual vs automated", size=16, weight="700")
    svg += txt(W - MARGIN_R, 28, "each dot = one sample; 1:1 line = perfect agreement",
               size=11, color=INK3, anchor="end")

    for idx, aname in enumerate(ds.aroclor_names):
        col = idx % COLS
        row = idx // COLS
        ox = MARGIN_L + col * (CELL_W + GAP_X)
        oy = MARGIN_T + row * (CELL_H + GAP_Y)
        color = AROCLOR_COLORS.get(aname, INK2)
        cms = [m.manual[aname] for m in found.values()]
        cas = [m.auto[aname] for m in found.values()]

        svg += f'<rect x="{ox}" y="{oy}" width="{CELL_W}" height="{CELL_H}" fill="white" rx="3"/>\n'

        def px(v, ox=ox):
            return ox + v * CELL_W

        def py(v, oy=oy):
            return oy + CELL_H - v * CELL_H

        for g in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
            svg += f'<line x1="{px(g)}" y1="{py(0)}" x2="{px(g)}" y2="{py(1)}" stroke="{RULE}" stroke-width="0.5"/>\n'
            svg += f'<line x1="{px(0)}" y1="{py(g)}" x2="{px(1)}" y2="{py(g)}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += (f'<line x1="{px(0)}" y1="{py(0)}" x2="{px(1)}" y2="{py(1)}" stroke="{INK3}"'
                f' stroke-width="1" stroke-dasharray="4,3"/>\n')

        for cm, ca in zip(cms, cas):
            svg += _circle(px(cm), py(ca), 2.5, color, 0.5)

        n = len(cms)
        r_val = 0
        med_diff = 0
        if n > 1:
            mean_m = sum(cms) / n
            mean_a = sum(cas) / n
            cov = sum((cms[i] - mean_m) * (cas[i] - mean_a) for i in range(n))
            var_m = sum((v - mean_m) ** 2 for v in cms)
            var_a = sum((v - mean_a) ** 2 for v in cas)
            r_val = cov / (math.sqrt(var_m) * math.sqrt(var_a)) if var_m > 0 and var_a > 0 else 0
            med_diff = sorted(abs(cms[i] - cas[i]) for i in range(n))[n // 2]

        svg += txt(ox + CELL_W / 2, oy - 6, f"Aroclor {aname}", size=12, weight="600", anchor="middle", color=color)

        for tv in [0, 0.5, 1.0]:
            lbl = f"{tv:.1f}" if tv == 0.5 else str(int(tv))
            if col == 0:
                svg += txt(ox - 4, py(tv) + 4, lbl, size=9, color=INK3, anchor="end", fam=MONO)
            if row == ROWS - 1:
                svg += txt(px(tv), oy + CELL_H + 13, lbl, size=9, color=INK3, anchor="middle", fam=MONO)

        stat_y = oy + CELL_H - 24
        svg += f'<rect x="{ox + CELL_W - 96}" y="{stat_y - 2}" width="94" height="26" fill="white" opacity="0.85" rx="2"/>\n'
        svg += txt(ox + CELL_W - 4, stat_y + 10, f"r = {r_val:.4f}", size=9, color=INK2, anchor="end", fam=MONO)
        svg += txt(ox + CELL_W - 4, stat_y + 22, f"med |Δ| = {med_diff:.4f}", size=9, color=INK2, anchor="end", fam=MONO)

        if col == 0:
            ry = oy + CELL_H / 2
            svg += txt(ox - 28, ry, "auto cos θ", size=10, color=INK2, anchor="middle", weight="500", rotate=-90)
        if row == ROWS - 1:
            svg += txt(ox + CELL_W / 2, oy + CELL_H + 28, "manual cos θ", size=10, color=INK2, anchor="middle", weight="500")

    svg += "</svg>"
    return "42-costheta-regression.svg", svg


def fig_43(ds: Dataset) -> tuple[str, str]:
    found = [m for m in matches(ds).values() if m.defined]
    n_agree = sum(1 for m in found if m.agrees)
    n_total = len(found)
    pct = 100 * n_agree / n_total if n_total else 0
    disagree = [m for m in found if not m.agrees]

    pair_counts = defaultdict(int)
    for m in disagree:
        pair_counts[tuple(sorted([m.best_manual, m.best_auto]))] += 1

    W = 720
    H = 560
    PLOT_L = 70
    PLOT_R = 260
    PLOT_T = 50
    PLOT_B = 60
    PW = W - PLOT_L - PLOT_R
    PH = H - PLOT_T - PLOT_B

    svg = header(W, H)
    svg += txt(PLOT_L, 30, "Aroclor best-match agreement", size=16, weight="700")
    svg += txt(PLOT_L, 44, "highest cos θ selects the aroclor; dot color = same answer or not",
               size=11, color=INK3)
    svg += f'<rect x="{PLOT_L}" y="{PLOT_T}" width="{PW}" height="{PH}" fill="white" rx="3"/>\n'

    def px(v):
        return PLOT_L + v * PW

    def py(v):
        return PLOT_T + PH - v * PH

    for g in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
        svg += f'<line x1="{px(g)}" y1="{py(0)}" x2="{px(g)}" y2="{py(1)}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += f'<line x1="{px(0)}" y1="{py(g)}" x2="{px(1)}" y2="{py(g)}" stroke="{RULE}" stroke-width="0.5"/>\n'
        svg += txt(PLOT_L - 6, py(g) + 4, f"{g:.1f}", size=9, color=INK3, anchor="end", fam=MONO)
        svg += txt(px(g), PLOT_T + PH + 14, f"{g:.1f}", size=9, color=INK3, anchor="middle", fam=MONO)
    svg += (f'<line x1="{px(0)}" y1="{py(0)}" x2="{px(1)}" y2="{py(1)}" stroke="{INK3}"'
            f' stroke-width="1" stroke-dasharray="4,3"/>\n')

    # Disagreements first, so the agreeing majority is drawn over them.
    for m in disagree:
        svg += _circle(px(m.cos_manual), py(m.cos_auto), 4, RED, 0.65)
    for m in found:
        if m.agrees:
            svg += _circle(px(m.cos_manual), py(m.cos_auto), 3.5, ACCENT, 0.5)

    svg += txt(PLOT_L + PW / 2, H - 8, "manual best-match cos θ", size=12, color=INK2, anchor="middle", weight="500")
    svg += txt(PLOT_L - 42, PLOT_T + PH / 2, "auto best-match cos θ", size=12, color=INK2, anchor="middle",
               weight="500", rotate=-90)

    sx = PLOT_L + PW + 24
    sy = PLOT_T + 10
    svg += txt(sx, sy, "Summary", size=13, weight="700")
    sy += 22
    svg += txt(sx, sy, f"{n_total}", size=22, weight="700", color=ACCENT)
    svg += txt(sx + 42, sy, "samples tested", size=11, color=INK2)
    sy += 28
    svg += f'<circle cx="{sx + 6}" cy="{sy - 4}" r="5" fill="{ACCENT}"/>\n'
    svg += txt(sx + 16, sy, f"Agree: {n_agree} ({pct:.1f}%)", size=11, color=INK, weight="600")
    sy += 18
    svg += f'<circle cx="{sx + 6}" cy="{sy - 4}" r="5" fill="{RED}"/>\n'
    svg += txt(sx + 16, sy, f"Disagree: {n_total - n_agree} ({100 - pct:.1f}%)", size=11, color=INK, weight="600")
    sy += 28

    svg += txt(sx, sy, "Disagreements by cos θ", size=11, weight="600", color=INK2)
    sy += 16
    bins_d = {"< 0.4": 0, "0.4–0.6": 0, "0.6–0.8": 0, "0.8–0.95": 0, "≥ 0.95": 0}
    for m in disagree:
        v = max(m.cos_manual, m.cos_auto)
        if v < 0.4:
            bins_d["< 0.4"] += 1
        elif v < 0.6:
            bins_d["0.4–0.6"] += 1
        elif v < 0.8:
            bins_d["0.6–0.8"] += 1
        elif v < 0.95:
            bins_d["0.8–0.95"] += 1
        else:
            bins_d["≥ 0.95"] += 1
    for label, cnt in bins_d.items():
        svg += txt(sx + 8, sy, f"{label}:", size=10, color=INK2, fam=MONO)
        svg += txt(sx + 100, sy, str(cnt), size=10, color=INK, anchor="end", fam=MONO, weight="600")
        sy += 15

    sy += 12
    svg += txt(sx, sy, "Top confusion pairs", size=11, weight="600", color=INK2)
    sy += 16
    for pair, cnt in sorted(pair_counts.items(), key=lambda x: -x[1])[:6]:
        svg += txt(sx + 8, sy, f"{pair[0]} / {pair[1]}", size=10, color=INK2, fam=MONO)
        svg += txt(sx + 180, sy, str(cnt), size=10, color=INK, anchor="end", fam=MONO, weight="600")
        sy += 15

    sy += 12
    low_disagree = sum(1 for m in disagree if max(m.cos_manual, m.cos_auto) < 0.8)
    svg += txt(sx, sy, f"{low_disagree} of {len(disagree)} disagree below 0.8", size=10, color=INK3)
    sy += 14
    svg += txt(sx, sy, "(weak fingerprints, less reliable)", size=10, color=INK3)

    svg += "</svg>"
    return "43-costheta-bestmatch.svg", svg


fig_42.needs_aroclors = True
fig_43.needs_aroclors = True
FIGURES = [fig_42, fig_43]
