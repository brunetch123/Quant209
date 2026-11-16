"""Figures 17, 18 and 26: agreement by boundary method, by flag, and for the worst congeners.

All three are scored on the detectable native peaks, the population the headline figures use:
native analytes whose manual area clears 1000 counts*s.
"""

from __future__ import annotations

import math
import statistics
from collections import Counter, defaultdict

from .data import Dataset
from .svg import (ACCENT, BLUE, HOMOLOGS, INK, INK2, INK3, MONO, PURPLE, RED, SANS, WARN,
                  header, line, txt)

MIN_AREA = 1000.0
BAND = 25.0
SURFACE = "#FFFFFF"
ACCENT_SOFT = "#E6F0EF"
RULE_SOFT = "#EAECEF"

#: The flag combinations the September figure showed, in its order, and their colours.
FLAG_ROWS = [("-", ACCENT), ("clustered", BLUE), ("shoulder", WARN), ("absent", RED),
             ("low_sn,shoulder", "#A85797"), ("absent,low_sn", RED)]
WORST_COLORS = [RED, "#C44D6C", WARN, "#D4874A", "#A85797"]
COUNT_WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five"}


def detectable(ds: Dataset) -> list[dict]:
    """Detectable native area rows, each with its absolute RPD in percent as ``_rpd``."""
    rows = []
    for r in ds.areas:
        if r["role"] != "native" or r["manual_f"] < MIN_AREA:
            continue
        mean = (r["auto_f"] + r["manual_f"]) / 2
        r["_rpd"] = abs(r["auto_f"] - r["manual_f"]) / mean * 100 if mean else 0.0
        rows.append(r)
    return rows


def canonical_flags(flags: str) -> str:
    """One spelling per flag set, whatever order the flags were raised in."""
    if not flags or flags == "-":
        return "-"
    return ",".join(sorted(flags.split(",")))


def _summary(rpds: list[float]) -> tuple[int, float, float]:
    """Rows, median |RPD| and share within the band, the last two to one decimal."""
    if not rpds:
        return 0, 0.0, 0.0
    return (len(rpds), round(statistics.median(rpds), 1),
            round(sum(v <= BAND for v in rpds) / len(rpds) * 100, 1))


def _text(x, y, s, font=SANS, size=12, fill=INK2, anchor="middle", weight=400):
    return txt(x, y, s, size=size, color=fill, anchor=anchor, fam=font, weight=weight)


def _box(x, y, w, h, fill, stroke, sw=1.2, rx=4):
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>\n'


def fig_17(ds: Dataset) -> tuple[str, str]:
    rows = detectable(ds)
    by_method = defaultdict(list)
    for r in rows:
        by_method[r["method"]].append(r["_rpd"])
    methods = [(name, *_summary(by_method[name]), col) for name, col in
               [("valley_drop", ACCENT), ("deconvolved", WARN), ("rt_window", BLUE), ("tangent_skim", RED)]]

    w, h = 650, 400
    s = header(w, h)
    s += _text(325, 24, "Integration Method Performance", size=14, weight=600, fill=INK)

    ox, oy = 60, 55
    bh = 24
    gap = 10
    mx_rpd = 100
    # Centred over the bars, the first heading sat on top of the 25% reference label and hid it.
    heading = {"anchor": "start"} if ds.extras else {}
    heading_x = ox if ds.extras else 220

    s += _text(heading_x, 46, "Median |RPD|", size=12, weight=600, fill=INK, **heading)
    for i, (name, n, rpd, cov, col) in enumerate(methods):
        y = oy + i * (bh + gap)
        bw = (min(rpd, mx_rpd) / mx_rpd) * 280
        s += f'<rect x="{ox + 100}" y="{y}" width="{max(bw, 2)}" height="{bh}" rx="2" fill="{col}" opacity="0.7"/>\n'
        s += _text(ox + 95, y + bh / 2 + 4, name, size=10, fill=INK, anchor="end", font=MONO, weight=500)
        s += _text(ox + 100 + bw + 6, y + bh / 2 + 4, f"{rpd:.1f}%", size=10, fill=INK, anchor="start", font=MONO)
        s += _text(ox + 430, y + bh / 2 + 4, f"n={n:,}", size=9, fill=INK3, anchor="start", font=MONO)

    ref_x = ox + 100 + (25 / mx_rpd) * 280
    s += line(ref_x, oy - 5, ref_x, oy + 4 * (bh + gap) - gap, WARN, 1, "3,2")
    s += _text(ref_x, oy - 8, "25%", size=8, fill=WARN, font=MONO)

    oy2 = 210
    s += _text(heading_x, oy2 - 10, "Within ±25%", size=12, weight=600, fill=INK, **heading)
    for i, (name, n, rpd, cov, col) in enumerate(methods):
        y = oy2 + i * (bh + gap)
        bw = (cov / 100) * 280
        s += f'<rect x="{ox + 100}" y="{y}" width="{max(bw, 2)}" height="{bh}" rx="2" fill="{col}" opacity="0.7"/>\n'
        s += _text(ox + 95, y + bh / 2 + 4, name, size=10, fill=INK, anchor="end", font=MONO, weight=500)
        s += _text(ox + 100 + bw + 6, y + bh / 2 + 4, f"{cov:.1f}%", size=10, fill=INK, anchor="start", font=MONO)

    ref_x2 = ox + 100 + (80 / 100) * 280
    s += line(ref_x2, oy2 - 5, ref_x2, oy2 + 4 * (bh + gap) - gap, ACCENT, 1, "3,2")
    s += _text(ref_x2, oy2 - 8, "80%", size=8, fill=ACCENT, font=MONO)

    s += "</svg>"
    return "17-method-performance.svg", s


def fig_18(ds: Dataset) -> tuple[str, str]:
    rows = detectable(ds)
    by_flag = defaultdict(list)
    for r in rows:
        by_flag[canonical_flags(r["flags"])].append(r["_rpd"])

    # The September rows keep their places; any left empty go to the commonest other combinations.
    shown = [(flag, col) for flag, col in FLAG_ROWS if by_flag.get(flag)]
    listed = {flag for flag, _ in FLAG_ROWS}
    others = sorted((f for f in by_flag if f not in listed), key=lambda f: (-len(by_flag[f]), f))
    shown += [(flag, INK2) for flag in others[:len(FLAG_ROWS) - len(shown)]]
    if ds.extras:
        # A flag the September engine did not have gets its own row rather than competing for
        # a September slot: ``forced`` marks every peak a named rule produced.
        taken = {flag for flag, _ in shown}
        shown += [(flag, PURPLE) for flag in others
                  if "forced" in flag.split(",") and flag not in taken and len(by_flag[flag]) >= 20]
    flags = [("(none)" if flag == "-" else flag, *_summary(by_flag[flag]), col) for flag, col in shown]

    w, h = 650, 365 + 40 * max(0, len(flags) - len(FLAG_ROWS))
    s = header(w, h)
    s += _text(325, 24, "Agreement by Detection Flag", size=14, weight=600, fill=INK)

    ox, oy = 50, 55
    bh = 30
    gap = 10
    pw = 330

    s += _text(ox + 140 + pw / 2, 45, "Fraction within ±25%", size=12, weight=600, fill=INK)
    for i, (name, n, rpd, cov, col) in enumerate(flags):
        y = oy + i * (bh + gap)
        bw = (cov / 100) * pw
        s += f'<rect x="{ox + 140}" y="{y}" width="{pw}" height="{bh}" rx="2" fill="{RULE_SOFT}"/>\n'
        s += f'<rect x="{ox + 140}" y="{y}" width="{max(bw, 2)}" height="{bh}" rx="2" fill="{col}" opacity="0.7"/>\n'
        s += _text(ox + 135, y + bh / 2 + 4, name, size=10, fill=INK, anchor="end", font=MONO, weight=500)
        inside = bw >= pw - 50
        vx = ox + 140 + bw - 6 if inside else ox + 140 + bw + 6
        s += _text(vx, y + bh / 2 + 4, f"{cov:.1f}%", size=10, fill=SURFACE if inside else INK,
                   anchor="end" if inside else "start", font=MONO)
        s += _text(ox + 140 + pw + 10, y + bh / 2 + 4, f"n={n:,}", size=9, fill=INK3, anchor="start", font=MONO)

    # What the flags do and do not explain, stated from the rows themselves.
    clean = [v for flag in ("-", "clustered") for v in by_flag.get(flag, [])]
    hard = [r["_rpd"] for r in rows if "shoulder" in r["flags"] or "absent" in r["flags"]]
    failures = sum(r["_rpd"] > BAND for r in rows)
    hard_share = len(hard) / len(rows) * 100 if rows else 0.0
    hard_failures = sum(v > BAND for v in hard) / failures * 100 if failures else 0.0
    floor = min((_summary(by_flag[flag])[2] for flag in ("-", "clustered") if by_flag.get(flag)), default=0.0)
    y_ann = oy + max(len(flags), len(FLAG_ROWS)) * (bh + gap) + 10
    s += _box(50, y_ann, 550, 45, fill=ACCENT_SOFT, stroke=ACCENT)
    if clean:
        s += _text(325, y_ann + 18,
                   f"Unflagged and clustered peaks agree at >{math.floor(floor)}%. Shoulder and absent", size=10, fill=INK2)
        s += _text(325, y_ann + 35,
                   f"peaks account for {hard_share:.0f}% of comparisons but {hard_failures:.0f}% of failures outside ±25%.",
                   size=10, fill=INK2)

    s += "</svg>"
    return "18-flag-performance.svg", s


def _nice_ceiling(value: float) -> float:
    """The next 1, 2 or 5 times a power of ten at or above ``value``."""
    if value <= 0:
        return 1.0
    magnitude = 10 ** math.floor(math.log10(value))
    for step in (1, 2, 5, 10):
        if value <= step * magnitude:
            return float(step * magnitude)
    return float(10 * magnitude)


def worst_congeners(ds: Dataset, count: int = 5, min_rows: int = 50) -> list[tuple]:
    """The congeners with the highest median |RPD| among those with enough detectable rows."""
    by_label = defaultdict(list)
    for r in detectable(ds):
        by_label[r["label"]].append(r)
    ranked = []
    for label, rows in by_label.items():
        if len(rows) < min_rows:
            continue
        near_zero = sum(1 for r in rows if r["auto_f"] <= 0.01 * r["manual_f"])
        ranked.append((statistics.median(r["_rpd"] for r in rows), len(rows), label, rows[0]["homolog"],
                       round(statistics.mean(r["auto_f"] for r in rows)),
                       round(statistics.mean(r["manual_f"] for r in rows)), near_zero))
    ranked.sort(key=lambda t: (-t[0], -t[1], t[2]))
    return ranked[:count]


def fig_26(ds: Dataset) -> tuple[str, str]:
    congeners = worst_congeners(ds)
    count_word = COUNT_WORDS.get(len(congeners), str(len(congeners))).capitalize()

    w, h = 700, 380
    s = header(w, h)
    s += _text(350, 24, f"Automated vs Manual Area: {count_word} Worst Congeners", size=14, weight=600, fill=INK)

    ox, oy = 80, 65
    bh = 40
    gap = 18
    pw = 420

    s += _text(ox + 150, 55, "Mean area (counts·s)", size=10, fill=INK3)
    max_val = _nice_ceiling(max((max(c[4], c[5]) for c in congeners), default=1.0))

    def label_of(value):
        return f"{value:,.0f}" if value < 1000000 else f"{value / 1000000:.1f}M"

    for i, (_median, n, label, _hom, mean_auto, mean_manual, nz) in enumerate(congeners):
        y = oy + i * (bh + gap)
        col = WORST_COLORS[i]
        s += _text(ox - 5, y + bh / 2 + 4, f"PCB {label}", size=11, fill=col, anchor="end", font=MONO, weight=600)

        m_bw = min((mean_manual / max_val) * pw, pw)
        s += f'<rect x="{ox}" y="{y}" width="{max(m_bw, 2)}" height="{bh / 2 - 1}" rx="2" fill="{INK2}" opacity="0.5"/>\n'
        s += _text(ox + m_bw + 5, y + bh / 4 + 3, f"manual: {label_of(mean_manual)}", size=8, fill=INK2, anchor="start", font=MONO)

        a_bw = min((mean_auto / max_val) * pw, pw)
        s += f'<rect x="{ox}" y="{y + bh / 2 + 1}" width="{max(a_bw, 2)}" height="{bh / 2 - 1}" rx="2" fill="{col}" opacity="0.7"/>\n'
        s += _text(ox + max(a_bw, 2) + 5, y + 3 * bh / 4 + 3, f"auto: {label_of(mean_auto)}", size=8, fill=col, anchor="start", font=MONO)

        s += _text(650, y + bh / 2 + 4, f"{nz}/{n} near zero", size=8, fill=INK3, anchor="end", font=MONO)

    s += _box(50, 360, 600, 8, fill=RED, stroke="none", sw=0, rx=2)
    if congeners:
        homologs = Counter(c[3] for c in congeners)
        if len(homologs) == 1:
            lead = f"All {count_word.lower()} are {next(iter(homologs))}chlorobiphenyls."
        else:
            lead = "Homologs: " + ", ".join(f"{homologs[hom]} {hom}" for hom in HOMOLOGS if homologs[hom]) + "."
        shares = [c[6] / c[1] * 100 for c in congeners]
        s += _text(350, 377,
                   f"{lead} The software returns near-zero for {min(shares):.0f}–{max(shares):.0f}% of comparisons.",
                   size=9, fill=INK3, font=SANS)

    s += "</svg>"
    return "26-worst5-distribution.svg", s


FIGURES = [fig_17, fig_18, fig_26]
