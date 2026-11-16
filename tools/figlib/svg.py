"""The palette, the type, and the handful of SVG primitives every validation figure shares.

A colour means the same thing in every figure. A boundary method, a homolog and an agreement
band each keep one colour throughout, so a reader who has learned one figure can read the next.
"""

from __future__ import annotations

BG = "#F8F9FB"
INK = "#1A1D21"
INK2 = "#4A5058"
INK3 = "#7A828A"
RULE = "#D4D9DE"
ACCENT = "#0E6B68"
WARN = "#A85723"
BLUE = "#2B5D8C"
RED = "#C43D3D"
PURPLE = "#7C3AED"
SANS = '"IBM Plex Sans", sans-serif'
MONO = '"IBM Plex Mono", monospace'

#: Boundary methods, in the order the legends list them.
METHOD_COLORS = {
    "valley_drop": ACCENT,
    "rt_window": BLUE,
    "deconvolved": WARN,
    "tangent_skim": RED,
}
METHOD_LABELS = {
    "valley_drop": "valley drop",
    "rt_window": "RT window",
    "deconvolved": "deconvolved",
    "tangent_skim": "tangent skim",
}

HOMOLOGS = ["mono", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona", "deca"]
HOMOLOG_COLORS = {
    "mono": "#94A3B8", "di": "#6366F1", "tri": "#2B5D8C",
    "tetra": "#0E6B68", "penta": "#059669", "hexa": "#A85723",
    "hepta": "#C43D3D", "octa": "#7C3AED", "nona": "#B45309", "deca": "#4A5058",
}

#: Agreement bands on a sample's total mass: under 10 %, 10 to 25 %, and 25 % or more.
GOOD, MARGINAL, BAD = ACCENT, WARN, RED
RPD_BANDS = [("RPD < 10%", GOOD), ("RPD 10–25%", MARGINAL), ("RPD ≥ 25%", BAD)]


def band_color(rpd: float) -> str:
    """The agreement-band colour of a sample whose total-mass RPD is ``rpd`` percent."""
    if rpd >= 25:
        return BAD
    if rpd >= 10:
        return MARGINAL
    return GOOD


def esc(s: object) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def header(w: float, h: float) -> str:
    """Open an SVG of ``w`` by ``h`` pixels on the house background."""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">\n'
            f'<rect width="{w}" height="{h}" fill="{BG}" rx="6"/>\n')


def txt(x, y, s, size=11, color=INK, anchor="start", fam=SANS, weight="400", rotate=None) -> str:
    """One line of text. ``rotate`` turns it about its own anchor point, in degrees."""
    rot = f' transform="rotate({rotate},{x},{y})"' if rotate else ""
    return (f'<text x="{x}" y="{y}" text-anchor="{anchor}" fill="{color}"'
            f' font-family=\'{fam}\' font-size="{size}" font-weight="{weight}"{rot}>'
            f'{esc(s)}</text>\n')


def line(x1, y1, x2, y2, color=RULE, width=0.5, dash="", opacity=None) -> str:
    extra = f' stroke-dasharray="{dash}"' if dash else ""
    if opacity is not None:
        extra += f' opacity="{opacity}"'
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}"{extra}/>\n'


def swatch(x, y, color, size=10, opacity=None, rx=2) -> str:
    """A square legend key with its top-left corner at ``(x, y)``."""
    op = f' opacity="{opacity}"' if opacity is not None else ""
    return f'<rect x="{x}" y="{y}" width="{size}" height="{size}" rx="{rx}" fill="{color}"{op}/>\n'


def upper_median(values: list[float]) -> float:
    """The element at index ``len // 2`` of the sorted values: the median, or the upper of its two."""
    ordered = sorted(values)
    return ordered[len(ordered) // 2]
