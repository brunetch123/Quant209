"""Draw the Quant 209 mark and write every asset the application and the docs need.

The mark is a skeletal formula of PCB 77, 3,3',4,4'-tetrachlorobiphenyl. Two aromatic rings
share a single bond, each ring carries three inner lines for its double bonds, and four
chlorines sit on the meta and para positions. Both carbons either side of the ring bond are
left bare, which is what keeps the middle of the mark open when it is scaled down.

Detail is removed in a fixed order as the mark gets smaller, so every size is the same drawing
with less of it rather than a different drawing. The lockup keeps the molecule horizontal.
Square tiles turn it 30 degrees anticlockwise, because a biphenyl lying flat is three times
wider than it is tall and would otherwise shrink to a sliver in a square.

    python tools/make_logo.py

Writes into ``quant209/app/assets/`` and ``docs/images/``. Needs PySide6 for rendering and
Pillow for the multi-resolution Windows icon, both of which come with the ``app`` and
``report`` extras.
"""

from __future__ import annotations

import argparse
import io
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

TEAL = "#0F6E6B"
TEAL_DARK = "#43A29E"
INK = "#14181A"
PAPER = "#FBFCFC"

#: Ring geometry. The circumradius equals the carbon to carbon bond length, so the hexagon
#: side and the bond joining the two rings are the same length.
R = 16.0
CY = 46.0
CX_L, CX_R = 56.0, 104.0

#: How far an aromatic inner line sits inside its edge, and how much it is trimmed at each
#: end. At these values the inner line runs 74 per cent of the bond, which reads as a double
#: bond rather than as a dash.
INSET = 3.5
TRIM = 2.1

CL_STICK = 7.0
CL_GAP = 11.5

#: 3,3',4,4'-tetrachlorobiphenyl, as IUPAC ring positions on each ring.
PCB77_LEFT = [3, 4]
PCB77_RIGHT = [3, 4]

ICON_ROTATION = -30.0

#: Position number to vertex index. The left ring points its ipso carbon right and the right
#: ring points its own left, so the two rings mirror about the central bond.
LEFT_POS = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5}
RIGHT_POS = {1: 3, 2: 2, 3: 1, 4: 0, 5: 5, 6: 4}


def vertices(cx: float) -> list[tuple[float, float]]:
    return [(cx + R * math.cos(math.radians(60 * k)),
             CY + R * math.sin(math.radians(60 * k))) for k in range(6)]


def aromatic_lines(cx: float) -> list[tuple[float, float, float, float]]:
    """Three inner lines on alternating edges, which is one Kekule structure of the ring."""
    v = vertices(cx)
    out = []
    for e in (0, 2, 4):
        p1, p2 = v[e], v[(e + 1) % 6]
        mx, my = (p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2
        nx, ny = cx - mx, CY - my
        n = math.hypot(nx, ny)
        nx, ny = nx / n, ny / n
        ex, ey = p2[0] - p1[0], p2[1] - p1[1]
        el = math.hypot(ex, ey)
        ex, ey = ex / el, ey / el
        out.append((p1[0] + nx * INSET + ex * TRIM, p1[1] + ny * INSET + ey * TRIM,
                    p2[0] + nx * INSET - ex * TRIM, p2[1] + ny * INSET - ey * TRIM))
    return out


def chlorine(cx: float, vi: int):
    """The bond drawn out to a chlorine, and where its label sits."""
    v = vertices(cx)[vi]
    ux, uy = (v[0] - cx) / R, (v[1] - CY) / R
    stick = (v[0], v[1], v[0] + ux * CL_STICK, v[1] + uy * CL_STICK)
    lx, ly = v[0] + ux * CL_GAP, v[1] + uy * CL_GAP
    anchor = "start" if ux > 0.4 else "end" if ux < -0.4 else "middle"
    return stick, (lx, ly, anchor)


def f(x: float) -> str:
    return f"{x:.2f}".rstrip("0").rstrip(".")


def molecule(*, stroke: str, width: float, aromatic: bool, sticks: bool, labels: bool,
             label_size: float = 12.0) -> str:
    """The mark itself, as an SVG group, drawn flat."""
    left = PCB77_LEFT if (sticks or labels) else []
    right = PCB77_RIGHT if (sticks or labels) else []
    parts = []

    for cx in (CX_L, CX_R):
        parts.append('<polygon points="'
                     + " ".join(f"{f(x)},{f(y)}" for x, y in vertices(cx)) + '"/>')
    parts.append(f'<line x1="{f(CX_L + R)}" y1="{f(CY)}" x2="{f(CX_R - R)}" y2="{f(CY)}"/>')

    if aromatic:
        seg = "".join(f'<line x1="{f(a)}" y1="{f(b)}" x2="{f(c)}" y2="{f(d)}"/>'
                      for cx in (CX_L, CX_R) for a, b, c, d in aromatic_lines(cx))
        parts.append(f'<g stroke-width="{f(width * 0.82)}">{seg}</g>')

    if sticks:
        seg = []
        for cx, pm, want in ((CX_L, LEFT_POS, left), (CX_R, RIGHT_POS, right)):
            for p in want:
                (x1, y1, x2, y2), _ = chlorine(cx, pm[p])
                seg.append(f'<line x1="{f(x1)}" y1="{f(y1)}" x2="{f(x2)}" y2="{f(y2)}"/>')
        parts.append(f'<g stroke-width="{f(width * 0.88)}">' + "".join(seg) + "</g>")

    if labels:
        seg = []
        for cx, pm, want in ((CX_L, LEFT_POS, left), (CX_R, RIGHT_POS, right)):
            for p in want:
                _, (lx, ly, anchor) = chlorine(cx, pm[p])
                seg.append(f'<text x="{f(lx)}" y="{f(ly)}" text-anchor="{anchor}" '
                           f'dominant-baseline="central" font-size="{f(label_size)}">Cl</text>')
        parts.append(f'<g stroke="none" fill="{stroke}" font-family="IBM Plex Sans, '
                     f'Segoe UI, sans-serif" font-weight="500">' + "".join(seg) + "</g>")

    return (f'<g fill="none" stroke="{stroke}" stroke-width="{f(width)}" '
            f'stroke-linecap="round" stroke-linejoin="miter">' + "".join(parts) + "</g>")


def flat_bounds(*, sticks: bool, labels: bool, pad: float = 4.0):
    xs, ys = [], []
    for cx in (CX_L, CX_R):
        for x, y in vertices(cx):
            xs.append(x)
            ys.append(y)
    if sticks or labels:
        for cx, pm, want in ((CX_L, LEFT_POS, PCB77_LEFT), (CX_R, RIGHT_POS, PCB77_RIGHT)):
            for p in want:
                (_, _, sx, sy), (lx, ly, anchor) = chlorine(cx, pm[p])
                xs += [sx]
                ys += [sy]
                if labels:
                    xs += [lx - (13 if anchor in ("end", "middle") else 0),
                           lx + (13 if anchor in ("start", "middle") else 0)]
                    ys += [ly - 7, ly + 7]
    return (min(xs) - pad, min(ys) - pad,
            max(xs) - min(xs) + 2 * pad, max(ys) - min(ys) + 2 * pad)


def flat_svg(*, stroke=TEAL, width=3.4, aromatic=True, sticks=True, labels=True) -> str:
    x, y, w, h = flat_bounds(sticks=sticks, labels=labels)
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{f(x)} {f(y)} {f(w)} {f(h)}" '
            f'width="{f(w)}" height="{f(h)}" role="img" aria-label="PCB 77">'
            + molecule(stroke=stroke, width=width, aromatic=aromatic, sticks=sticks,
                       labels=labels) + "</svg>")


def tile_svg(*, stroke=TEAL, width: float, aromatic: bool, sticks: bool,
             background: str | None = None, radius_frac: float = 0.0) -> str:
    """A square icon. The viewBox is fitted after rotation, plus half the stroke width."""
    pts: list[tuple[float, float]] = []
    for cx in (CX_L, CX_R):
        pts.extend(vertices(cx))
    pts += [(CX_L + R, CY), (CX_R - R, CY)]
    if sticks:
        for cx, pm, want in ((CX_L, LEFT_POS, PCB77_LEFT), (CX_R, RIGHT_POS, PCB77_RIGHT)):
            for p in want:
                (_, _, sx, sy), _ = chlorine(cx, pm[p])
                pts.append((sx, sy))

    cx0 = (min(p[0] for p in pts) + max(p[0] for p in pts)) / 2
    cy0 = (min(p[1] for p in pts) + max(p[1] for p in pts)) / 2
    a = math.radians(ICON_ROTATION)
    ca, sa = math.cos(a), math.sin(a)
    rot = [((x - cx0) * ca - (y - cy0) * sa, (x - cx0) * sa + (y - cy0) * ca) for x, y in pts]
    half = max(max(abs(x) for x, _ in rot), max(abs(y) for _, y in rot)) + width / 2 + 2.0

    # A filled tile needs more room around the molecule than a transparent one.
    if background:
        half *= 1.22
    side = half * 2
    x0, y0 = cx0 - half, cy0 - half

    body = molecule(stroke=stroke, width=width, aromatic=aromatic, sticks=sticks, labels=False)
    bg = ""
    if background:
        r = side * radius_frac
        bg = (f'<rect x="{f(x0)}" y="{f(y0)}" width="{f(side)}" height="{f(side)}" '
              f'rx="{f(r)}" ry="{f(r)}" fill="{background}"/>')
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{f(x0)} {f(y0)} {f(side)} '
            f'{f(side)}" width="{f(side)}" height="{f(side)}" role="img" '
            f'aria-label="Quant 209">{bg}'
            f'<g transform="rotate({f(ICON_ROTATION)} {f(cx0)} {f(cy0)})">{body}</g></svg>')


#: What the mark carries at each icon size.
BANDS = {
    256: dict(width=3.4, aromatic=True, sticks=True),
    128: dict(width=3.8, aromatic=True, sticks=True),
    64: dict(width=4.4, aromatic=True, sticks=True),
    48: dict(width=5.6, aromatic=True, sticks=False),
    32: dict(width=7.0, aromatic=True, sticks=False),
    24: dict(width=9.0, aromatic=False, sticks=False),
    16: dict(width=12.0, aromatic=False, sticks=False),
}


# ----------------------------------------------------------------------------------------
# The lockup, with the wordmark converted to outlines so no font has to be installed
# ----------------------------------------------------------------------------------------

def path_to_svg(path) -> str:
    """Convert a ``QPainterPath`` to SVG path data."""
    from PySide6.QtGui import QPainterPath

    out: list[str] = []
    i = 0
    n = path.elementCount()
    while i < n:
        e = path.elementAt(i)
        if e.type == QPainterPath.ElementType.MoveToElement:
            out.append(f"M{e.x:.2f},{e.y:.2f}")
            i += 1
        elif e.type == QPainterPath.ElementType.LineToElement:
            out.append(f"L{e.x:.2f},{e.y:.2f}")
            i += 1
        elif e.type == QPainterPath.ElementType.CurveToElement:
            c1, c2, to = e, path.elementAt(i + 1), path.elementAt(i + 2)
            out.append(f"C{c1.x:.2f},{c1.y:.2f} {c2.x:.2f},{c2.y:.2f} {to.x:.2f},{to.y:.2f}")
            i += 3
        else:
            i += 1
    return "".join(out) + "Z"


def wordmark_paths(app_font: str = "IBM Plex Sans"):
    """Outlines for QUANT and 209, plus the width each occupies."""
    from PySide6.QtGui import QFont, QPainterPath

    def build(text: str, px: int, weight, spacing: float):
        font = QFont(app_font)
        font.setPixelSize(px)
        font.setWeight(weight)
        font.setLetterSpacing(QFont.SpacingType.AbsoluteSpacing, spacing)
        p = QPainterPath()
        p.addText(0.0, 0.0, font, text)
        return p

    quant = build("QUANT", 26, QFont.Weight.Normal, 5.2)
    num = build("209", 34, QFont.Weight.DemiBold, 3.4)
    return quant, num


def lockup_svg(*, dark: bool = False) -> str:
    stroke = TEAL_DARK if dark else TEAL
    ink = "#F2F5F5" if dark else INK
    rule = "#2E3A3A" if dark else "#DCE4E4"

    quant, num = wordmark_paths()
    qw = quant.boundingRect().width()
    nw = num.boundingRect().width()

    mx, my, mw, mh = flat_bounds(sticks=True, labels=True)
    scale = 1.22
    height = 104.0
    left_pad = 8.0
    tx = left_pad - mx * scale
    ty = (height - mh * scale) / 2 - my * scale
    div = left_pad + mw * scale + 24.0
    text_x = div + 22.0
    width = text_x + max(qw, nw) + 10.0

    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {f(width)} {f(height)}" \
width="{f(width)}" height="{f(height)}" role="img" aria-label="Quant 209">
<g transform="translate({f(tx)} {f(ty)}) scale({f(scale)})">\
{molecule(stroke=stroke, width=3.4 / scale, aromatic=True, sticks=True, labels=True)}</g>
<line x1="{f(div)}" y1="26" x2="{f(div)}" y2="78" stroke="{rule}" stroke-width="1.4"/>
<g transform="translate({f(text_x)} 47)" fill="{ink}"><path d="{path_to_svg(quant)}"/></g>
<g transform="translate({f(text_x)} 83)" fill="{stroke}"><path d="{path_to_svg(num)}"/></g>
</svg>'''


# ----------------------------------------------------------------------------------------
# Rendering
# ----------------------------------------------------------------------------------------

def render_png(svg: str, px_w: int, px_h: int, background: str | None = None):
    from PySide6.QtCore import QRectF
    from PySide6.QtGui import QColor, QImage, QPainter
    from PySide6.QtSvg import QSvgRenderer

    img = QImage(px_w, px_h, QImage.Format.Format_ARGB32)
    img.fill(QColor(background) if background else QColor(0, 0, 0, 0))
    p = QPainter(img)
    p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    QSvgRenderer(svg.encode("utf-8")).render(p, QRectF(0, 0, px_w, px_h))
    p.end()
    return img


def qimage_to_pil(img):
    from PIL import Image

    buf = io.BytesIO()
    img.save_bytes = None
    from PySide6.QtCore import QBuffer, QByteArray

    ba = QByteArray()
    qb = QBuffer(ba)
    qb.open(QBuffer.OpenModeFlag.WriteOnly)
    img.save(qb, "PNG")
    qb.close()
    buf.write(bytes(ba))
    buf.seek(0)
    return Image.open(buf).convert("RGBA")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Write the Quant 209 logo assets.")
    ap.add_argument("--assets", default=str(ROOT / "quant209" / "app" / "assets"))
    ap.add_argument("--images", default=str(ROOT / "docs" / "images"))
    args = ap.parse_args(argv)

    from PySide6.QtWidgets import QApplication

    QApplication.instance() or QApplication([])

    assets = Path(args.assets)
    images = Path(args.images)
    assets.mkdir(parents=True, exist_ok=True)
    images.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    def write(path: Path, text: str) -> None:
        path.write_text(text, encoding="utf-8", newline="\n")
        written.append(path)

    # The mark lying flat, for the sidebar.
    write(assets / "mark.svg", flat_svg())

    # The turned mark on its own, for use inside the interface.
    write(assets / "icon.svg", tile_svg(width=4.4, aromatic=True, sticks=True))

    # The application icon: the mark reversed out of a teal tile.
    write(assets / "icon-app.svg",
          tile_svg(stroke=PAPER, width=4.4, aromatic=True, sticks=True,
                   background=TEAL, radius_frac=0.16))

    # The Windows icon, one image per size so each carries the right amount of detail.
    from PIL import Image

    frames = []
    for size in sorted(BANDS, reverse=True):
        band = BANDS[size]
        svg = tile_svg(stroke=PAPER, background=TEAL, radius_frac=0.16, **band)
        frames.append(qimage_to_pil(render_png(svg, size, size)))
    ico = assets / "quant209.ico"
    frames[0].save(ico, format="ICO",
                   sizes=[(im.width, im.height) for im in frames],
                   append_images=frames[1:])
    written.append(ico)

    # The lockup, for the documentation.
    write(images / "logo.svg", lockup_svg())
    write(images / "logo-dark.svg", lockup_svg(dark=True))
    for name, dark in (("logo.png", False), ("logo-dark.png", True)):
        svg = lockup_svg(dark=dark)
        w = float(svg.split('viewBox="0 0 ')[1].split(" ")[0])
        scale = 2.0
        img = render_png(svg, int(w * scale), int(104 * scale))
        img.save(str(images / name))
        written.append(images / name)

    for p in written:
        print(f"  {p.relative_to(ROOT)}  ({p.stat().st_size:,} bytes)")
    print(f"\n{len(written)} files written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
