"""Every chart in the application, and the only module that imports pyqtgraph.

Why this module exists at all
------------------------------
Two reasons, and the second is the one that matters more.

**Cold start.** pyqtgraph costs ~400 ms to import (it pulls numpy and walks the Qt binding
detection on the way in). The application's budget is three seconds from double-click to a
usable window, so nothing on the startup path may touch it. Keeping *every* ``import
pyqtgraph`` in one module makes that rule checkable rather than aspirational: the guard test
asserts this module is absent from ``sys.modules`` after importing the app, and the only way
to break it is to import :mod:`quant209.app.charts` from somewhere on the startup path.

**Two charts drawn by two pages drift.** The Review page draws a chromatogram in four
places: the sample browser, the calibration preview, the recovery preview, and whatever
comes next. If the shading, the baseline pen and the apex marker were re-specified at
each call site, the fourth one would look like a different application. So the whole recipe
lives in :class:`ChromatogramView` and the pages only say *which* trace and *which* peaks.

The overlay recipe, and what each mark means scientifically
------------------------------------------------------------
The chromatogram browser is the audit centrepiece: it is where an analyst goes to ask "is
that number a peak?". So every mark on it is a persisted fact from :class:`~quant209.core.
peaks.PeakResult`, never a re-derivation:

``the trace``
    the acquired transition, thin and near-black. One pixel wide and dark so the shading
    reads as *under* the data rather than as more data.
``the shaded region``
    between the persisted straight baseline (``baseline_start`` at ``i_start`` to
    ``baseline_end`` at ``i_end``) and the trace, clipped to the baseline exactly as
    ``core.peaks.integrate_region`` clips it. **This is the integrated area, drawn.** If the
    shading does not look like the peak, the area is wrong. That is the single question
    this screen exists to answer.
``its hue``
    accent teal for an ordinary chord baseline, violet for a peak whose area was measured
    above a *cluster* baseline (the ``clustered`` flag). Not a warning colour, because
    clustering is not a defect. It says which baseline the two persisted numbers describe,
    and anyone recomputing the area by hand needs to know.
``a paler fill and a hollow apex``
    the presence test rejected this analyte (``absent``). Its region is a bounded estimate of
    whatever occupies the analyte's retention-time territory, not a peak, so it must not be
    painted with the same confidence as one. Its apex may legitimately sit outside the
    shaded region, which is exactly what the hollow marker lets you see.
``a hairline at each boundary``
    the perpendicular drop at ``i_start`` and ``i_end``, from the baseline up to the trace.
    Without it a run of five clustered congeners is one continuous violet blob and the
    screen cannot answer "where was 102 split from 98?" That is the question a cluster
    baseline exists to raise.

Labels thin out by prominence, not by position
------------------------------------------------
A crowded calibration transition carries 37 analytes across four minutes, so most of their
names cannot be drawn. **Which ones give way is the whole design.** A left-to-right pass
keeps whichever label comes first, which on real data means a near-zero ``absent`` estimate
suppresses the name of the dominant peak two seconds later. That is the exact opposite of what an
analyst needs. So candidates are ranked *detected first, then tallest first*, placed
greedily, and a label is drawn only where its own measured text box clears every label
already placed by :data:`LABEL_MIN_GAP_PX` pixels. Measured, not assumed: ``86+97+109+119``
is twice the width of ``84``, and a constant gap that fits the second guarantees the first
overlaps its neighbour.

The test is redone from the viewport on every range change, so zooming in reveals more names
continuously rather than snapping at a threshold. ``QGraphicsItem.setVisible`` is cheap.
Rebuilding the items would not be, and this runs on every wheel event.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import QEasingCurve, QVariantAnimation, Qt
from PySide6.QtGui import QColor, QFont, QFontMetricsF
from PySide6.QtWidgets import QApplication, QWidget

from quant209.app import theme

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------

#: The acquired trace: 1 px, near-black. See the module docstring.
TRACE_WIDTH = 1

#: Alpha of the integrated-region fill, out of 255. 25%: enough to read as a filled area
#: against white, light enough that the trace stays the darkest thing in the picture.
FILL_ALPHA = 64
#: The same fill for an analyte the presence test rejected: present, but visibly not a claim.
ABSENT_FILL_ALPHA = 26
#: Alpha of the perpendicular-drop hairlines at each region's boundaries. Strong enough to
#: separate two abutting regions of the same hue, light enough not to compete with the trace.
BOUNDARY_ALPHA = 120

#: The two integration-baseline hues. Violet is ``theme.CHART_SERIES``'s fifth entry, chosen
#: over the palette's amber/red because ``clustered`` is a statement about which baseline was
#: drawn, not a warning about the number.
CLUSTER_COLOR = theme.CHART_SERIES[4]

#: Padding either side of a peak when zooming to it, in minutes. Wide enough to show the
#: peak's own shoulders and its nearest neighbour, which is what "is that boundary right?"
#: is actually asked against.
ZOOM_PAD_MIN = 0.2
#: Headroom above the tallest point in a zoomed window, as a fraction.
ZOOM_HEAD_FRAC = 0.15
#: How long the zoom takes. A jump is indistinguishable from a redraw. A slide shows the
#: analyst *where* in the chromatogram they landed.
ZOOM_DURATION_MS = 220

#: Clear space required between the *text boxes* of two visible analyte labels, in device
#: pixels. A gap between boxes rather than between apexes: see the module docstring on why a
#: constant centre-to-centre distance cannot work when one label is ``84`` and the next is
#: ``86+97+109+119``.
LABEL_MIN_GAP_PX = 12

#: Axis tick/label size. The floor is 10 pt: chart type smaller than the body text
#: is the most common way a scientific plot ends up unreadable on a lab monitor at 125%.
CHART_FONT_PT = 10.0

#: The surrogate-recovery acceptance band, as fractions. Mirrors
#: ``core.qaqc.QaqcParams.recovery_low``/``recovery_high``. Pinned to them by a test.
RECOVERY_LOW = 0.50
RECOVERY_HIGH = 1.50

_CONFIGURED = False


# --------------------------------------------------------------------------------------
# Colour helpers
# --------------------------------------------------------------------------------------


def blend(color: str | QColor, onto: str | QColor, weight: float) -> QColor:
    """``weight`` parts ``color`` over ``onto``, in sRGB. Used to derive tints from PALETTE.

    Deriving the chart tints rather than writing new hex values keeps a single source of
    colour truth: re-tune ``PALETTE["accent"]`` and the recovery heatmap follows it.
    """
    a, b = QColor(color), QColor(onto)
    weight = min(max(weight, 0.0), 1.0)
    return QColor(
        round(a.red() * weight + b.red() * (1 - weight)),
        round(a.green() * weight + b.green() * (1 - weight)),
        round(a.blue() * weight + b.blue() * (1 - weight)),
    )


def with_alpha(color: str | QColor, alpha: int) -> QColor:
    out = QColor(color)
    out.setAlpha(alpha)
    return out


#: The recovery heatmap's diverging ramp, as ``(fraction, colour)`` stops. Centred on 100%
#: rather than on the middle of the data: 100% is what the reference injection *defines*, so
#: it is the only meaningful centre, and a scale centred on the run's own mean would make a
#: uniformly bad run look fine.
#:
#: Palette-derived, not rainbow (see ``docs`` dataviz conventions): amber below the band,
#: a pale green at the centre, teal above. The two extremes are the palette hues blended
#: halfway onto white, which is dark enough to read as "far out" and still light enough for
#: the palette's ordinary dark text. One text colour for every cell keeps the numbers in the
#: grid crisp, and the eye reads the *column*, not the contrast.
_RECOVERY_STOPS: tuple[tuple[float, QColor], ...] = ()


def _recovery_stops() -> tuple[tuple[float, QColor], ...]:
    global _RECOVERY_STOPS
    if not _RECOVERY_STOPS:
        surface = theme.PALETTE["surface"]
        _RECOVERY_STOPS = (
            (0.00, blend(theme.PALETTE["warning"], surface, 0.55)),
            (RECOVERY_LOW, QColor(theme.PALETTE["warning_soft"])),
            (1.00, QColor(theme.PALETTE["success_soft"])),
            (RECOVERY_HIGH, QColor(theme.PALETTE["accent_soft"])),
            (2.00, blend(theme.PALETTE["accent"], surface, 0.45)),
        )
    return _RECOVERY_STOPS


def recovery_color(fraction: float) -> QColor:
    """The heatmap fill for one surrogate recovery (``1.0`` == 100%).

    Piecewise-linear over :data:`_RECOVERY_STOPS`, clamped at both ends. A 340% recovery
    and a 200% one are both "far too high", and letting the ramp run on would spend the whole
    scale on one broken cell.
    """
    stops = _recovery_stops()
    if not np.isfinite(fraction):
        return QColor(theme.PALETTE["surface"])
    if fraction <= stops[0][0]:
        return QColor(stops[0][1])
    for (x0, c0), (x1, c1) in zip(stops, stops[1:]):
        if fraction <= x1:
            t = (fraction - x0) / (x1 - x0) if x1 > x0 else 0.0
            return blend(c1, c0, t)
    return QColor(stops[-1][1])


def recovery_in_band(fraction: float) -> bool:
    """Whether a recovery is inside the 50-150% acceptance band."""
    return bool(np.isfinite(fraction)) and RECOVERY_LOW <= fraction <= RECOVERY_HIGH


# --------------------------------------------------------------------------------------
# Labels
# --------------------------------------------------------------------------------------


@dataclass
class LabelMark:
    """One analyte name over its apex, and everything the thinning pass needs to rank it.

    ``half_px`` is measured once, at construction, from the label's own text: the font never
    changes afterwards and re-measuring 37 strings on every wheel event would be the one
    expensive thing in a repaint that has to stay under a frame.
    """

    x: float
    item: pg.TextItem
    #: The analyte this names, so :meth:`ChromatogramView.pin_label` can find it.
    text: str
    #: The presence test's verdict. Detected analytes are named before rejected ones, whose
    #: "peak" is a bounded estimate rather than a measurement.
    detected: bool
    #: Baseline-corrected apex height: the peak's own persisted number, which is the
    #: prominence the eye ranks the picture by. Not the raw trace value at the apex: on a
    #: transition riding a large baseline that would be almost the same for every peak.
    height: float
    #: Half the rendered text width, in device pixels.
    half_px: float

    @property
    def rank(self) -> tuple[int, float, float]:
        """Sort key: detected first, then tallest, then earliest. Lower sorts first."""
        return (0 if self.detected else 1, -self.height, self.x)


# --------------------------------------------------------------------------------------
# pyqtgraph configuration
# --------------------------------------------------------------------------------------


def chart_font(point_size: float = CHART_FONT_PT) -> QFont:
    """The platform UI font at chart size: never a named family (see :mod:`theme`)."""
    app = QApplication.instance()
    font = QFont(app.font()) if app is not None else QFont()
    font.setPointSizeF(point_size)
    return font


def configure_pyqtgraph() -> None:
    """Set pyqtgraph's global options to the design system. Idempotent. Called on import.

    Global rather than per-widget because several of these (``background``, ``foreground``)
    are read by ``PlotWidget.__init__`` and by axis items at construction, so a widget
    created before they are set keeps both Qt's and pyqtgraph's own defaults: a black
    background with white axes, which is precisely the "half-styled" look the theme module
    exists to prevent.
    """
    global _CONFIGURED
    if _CONFIGURED:
        return
    pg.setConfigOptions(
        background=theme.PALETTE["surface"],
        foreground=theme.PALETTE["text_secondary"],
        antialias=True,
        useOpenGL=False,
        imageAxisOrder="row-major",
    )
    _CONFIGURED = True


def apply_chart_theme(plot_widget: pg.PlotWidget) -> None:
    """Style one plot: white ground, hairline grid, palette axes, ≥10 pt type, no chrome.

    The "no chrome" half is not cosmetic. pyqtgraph ships a right-click menu of thirty
    developer-facing options (log scaling, FFT, alpha, downsample mode) and a small ``A``
    button over the corner of every plot. Both are debugging affordances. Neither answers a
    question an analyst has.
    """
    plot_item = plot_widget.getPlotItem()
    plot_widget.setBackground(theme.PALETTE["surface"])
    plot_widget.setMenuEnabled(False)
    plot_item.hideButtons()
    plot_item.setMenuEnabled(False)
    plot_item.showGrid(x=True, y=True, alpha=0.12)

    font = chart_font()
    axis_pen = pg.mkPen(theme.PALETTE["border_strong"], width=1)
    text_color = QColor(theme.PALETTE["text_secondary"])
    for name in ("bottom", "left"):
        axis = plot_item.getAxis(name)
        axis.setPen(axis_pen)
        axis.setTextPen(text_color)
        axis.setStyle(tickFont=font, tickTextOffset=6)
    for name in ("top", "right"):
        # Shown but empty: the frame closes the plotting area on all four sides, which is
        # what stops a tall peak from looking as though it runs off the top of the picture.
        axis = plot_item.getAxis(name)
        axis.setPen(axis_pen)
        axis.setStyle(showValues=False)
        axis.show()

    label_style = {"color": theme.PALETTE["text_secondary"], "font-size": f"{CHART_FONT_PT:g}pt"}
    plot_item.setLabel("bottom", "Retention time (min)", **label_style)
    # ``units`` rather than a parenthesised word, because it turns on pyqtgraph's SI prefix:
    # a calibration transition peaks near 100 000 counts, and six-digit tick labels cost 40 px
    # of plot width apiece and are read digit-by-digit. "100 k counts" is one glance.
    plot_item.setLabel("left", "Intensity", units="counts", **label_style)
    plot_item.vb.setDefaultPadding(0.045)
    # The two performance switches used here: pyqtgraph decides the stride from the
    # viewport ("peak" keeps the extremes, so a spike is never downsampled away) and draws
    # only the samples inside the current x range. A 30 000-point transition then costs a
    # viewport's worth of line segments per repaint rather than 30 000.
    plot_item.setDownsampling(auto=True, mode="peak")
    plot_item.setClipToView(True)


# --------------------------------------------------------------------------------------
# The chromatogram view
# --------------------------------------------------------------------------------------


class ChromatogramView(pg.PlotWidget):
    """One transition's trace with every analyte's integration drawn over it.

    The page hands over a trace and the peaks measured on it. Everything else is this
    class's: the shading, the baselines, the apex markers, the label thinning, the zoom.
    Reused at three sizes: the sample browser (full), the calibration preview and
    the recovery preview (``compact=True``, which drops the axis titles that would eat half
    the height of a 180 px plot).
    """

    def __init__(self, parent: QWidget | None = None, compact: bool = False) -> None:
        configure_pyqtgraph()
        super().__init__(parent=parent, background=theme.PALETTE["surface"])
        self.compact = compact
        apply_chart_theme(self)
        if compact:
            # A left-axis title is rotated, so it spends *width*, about fourteen pixels, not
            # the height a 200 px preview cannot spare. That is the bottom axis's problem, and
            # the bottom axis is the one shortened to "min".
            self.getPlotItem().setLabel(
                "left", "Intensity", units="counts",
                color=theme.PALETTE["text_secondary"], **{"font-size": "9pt"},
            )
            self.getPlotItem().setLabel(
                "bottom", "min", color=theme.PALETTE["text_secondary"], **{"font-size": "9pt"}
            )
            self.getPlotItem().getAxis("left").setWidth(56)

        self.plot_item = self.getPlotItem()
        self.view_box = self.plot_item.vb

        #: The acquired trace. ``None`` until :meth:`show_trace` is called.
        self.trace = None
        #: One :class:`pyqtgraph.FillBetweenItem` per drawn peak: the integrated regions.
        self.region_items: list[pg.FillBetweenItem] = []
        #: The straight baselines, one per region, drawn so the shading's lower edge is a
        #: visible line rather than an implied one.
        self.baseline_items: list[pg.PlotCurveItem] = []
        #: All apexes in one scatter (N spots): one graphics item repaints far faster than N.
        self.apex_scatter: pg.ScatterPlotItem | None = None
        #: The perpendicular drops at every region boundary, one item per hue.
        self.boundary_items: list[pg.PlotCurveItem] = []
        #: One :class:`LabelMark` per peak, sorted by x. See the module docstring on thinning.
        self.labels: list[LabelMark] = []
        #: The analyte whose label is drawn whatever the crowding: the one the analyst
        #: asked for. See :meth:`pin_label`.
        self.pinned: str | None = None
        #: The peaks currently drawn, in the order they were handed over.
        self.peaks: list = []
        #: The x/y range the last :meth:`zoom_to_peak` is heading for. Read by the tests, so
        #: they never have to wait for an animation to assert that the right thing happened.
        self.zoom_target: tuple[tuple[float, float], tuple[float, float]] | None = None

        self._times: np.ndarray = np.empty(0)
        self._intensities: np.ndarray = np.empty(0)
        self._animation: QVariantAnimation | None = None
        self.view_box.sigRangeChanged.connect(self._refresh_labels)

    # ---- drawing ------------------------------------------------------------------------

    def clear_plot(self) -> None:
        """Drop every item. Called before each redraw and when a page has nothing to show."""
        self._stop_animation()
        self.plot_item.clear()
        self.trace = None
        self.region_items = []
        self.baseline_items = []
        self.apex_scatter = None
        self.boundary_items = []
        self.labels = []
        self.pinned = None
        self.peaks = []
        self.zoom_target = None
        self._times = np.empty(0)
        self._intensities = np.empty(0)

    def show_trace(
        self,
        times: Sequence[float] | np.ndarray,
        intensities: Sequence[float] | np.ndarray,
        peaks: Sequence = (),
        *,
        labelled: bool = True,
    ) -> None:
        """Draw one transition and overlay ``peaks``' integrations on it.

        ``peaks`` must be the peaks measured *on this trace*. The caller selects them by
        transition. A peak whose recorded indices do not address this array (a re-plot
        against a different conversion of the file, say) is skipped rather than drawn
        somewhere plausible and wrong.
        """
        self.clear_plot()
        self._times = np.asarray(times, dtype=float)
        self._intensities = np.asarray(intensities, dtype=float)
        if self._times.size == 0:
            return

        self.trace = self.plot_item.plot(
            self._times,
            self._intensities,
            pen=pg.mkPen(theme.PALETTE["text"], width=TRACE_WIDTH),
            name="trace",
        )

        spots: list[dict] = []
        drops: dict[str, list[list[float]]] = {}
        for peak in peaks:
            if not self._draw_region(peak, spots, drops, labelled=labelled):
                continue
            self.peaks.append(peak)
        if spots:
            self.apex_scatter = pg.ScatterPlotItem(spots=spots, pxMode=True, size=7)
            self.plot_item.addItem(self.apex_scatter)
        self._draw_boundaries(drops)
        self.labels.sort(key=lambda mark: mark.x)
        self.view_box.enableAutoRange()
        self.view_box.autoRange()
        self._refresh_labels()

    def _draw_boundaries(self, drops: dict[str, list[list[float]]]) -> None:
        """The perpendicular drops, one ``connect="pairs"`` curve per hue.

        One item per hue rather than two per peak: a 37-analyte calibration transition would
        otherwise add 74 graphics items to a plot that is re-ranged on every wheel event.
        """
        for hue, segments in drops.items():
            if not segments:
                continue
            xs = np.array([value for segment in segments for value in segment[:2]], dtype=float)
            ys = np.array([value for segment in segments for value in segment[2:]], dtype=float)
            curve = pg.PlotCurveItem(
                xs,
                ys,
                connect="pairs",
                pen=pg.mkPen(with_alpha(hue, BOUNDARY_ALPHA), width=1),
            )
            self.plot_item.addItem(curve)
            self.boundary_items.append(curve)

    def _draw_region(
        self, peak, spots: list[dict], drops: dict[str, list[list[float]]], *, labelled: bool
    ) -> bool:
        """One peak's shaded region, baseline, apex spot and label. False if unplottable."""
        i_start, i_end = int(getattr(peak, "i_start", -1)), int(getattr(peak, "i_end", -1))
        n = self._times.size
        if not (0 <= i_start < i_end < n):
            return False

        x = self._times[i_start : i_end + 1]
        y = self._intensities[i_start : i_end + 1]
        baseline = np.linspace(
            float(peak.baseline_start), float(peak.baseline_end), x.size
        )
        # Clipped to the baseline exactly as ``core.peaks.integrate_region`` clips it, so the
        # shaded area *is* the reported area rather than a picture that resembles it.
        top = np.maximum(y, baseline)

        flags = list(getattr(peak, "flags", ()) or ())
        absent = "absent" in flags
        clustered = "clustered" in flags
        hue = CLUSTER_COLOR if clustered else theme.PALETTE["accent"]
        alpha = ABSENT_FILL_ALPHA if absent else FILL_ALPHA

        base_curve = pg.PlotCurveItem(
            x,
            baseline,
            pen=pg.mkPen(
                with_alpha(hue, 150),
                width=1,
                style=Qt.PenStyle.DashLine if absent else Qt.PenStyle.SolidLine,
            ),
        )
        top_curve = pg.PlotCurveItem(x, top, pen=pg.mkPen(None))
        fill = pg.FillBetweenItem(top_curve, base_curve, brush=pg.mkBrush(with_alpha(hue, alpha)))
        for item in (base_curve, top_curve, fill):
            self.plot_item.addItem(item)
        self.region_items.append(fill)
        self.baseline_items.append(base_curve)

        # The two perpendicular drops, collected rather than drawn: see _draw_boundaries.
        drops.setdefault(hue, []).extend(
            [
                [float(x[0]), float(x[0]), float(baseline[0]), float(top[0])],
                [float(x[-1]), float(x[-1]), float(baseline[-1]), float(top[-1])],
            ]
        )

        i_apex = int(getattr(peak, "i_apex", -1))
        if 0 <= i_apex < n:
            apex_x, apex_y = float(self._times[i_apex]), float(self._intensities[i_apex])
        else:  # an apex outside the trace cannot be marked. The region still can be drawn
            apex_x, apex_y = float(x[int(np.argmax(top))]), float(top.max())
        spots.append(
            {
                "pos": (apex_x, apex_y),
                "symbol": "o",
                "pen": pg.mkPen(hue, width=1.4),
                "brush": pg.mkBrush(None) if absent else pg.mkBrush(hue),
            }
        )

        if labelled:
            text = str(peak.label)
            font = chart_font(9)
            label = pg.TextItem(
                text,
                color=theme.PALETTE["text_secondary"] if absent else theme.PALETTE["text"],
                anchor=(0.5, 1.25),
                # A near-opaque plate of the plot's own ground. A label that lands on a
                # neighbouring peak's trace is otherwise read letter-by-letter against a
                # black line, which is the difference between a name and a smudge.
                fill=pg.mkBrush(with_alpha(theme.PALETTE["surface"], 205)),
            )
            label.setFont(font)
            label.setPos(apex_x, apex_y)
            self.plot_item.addItem(label)
            self.labels.append(
                LabelMark(
                    x=apex_x,
                    item=label,
                    text=text,
                    detected=not absent,
                    height=abs(float(getattr(peak, "height", 0.0) or 0.0)) or abs(apex_y),
                    half_px=QFontMetricsF(font).horizontalAdvance(text) / 2.0,
                )
            )
        return True

    # ---- label thinning -------------------------------------------------------------------

    def _refresh_labels(self, *_args) -> None:
        """Draw the most prominent labels that still fit. See the module docstring.

        Greedy over :attr:`LabelMark.rank`: the tallest detected analyte in the window is
        placed first and never gives way, and each later candidate is drawn only where its
        own text box clears everything already placed.
        """
        if not self.labels:
            return
        (x0, x1), _y = self.view_box.viewRange()
        width_px = max(float(self.view_box.width()), 1.0)
        if x1 <= x0:
            return
        per_unit = width_px / (x1 - x0)
        placed: list[tuple[float, float]] = []
        shown: set[int] = set()
        visible = [m for m in self.labels if x0 <= m.x <= x1]
        # The pinned analyte first, ahead of every ranking rule: the analyst clicked its row
        # and the view flew to it, and a picture that then names only its taller neighbour is
        # answering a question nobody asked.
        order = sorted(visible, key=lambda m: (m.text != self.pinned, *m.rank))
        for mark in order:
            centre = (mark.x - x0) * per_unit
            left, right = centre - mark.half_px, centre + mark.half_px
            if any(
                left - LABEL_MIN_GAP_PX < taken_right and right + LABEL_MIN_GAP_PX > taken_left
                for taken_left, taken_right in placed
            ):
                continue
            placed.append((left, right))
            shown.add(id(mark.item))
        for mark in self.labels:
            mark.item.setVisible(id(mark.item) in shown)

    def pin_label(self, text: str | None) -> None:
        """Always draw ``text``'s label, whatever else wants the space. ``None`` clears it."""
        self.pinned = text
        self._refresh_labels()

    def visible_label_count(self) -> int:
        """How many analyte labels are actually painted. The clutter test reads this."""
        return sum(1 for mark in self.labels if mark.item.isVisible())

    # ---- zoom ------------------------------------------------------------------------------

    def peak_range(self, peak) -> tuple[tuple[float, float], tuple[float, float]]:
        """The ``(xRange, yRange)`` that frames one peak, with its shoulders and headroom."""
        x0 = float(peak.rt_start) - ZOOM_PAD_MIN
        x1 = float(peak.rt_end) + ZOOM_PAD_MIN
        y_hi = float(getattr(peak, "height", 0.0)) or 1.0
        y_lo = 0.0
        if self._times.size:
            window = (self._times >= x0) & (self._times <= x1)
            if window.any():
                inside = self._intensities[window]
                y_lo, y_hi = float(inside.min()), float(inside.max())
        span = (y_hi - y_lo) or abs(y_hi) or 1.0
        return (x0, x1), (y_lo - span * 0.05, y_hi + span * ZOOM_HEAD_FRAC)

    def zoom_to_peak(self, peak, *, animate: bool = True) -> None:
        """Slide the view to frame ``peak``. See :data:`ZOOM_DURATION_MS` for why it slides."""
        target = self.peak_range(peak)
        self.zoom_target = target
        self.pin_label(str(getattr(peak, "label", "")) or None)
        self._stop_animation()
        self.view_box.disableAutoRange()
        if not animate or ZOOM_DURATION_MS <= 0:
            self._apply_range(target)
            return
        start = tuple(tuple(float(v) for v in axis) for axis in self.view_box.viewRange())
        animation = QVariantAnimation(self)
        animation.setDuration(ZOOM_DURATION_MS)
        animation.setStartValue(0.0)
        animation.setEndValue(1.0)
        animation.setEasingCurve(QEasingCurve.Type.OutCubic)
        animation.valueChanged.connect(
            lambda value: self._apply_range(_lerp_range(start, target, float(value)))
        )
        # Held on the instance: a QVariantAnimation whose only reference is a local goes out
        # of scope at the end of this call and is garbage-collected mid-flight, which shows
        # up as a zoom that stops a third of the way there.
        self._animation = animation
        animation.start()

    def reset_zoom(self) -> None:
        self._stop_animation()
        self.zoom_target = None
        self.pinned = None
        self.view_box.enableAutoRange()
        self.view_box.autoRange()

    def _apply_range(self, ranges: tuple[tuple[float, float], tuple[float, float]]) -> None:
        (x0, x1), (y0, y1) = ranges
        self.view_box.setRange(xRange=(x0, x1), yRange=(y0, y1), padding=0)

    def _stop_animation(self) -> None:
        if self._animation is not None:
            self._animation.stop()
            self._animation = None


def _lerp_range(
    start: tuple[tuple[float, float], tuple[float, float]],
    end: tuple[tuple[float, float], tuple[float, float]],
    t: float,
) -> tuple[tuple[float, float], tuple[float, float]]:
    return tuple(  # type: ignore[return-value]
        tuple(a + (b - a) * t for a, b in zip(axis_start, axis_end))
        for axis_start, axis_end in zip(start, end)
    )


configure_pyqtgraph()
