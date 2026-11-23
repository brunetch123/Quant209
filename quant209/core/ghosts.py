"""Ghost peaks: a congener's in-source fragments, recorded on lower homologs' channels.

What a ghost peak is
--------------------
A PCB with n chlorines can lose HCl or Cl2 in the ion source. The fragment has exactly the
precursor mass of the transition for n-1 chlorines (HCl lost) or n-2 chlorines (Cl2 lost), and
its own loss of Cl2 gives that transition's product ion. So the instrument records a copy of
the congener's peak on the lower channel, at the congener's own retention time and with its
shape: a ghost peak. Measured on the archive, the internal standard PCB 204 (an octa-CB) leaves
4-6% of its area on the hepta channel and 14-16% on the hexa channel in every injection that
holds it (2% on hexa in 2019, a different tune), and so do the other octa-, nona- and deca-CBs
and most penta- to hepta-CBs: 2-5% one step down, up to about 25% two steps down. Only steps
whose precursor and product masses both match can happen. A tetra-CB losing Cl2 cannot give the
di-CB product ion, and mono- to tri-CBs give none (:data:`STEPS`).

A copy is not an analyte, and an analyst integrating the channel leaves it out. Integrated as it
stands it becomes one. A copy that lands on a congener's retention time is reported as that
congener: PCB 204's copy as 156+157 in every blank. A copy that lands between congeners is
swallowed by a neighbour's valley-bounded window: PCB 204's under 172 and 171+173. The
calibration standard holds every congener, so copies sit on every channel there and inflate the
windows that response factors are measured in (39, 105, 156+157, 167 and 189 among them).

What is done about it
---------------------
:func:`integrate_with_ghosts` integrates a file as acquired, then looks for each copy:

1. **Predict.** Every peak of every congener in :func:`core.config.load_ghost_fractions`, and the
   internal standard 204, highest homolog first, has its copy predicted on each lower channel
   its masses allow, as ``ratio x the run's PCB 204 fraction`` (:class:`GhostScale`). ``ratio``
   is the congener's own propensity relative to 204's, measured once on the archive. The run's
   204 fraction carries the instrument's scale, measured in the run's own blanks.
2. **Confirm.** The lower channel's trace is fitted to the parent peak's shape over the parent's
   core, allowing a sub-cycle shift. The copy counts only when the trace follows that shape
   (correlation >= :data:`SHAPE_MIN_R`). A copy the trace cannot confirm, typically one under a
   real peak of the lower homolog, is left alone: the analysts integrate that peak whole.
3. **Keep it out.**

   * A copy more than :data:`COELUTION` from every analyte's expected retention time on its
     channel is left in the trace, and a placeholder analyte is anchored on it. Neighbouring
     windows then end at the valley before it, as an analyst draws them, and no analyte's apex
     search can take it. The placeholder's own result is discarded.
   * A copy within :data:`COELUTION` of an analyte's expected retention time is subtracted from
     the trace: the measured copy, never more than predicted, so a real peak underneath keeps
     everything above the copy, and never below the local baseline. Nothing is subtracted under
     an internal standard, which is integrated whole as the analysts integrate it, nor under
     PCB 58 or PCB 159 (:data:`NO_SUBTRACTION_UNDER`), whose reported peaks usually hold the
     copy of PCB 95 and PCB 174.

4. **Bound.** The channels a copy touched are integrated again. An analyte keeps the new result
   only if it is no larger than the first and smaller by no more than :data:`LOSS_BOUND` times
   the copies inside its first window (plus :data:`LOSS_TOL_REL` and :data:`LOSS_TOL_ABS`).
   Anything else, a window or a presence verdict that moved for some other reason, keeps the
   first result. A kept new result that differs carries the ``ghost`` flag (:data:`GHOST_FLAG`),
   and every copy found is recorded on the file (:class:`Ghost`).

Measured on the 142-run validation dataset and on the benchmark the published figures come from,
see ``docs/METHODS.md`` section 2.10.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

from .config import Analyte, load_ghost_fractions
from .integrate import integrate_run
from .mzml import Chromatogram, get_transition
from .peaks import IntegrationParams, PeakResult

#: Homolog channels, fewest chlorines first.
HOMOLOGS = ("mono", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona", "deca")

#: Steps down (1 = HCl lost, 2 = Cl2 lost) whose precursor AND product masses both match the
#: lower channel. A tetra-CB losing Cl2 cannot give the di-CB product ion, and a di- or tri-CB
#: losing either cannot give the lower channel's product ion.
STEPS: dict[str, tuple[int, ...]] = {
    "deca": (1, 2), "nona": (1, 2), "octa": (1, 2), "hepta": (1, 2), "hexa": (1, 2),
    "penta": (1, 2), "tetra": (1,),
}

#: PCB 204's copy as a share of its peak, per step, used only when neither the run's blanks and
#: references nor the file itself show a copy that stands clear. Medians over the archive's
#: blanks since 2021.
DEFAULT_IS_FRACTION: dict[int, float] = {1: 0.0416, 2: 0.147}

#: A copy is confirmed when the lower trace follows the parent's shape this closely (Pearson
#: correlation over the parent's core, best sub-cycle shift).
SHAPE_MIN_R = 0.95
#: PCB 204's own copy, which scales everything else, is read off a sample only where it stands
#: this clear (correlation; signal-to-noise): a neighbour's tail can sit under it there.
IS_SHAPE_MIN_R = 0.98
IS_MIN_SN = 10.0
#: In a blank or a reference injection nothing else elutes at 204's retention time on the hepta
#: and hexa channels, so the copy needs only to be there: IS_MIN_SN, and a shape that follows
#: 204's at least this loosely. Blank 204 peaks are small and their copies noisy: in one 2024
#: run's eight blanks the hepta copy correlated 0.50-0.86 with 204's shape and the hexa copy
#: 0.86-0.96, while the fractions agreed to within 10% (hepta 3.7-4.5%, hexa 15-21%).
SCALE_MIN_R = 0.5
#: A copy this close to an analyte's expected retention time (minutes) can be taken for that
#: analyte, so it is subtracted. One farther away is fenced off.
COELUTION = 0.10
#: Congeners a copy is never subtracted under, as under an internal standard. PCB 95's HCl-loss
#: copy elutes 0.09 min before PCB 58 and PCB 174's with PCB 159, and where either congener is
#: reported the reported peak usually holds the copy. A copy farther from them is still fenced.
NO_SUBTRACTION_UNDER = frozenset({"58", "159"})
#: A copy is fenced only where it shows: its measured height at least this many times the
#: lower trace's noise.
FENCE_MIN_SN = 5.0
#: A copy is subtracted only where its height reaches this many times the lower trace's noise.
SUBTRACT_MIN_SN = 2.0
#: The corrected trace may not go below the local baseline by more than this many noise units.
FLOOR_NOISE = 0.5
#: An analyte keeps the corrected result only if it lost no more than this many times the
#: copies in its window, plus a tolerance of LOSS_TOL_REL of its area and LOSS_TOL_ABS counts*s.
LOSS_BOUND = 1.5
LOSS_TOL_REL = 0.03
LOSS_TOL_ABS = 200.0
#: The parent shape is slid across one acquisition cycle either way in this many steps.
SHIFT_STEPS = 21

#: Flag on a peak whose reported result is the one with ghost peaks kept out.
GHOST_FLAG = "ghost"
#: Label prefix of the placeholder analytes that fence copies off. Never reported.
FENCE_PREFIX = "~ghost "

FENCED = "fenced"
SUBTRACTED = "subtracted"


@dataclass(frozen=True)
class Ghost:
    """One copy found and kept out of a file's integration."""

    source: str  # the congener (or internal standard) whose fragment it is
    channel: str  # the homolog channel it was found on
    step: int  # 1: HCl lost, 2: Cl2 lost
    rt: float  # its apex, minutes
    fraction: float  # its measured size, as a share of the source peak
    predicted: float  # ratio x the run's 204 fraction
    shape_r: float  # correlation of the lower trace with the source's shape
    area: float  # its area as used, counts*s
    action: str  # FENCED or SUBTRACTED


@dataclass(frozen=True)
class GhostScale:
    """The run's PCB 204 copy, as a share of its peak, per step, and where it was measured.

    A step the run's blanks and references did not measure is read off each file's own 204
    peak where its copy stands clear, and otherwise taken from :data:`DEFAULT_IS_FRACTION`.
    """

    fractions: dict[int, float] = field(default_factory=dict)
    #: Per step: which files the median is over, e.g. ``"3 blanks"``.
    measured_in: dict[int, str] = field(default_factory=dict)

    def note(self) -> str:
        """The run note: what scale the ghost-peak correction used."""
        names = {1: "hepta", 2: "hexa"}
        parts = []
        for step in (1, 2):
            if step in self.fractions:
                parts.append(f"on the {names[step]} channel {100 * self.fractions[step]:.1f}% of its "
                             f"peak (median of {self.measured_in[step]})")
            else:
                parts.append(f"on the {names[step]} channel not measured in this run's blanks or "
                             f"references, so read off each file, or taken as "
                             f"{100 * DEFAULT_IS_FRACTION[step]:.1f}% where a file shows none")
        return ("ghost peaks: the internal standard 204's copy " + "; ".join(parts)
                + ". Every congener's copy on lower homologs' channels is predicted from it and kept "
                "out of their integrations where the trace confirms it (flag 'ghost')")


@dataclass
class GhostPlan:
    """The copies found in one file and the traces with them dealt with."""

    chroms: dict[tuple[float, float], Chromatogram]
    fences: list[Analyte]
    fence_anchors: dict[str, float]
    ghosts: list[Ghost]
    #: Transition keys whose trace changed or that carry a fence.
    keys: set[tuple[float, float]]


# --------------------------------------------------------------------------------------
# Trace helpers
# --------------------------------------------------------------------------------------


def _key(q1: float, q3: float) -> tuple[float, float]:
    return (round(q1, 1), round(q3, 1))


def _hermite(xn: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """``y(x)`` at ``xn`` by cubic Hermite interpolation (Catmull-Rom slopes), 0 outside."""
    out = np.zeros_like(xn, dtype=float)
    inside = (xn >= x[0]) & (xn <= x[-1])
    if not inside.any() or x.size < 3:
        return np.interp(xn, x, y, left=0.0, right=0.0)
    m = np.empty_like(y, dtype=float)
    m[1:-1] = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
    m[0] = (y[1] - y[0]) / (x[1] - x[0])
    m[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
    xi = xn[inside]
    i = np.clip(np.searchsorted(x, xi) - 1, 0, x.size - 2)
    dx = x[i + 1] - x[i]
    t = (xi - x[i]) / dx
    t2, t3 = t * t, t * t * t
    out[inside] = ((2 * t3 - 3 * t2 + 1) * y[i] + (t3 - 2 * t2 + t) * dx * m[i]
                   + (-2 * t3 + 3 * t2) * y[i + 1] + (t3 - t2) * dx * m[i + 1])
    return out


def _sides(t: np.ndarray, y: np.ndarray, t_lo: float, t_hi: float) -> tuple[float | None, float | None]:
    """Local baseline and noise from 0.02-0.30 min either side of ``[t_lo, t_hi]``.

    The baseline is the lower side's median, the noise the quieter side's scaled MAD, so a
    neighbour's flank on one side does not raise either.
    """
    left = y[(t >= t_lo - 0.30) & (t <= t_lo - 0.02)]
    right = y[(t >= t_hi + 0.02) & (t <= t_hi + 0.30)]
    meds = [float(np.median(a)) for a in (left, right) if a.size >= 3]
    mads = [float(1.4826 * np.median(np.abs(a - np.median(a)))) for a in (left, right) if a.size >= 5]
    return (min(meds) if meds else None), (min(mads) if mads else None)


def _peak(t: np.ndarray, y: np.ndarray, anchor: float, window: float) -> dict | None:
    """The source peak near ``anchor``: apex, extent (to 2% of its height), baseline, height.

    The search runs a little wider than the apex window but accepts only a true local maximum
    within 1.2 windows of the anchor, so a neighbour's flank rising into the search edge is not
    taken for the source. A peak under 10x the local noise is not a source worth predicting.
    """
    w = np.where((t >= anchor - 1.5 * window) & (t <= anchor + 1.5 * window))[0]
    if w.size < 3:
        return None
    i0 = int(w[np.argmax(y[w])])
    if i0 in (w[0], w[-1]) or abs(t[i0] - anchor) > 1.2 * window:
        return None
    base, noise = _sides(t, y, t[max(i0 - 1, 0)], t[min(i0 + 1, len(t) - 1)])
    if base is None:
        return None
    h = y[i0] - base
    lo = i0
    while lo > 0 and (y[lo] - base) > 0.02 * h and t[i0] - t[lo] < 0.15:
        lo -= 1
    hi = i0
    while hi < len(t) - 1 and (y[hi] - base) > 0.02 * h and t[hi] - t[i0] < 0.15:
        hi += 1
    base, noise = _sides(t, y, t[lo], t[hi])
    if base is None:
        return None
    h = float(y[i0] - base)
    if h <= 0 or (noise and h < 10 * noise):
        return None
    return {"i0": i0, "lo": lo, "hi": hi, "base": base, "h": h, "t0": float(t[i0])}


def _fit(st: np.ndarray, sy: np.ndarray, pk: dict, tt: np.ndarray, ty: np.ndarray) -> dict | None:
    """Fit the lower trace to the source peak's shape over the source's core.

    The source segment (above its own baseline) is slid across one acquisition cycle either way,
    and at each shift the lower trace (above its local baseline) is regressed on it through the
    origin, over the points where the source stands above 20% of its height. Returns the slope
    (the copy's share of the source), the correlation and the shift at the best correlation,
    with the lower trace's local noise and baseline.
    """
    t_lo, t_hi = st[pk["lo"]], st[pk["hi"]]
    if tt.size < 10 or tt[0] > t_lo - 0.35 or tt[-1] < t_hi + 0.35:
        return None
    m = (tt >= t_lo) & (tt <= t_hi)
    if m.sum() < 5:
        return None
    base_t, noise_t = _sides(tt, ty, t_lo, t_hi)
    if base_t is None:
        return None
    seg = np.zeros_like(sy, dtype=float)
    seg[pk["lo"]:pk["hi"] + 1] = sy[pk["lo"]:pk["hi"] + 1] - pk["base"]
    yt = ty[m] - base_t
    dt = float(np.median(np.diff(tt)))
    # Every shift at once: row i is the source shape slid by shifts[i] cycles.
    shifts = np.linspace(-1.0, 1.0, SHIFT_STEPS) * dt
    xs = _hermite((tt[m][None, :] - shifts[:, None]).ravel(), st, seg).reshape(shifts.size, -1)
    core = xs >= 0.2 * pk["h"]
    n = core.sum(axis=1)
    x = np.where(core, xs, 0.0)
    y = np.where(core, yt[None, :], 0.0)
    sx, sy_, sxx, syy, sxy = x.sum(1), y.sum(1), (x * x).sum(1), (y * y).sum(1), (x * y).sum(1)
    usable = (n >= 3) & (sxx > 0)
    if not usable.any():
        return None
    with np.errstate(invalid="ignore", divide="ignore"):
        vx = n * sxx - sx * sx
        vy = n * syy - sy_ * sy_
        r = np.where((vx > 0) & (vy > 0), (n * sxy - sx * sy_) / np.sqrt(np.maximum(vx * vy, 1e-300)), 0.0)
        f = np.where(sxx > 0, sxy / np.where(sxx > 0, sxx, 1.0), 0.0)
    r = np.where(usable, r, -np.inf)
    best = int(np.argmax(r))  # the first of equal maxima, the most negative shift
    return {"f": float(f[best]), "r": float(r[best]), "shift": float(shifts[best]), "noise": noise_t or 1.0,
            "base": base_t, "seg": seg}


def _channels(analytes: Sequence[Analyte]) -> tuple[dict[str, tuple[float, float]], dict[str, Analyte]]:
    """Each homolog's channel (from its natives and PCB 204), and a native to copy for fences."""
    channel: dict[str, tuple[float, float]] = {}
    template: dict[str, Analyte] = {}
    for a in analytes:
        if a.role in ("native", "IS_204") and a.homolog in HOMOLOGS:
            channel.setdefault(a.homolog, (a.q1, a.q3))
            if a.role == "native":
                template.setdefault(a.homolog, a)
    return channel, template


# --------------------------------------------------------------------------------------
# The scale: PCB 204's copy
# --------------------------------------------------------------------------------------


def is_fraction(
    chroms: Mapping[tuple[float, float], Chromatogram], analytes: Sequence[Analyte],
    anchors: Mapping[str, float], params: IntegrationParams, min_r: float = IS_SHAPE_MIN_R,
) -> dict[int, float]:
    """This injection's PCB 204 copy as a share of its peak, per step, where it stands clear.

    A step is reported only where the lower trace follows 204's shape to ``min_r`` and the copy
    reaches :data:`IS_MIN_SN` times the local noise. In a concentrated sample a neighbour's tail
    can sit under the copy, which is why a run takes the median over its blanks
    (:func:`run_scale`, read with :data:`SCALE_MIN_R`) and a file falls back on its own value
    only when the run has none.
    """
    channel, _ = _channels(analytes)
    is204 = next((a for a in analytes if a.role == "IS_204"), None)
    if is204 is None or is204.label not in anchors or is204.homolog not in channel:
        return {}
    src = get_transition(chroms, *channel[is204.homolog])
    if src is None:
        return {}
    st, sy = src.times, src.intensities.astype(float)
    pk = _peak(st, sy, anchors[is204.label], params.rt_window)
    if pk is None:
        return {}
    out: dict[int, float] = {}
    for k in STEPS.get(is204.homolog, ()):
        lower = HOMOLOGS[HOMOLOGS.index(is204.homolog) - k]
        if lower not in channel:
            continue
        tgt = get_transition(chroms, *channel[lower])
        if tgt is None:
            continue
        res = _fit(st, sy, pk, tgt.times, tgt.intensities.astype(float))
        if res and res["r"] >= min_r and res["f"] > 0 and res["f"] * pk["h"] >= IS_MIN_SN * res["noise"]:
            out[k] = res["f"]
    return out


def run_scale(blanks: Iterable[Mapping[int, float]], references: Iterable[Mapping[int, float]] = ()) -> GhostScale:
    """The run's PCB 204 copy per step: the median over its blanks, else over its references.

    Blanks hold the internal standard and next to nothing else, so its copy stands clear there.
    Reference injections are the fallback for a run without a usable blank.
    """
    blanks, references = list(blanks), list(references)
    fractions: dict[int, float] = {}
    measured_in: dict[int, str] = {}
    for step in (1, 2):
        for files, kind in ((blanks, "blank"), (references, "reference injection")):
            values = [f[step] for f in files if step in f]
            if values:
                fractions[step] = float(np.median(values))
                measured_in[step] = f"{len(values)} {kind}{'s' if len(values) != 1 else ''}"
                break
    return GhostScale(fractions=fractions, measured_in=measured_in)


# --------------------------------------------------------------------------------------
# Finding the copies
# --------------------------------------------------------------------------------------


def find_ghosts(
    chroms: Mapping[tuple[float, float], Chromatogram], analytes: Sequence[Analyte],
    anchors: Mapping[str, float], params: IntegrationParams, scale: GhostScale,
) -> GhostPlan:
    """Find every confirmed copy in one file; fence or subtract it (see the module docstring).

    Sources are walked from the highest homolog down, each on its channel as already corrected
    for the homologs above it, so a copy of a copy is never counted twice.
    """
    window = params.rt_window
    ratios = load_ghost_fractions()
    channel, template = _channels(analytes)
    work: dict[str, list] = {}
    for hom, (q1, q3) in channel.items():
        c = get_transition(chroms, q1, q3)
        if c is not None:
            work[hom] = [c, c.times, c.intensities.astype(float).copy(), False]
    sources: dict[str, list[Analyte]] = {}
    targets: dict[str, list[tuple[float, str, str]]] = {}
    for a in analytes:
        if a.role in ("native", "IS_204") and a.homolog in work and a.label in anchors:
            targets.setdefault(a.homolog, []).append((anchors[a.label], a.role, a.label))
            sources.setdefault(a.homolog, []).append(a)

    f204 = dict(DEFAULT_IS_FRACTION)
    f204.update(is_fraction(chroms, analytes, anchors, params))
    f204.update(scale.fractions)

    ghosts: list[Ghost] = []
    fences: list[Analyte] = []
    fence_anchors: dict[str, float] = {}
    for hom in reversed(HOMOLOGS):
        if hom not in sources:
            continue
        _, st, sy, _ = work[hom]
        seen: set[int] = set()
        for a in sorted(sources[hom], key=lambda a: anchors[a.label]):
            # One apex is one source: the first analyte (by anchor) to find it owns it, whether
            # or not it has a ratio of its own.
            pk = _peak(st, sy, anchors[a.label], window)
            if pk is None or pk["i0"] in seen:
                continue
            seen.add(pk["i0"])
            own = {1: 1.0, 2: 1.0} if a.role == "IS_204" else ratios.get(a.label, {})
            for k in STEPS.get(hom, ()):
                ratio = own.get(k)
                if ratio is None or ratio <= 0:
                    continue
                lower = HOMOLOGS[HOMOLOGS.index(hom) - k]
                if lower not in work:
                    continue
                predicted = ratio * f204[k]
                _, tt, ty, _ = work[lower]
                res = _fit(st, sy, pk, tt, ty)
                if res is None or res["r"] < SHAPE_MIN_R or res["f"] <= 0:
                    continue  # the trace does not confirm a copy
                f_use = min(res["f"], predicted)
                if f_use * pk["h"] < SUBTRACT_MIN_SN * res["noise"]:
                    continue
                t_g = pk["t0"] + res["shift"]
                seg_area = float(np.sum((res["seg"][1:] + res["seg"][:-1]) * np.diff(st)) / 2.0) * 60.0
                near = min(targets.get(lower, []), key=lambda x: abs(t_g - x[0]), default=None)
                if near is None or abs(t_g - near[0]) > COELUTION:
                    if res["f"] * pk["h"] < FENCE_MIN_SN * res["noise"] or lower not in template:
                        continue
                    label = f"{FENCE_PREFIX}{a.label}>{lower}@{t_g:.3f}"
                    if label in fence_anchors:
                        continue
                    fences.append(dataclasses.replace(template[lower], label=label, cal_rt=t_g,
                                                      order=900 + len(fences)))
                    fence_anchors[label] = t_g
                    work[lower][3] = work[lower][3] or "fence"
                    ghosts.append(Ghost(a.label, lower, k, round(t_g, 4), res["f"], predicted,
                                        res["r"], f_use * seg_area, FENCED))
                    continue
                if near[1] != "native" or near[2] in NO_SUBTRACTION_UNDER:
                    continue  # never under an internal standard, nor under 58 or 159
                t_lo, t_hi = st[pk["lo"]], st[pk["hi"]]
                m = (tt >= t_lo - 0.01) & (tt <= t_hi + 0.01)
                copy = f_use * _hermite(tt[m] - res["shift"], st, res["seg"])
                before = ty[m].copy()
                floor = res["base"] - FLOOR_NOISE * res["noise"]
                ty[m] = np.maximum(before - copy, np.minimum(before, floor))
                work[lower][3] = True
                ghosts.append(Ghost(a.label, lower, k, round(t_g, 4), res["f"], predicted, res["r"],
                                    f_use * seg_area, SUBTRACTED))

    out = dict(chroms)
    keys: set[tuple[float, float]] = set()
    for hom, (c, t, y, changed) in work.items():
        if not changed:
            continue
        keys.add(_key(c.q1, c.q3))
        if changed is True:
            for key, val in chroms.items():
                if val is c:
                    out[key] = Chromatogram(q1=c.q1, q3=c.q3, times=t, intensities=y)
                    break
    return GhostPlan(chroms=out, fences=fences, fence_anchors=fence_anchors, ghosts=ghosts, keys=keys)


# --------------------------------------------------------------------------------------
# Integrating with the copies kept out
# --------------------------------------------------------------------------------------


def keeps_corrected(raw_area: float, corrected_area: float, ghost_area: float) -> bool:
    """Whether an analyte takes its corrected result: no larger, and smaller by no more than
    :data:`LOSS_BOUND` times the copies in its window (plus the tolerance)."""
    tol = LOSS_TOL_REL * abs(raw_area) + LOSS_TOL_ABS
    return raw_area - LOSS_BOUND * ghost_area - tol <= corrected_area <= raw_area + tol


def integrate_with_ghosts(
    chroms: Mapping[tuple[float, float], Chromatogram],
    analytes: list[Analyte],
    anchors: dict[str, float],
    params: IntegrationParams,
    unresolved: frozenset[tuple[str, str]] = frozenset(),
    scale: GhostScale | None = None,
) -> tuple[list[PeakResult], list[Ghost]]:
    """``integrate_run`` with ghost peaks kept out (see the module docstring).

    ``scale`` is the run's PCB 204 copy (:func:`run_scale`). ``None`` turns the correction off
    and returns exactly ``integrate_run``'s result. Returns the peaks, in ``integrate_run``'s
    order, and the copies found.
    """
    raw = integrate_run(chroms, analytes, anchors, params, unresolved=unresolved)
    if scale is None:
        return raw, []
    plan = find_ghosts(chroms, analytes, anchors, params, scale)
    if not plan.ghosts or not plan.keys:
        return raw, plan.ghosts
    redo = [a for a in analytes if _key(a.q1, a.q3) in plan.keys]
    corrected = integrate_run(plan.chroms, redo + plan.fences, {**anchors, **plan.fence_anchors},
                              params, unresolved=unresolved)
    by_label = {p.label: p for p in corrected if not p.label.startswith(FENCE_PREFIX)}
    homolog_of = {a.label: a.homolog for a in analytes}
    out: list[PeakResult] = []
    for q in raw:
        p = by_label.get(q.label)
        if p is None:
            out.append(q)
            continue
        hom = homolog_of.get(q.label)
        g = sum(x.area for x in plan.ghosts
                if x.channel == hom and q.rt_start - 0.02 <= x.rt <= q.rt_end + 0.02)
        if not keeps_corrected(q.area, p.area, g):
            out.append(q)
            continue
        if abs(p.area - q.area) > 1e-6 * max(1.0, abs(q.area)) or p.detected != q.detected:
            p.flags.append(GHOST_FLAG)
        out.append(p)
    return out, plan.ghosts
