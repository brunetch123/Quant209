"""Ghost peaks (``core.ghosts``): a congener's in-source fragment copies on lower homologs' channels.

Synthetic traces, built so each test isolates one rule: a copy under an analyte is subtracted, a
copy between analytes is fenced off from a neighbour's window, a copy the trace cannot confirm is
left alone, nothing is subtracted under an internal standard or under PCB 58 or PCB 159, a
fragment the masses do not allow is never looked for, and the bound keeps an analyte's first
result when the second moved for another reason. The validation against analysts is in
``docs/METHODS.md`` section 2.10.
"""

import numpy as np
import pytest

from quant209.core import ghosts as G
from quant209.core.config import Analyte, load_canonical_peaks, load_ghost_fractions
from quant209.core.integrate import integrate_run
from quant209.core.mzml import Chromatogram
from quant209.core.peaks import IntegrationParams

CHANNEL = {
    "di": (222.0, 152.1), "tri": (256.0, 186.0), "tetra": (291.9, 222.0), "penta": (325.9, 255.9),
    "hexa": (359.8, 289.9), "hepta": (393.8, 323.9), "octa": (429.8, 359.8), "nona": (463.7, 393.8),
}
#: A 0.4 s acquisition cycle over 20 minutes.
T = np.arange(40.0, 60.0, 0.4 / 60.0)
SIGMA = 0.02
SCALE = G.GhostScale(fractions={1: 0.04, 2: 0.15}, measured_in={1: "2 blanks", 2: "2 blanks"})


def gauss(t0: float, height: float, sigma: float = SIGMA) -> np.ndarray:
    return height * np.exp(-0.5 * ((T - t0) / sigma) ** 2)


def area(height: float, sigma: float = SIGMA) -> float:
    """A Gaussian's area in counts*s (times in minutes)."""
    return height * sigma * np.sqrt(2 * np.pi) * 60.0


def trace(*components: np.ndarray, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    y = 100.0 + rng.normal(0.0, 5.0, T.size)
    for c in components:
        y = y + c
    return np.clip(y, 0.0, None)


def chroms(**traces: np.ndarray) -> dict:
    return {CHANNEL[hom]: Chromatogram(q1=CHANNEL[hom][0], q3=CHANNEL[hom][1], times=T, intensities=y)
            for hom, y in traces.items()}


def analyte(label: str, homolog: str, rt: float, role: str = "native") -> Analyte:
    q1, q3 = CHANNEL[homolog]
    return Analyte(label=label, q1=q1, q3=q3, homolog=homolog, order=1, role=role,
                   is_norm="204", cal_mass=10.0, cal_rt=rt)


def setup(analytes: list[Analyte]) -> tuple[list[Analyte], dict[str, float]]:
    return analytes, {a.label: a.cal_rt for a in analytes}


@pytest.fixture
def params() -> IntegrationParams:
    return IntegrationParams()


def by_label(peaks):
    return {p.label: p for p in peaks}


# --------------------------------------------------------------------------------------
# The table
# --------------------------------------------------------------------------------------


def test_ghost_fractions_name_real_congeners_on_steps_their_masses_allow():
    table = load_ghost_fractions()
    homolog = {a.label: a.homolog for a in load_canonical_peaks() if a.role == "native"}
    assert table, "the packaged table is empty"
    for label, steps in table.items():
        assert label in homolog, f"{label} is not a native analyte"
        assert set(steps) <= set(G.STEPS.get(homolog[label], ())), f"{label}: {sorted(steps)}"
        assert all(r > 0 for r in steps.values())


# --------------------------------------------------------------------------------------
# The scale: PCB 204's copy
# --------------------------------------------------------------------------------------


def test_is_fraction_reads_204s_copy_on_both_lower_channels(params):
    analytes, anchors = setup([analyte("204", "octa", 50.0, role="IS_204"),
                               analyte("172", "hepta", 52.0), analyte("156+157", "hexa", 52.0)])
    c = chroms(octa=trace(gauss(50.0, 20000)), hepta=trace(gauss(50.0, 800), seed=1),
               hexa=trace(gauss(50.0, 3000), seed=2))
    got = G.is_fraction(c, analytes, anchors, params)
    assert got[1] == pytest.approx(0.04, rel=0.03)
    assert got[2] == pytest.approx(0.15, rel=0.03)


def test_run_scale_is_the_blanks_median_with_references_as_fallback():
    scale = G.run_scale([{1: 0.04, 2: 0.15}, {1: 0.05, 2: 0.16}, {1: 0.045}], [{1: 0.9, 2: 0.9}])
    assert scale.fractions == {1: pytest.approx(0.045), 2: pytest.approx(0.155)}
    assert scale.measured_in == {1: "3 blanks", 2: "2 blanks"}
    only_refs = G.run_scale([{1: 0.04}], [{2: 0.2}])
    assert only_refs.fractions == {1: 0.04, 2: 0.2}
    assert only_refs.measured_in[2] == "1 reference injection"
    assert G.run_scale([], []).fractions == {}


# --------------------------------------------------------------------------------------
# Keeping copies out
# --------------------------------------------------------------------------------------


def test_correction_off_is_exactly_integrate_run(params):
    analytes, anchors = setup([analyte("204", "octa", 50.0, role="IS_204"),
                               analyte("156+157", "hexa", 50.04)])
    c = chroms(octa=trace(gauss(50.0, 20000)), hexa=trace(gauss(50.0, 3000), seed=2))
    peaks, ghosts = G.integrate_with_ghosts(c, analytes, anchors, params, scale=None)
    plain = integrate_run(c, analytes, anchors, params)
    assert ghosts == []
    assert [(p.label, p.area, p.flags) for p in peaks] == [(p.label, p.area, p.flags) for p in plain]


def test_copy_under_an_analyte_is_subtracted(params):
    # PCB 204's Cl2-loss copy sits 0.04 min from 156+157's expected retention time, on a hexa
    # channel with nothing else there: integrated as acquired it is reported as 156+157.
    analytes, anchors = setup([analyte("204", "octa", 50.0, role="IS_204"),
                               analyte("156+157", "hexa", 50.04)])
    c = chroms(octa=trace(gauss(50.0, 20000)), hexa=trace(gauss(50.0, 3000), seed=2))
    raw = by_label(integrate_run(c, analytes, anchors, params))["156+157"]
    peaks, ghosts = G.integrate_with_ghosts(c, analytes, anchors, params, scale=SCALE)
    got = by_label(peaks)["156+157"]
    assert raw.detected and raw.area == pytest.approx(area(3000), rel=0.1)
    assert got.area < 0.1 * raw.area
    assert not got.detected
    assert G.GHOST_FLAG in got.flags
    sub = [g for g in ghosts if g.channel == "hexa"]
    assert len(sub) == 1 and sub[0].source == "204" and sub[0].action == G.SUBTRACTED and sub[0].step == 2
    assert sub[0].fraction == pytest.approx(0.15, rel=0.05)
    # The internal standard's own peak is untouched.
    assert by_label(peaks)["204"].area == by_label(integrate_run(c, analytes, anchors, params))["204"].area


def test_copy_between_analytes_is_fenced_off_a_neighbours_window(params):
    # A real 153+168 peak 0.12 min before the copy, whose next neighbour is absent: its window
    # runs through the absent neighbour's territory and takes the copy with it.
    analytes, anchors = setup([analyte("204", "octa", 50.0, role="IS_204"),
                               analyte("153+168", "hexa", 49.88), analyte("141", "hexa", 50.40)])
    sigma = 0.025
    c = chroms(octa=trace(gauss(50.0, 20000)),
               hexa=trace(gauss(49.88, 6000, sigma) + gauss(50.0, 3000), seed=3))
    true = area(6000, sigma)
    raw = by_label(integrate_run(c, analytes, anchors, params))["153+168"]
    peaks, ghosts = G.integrate_with_ghosts(c, analytes, anchors, params, scale=SCALE)
    got = by_label(peaks)
    assert raw.area > 1.2 * true
    assert got["153+168"].area == pytest.approx(true, rel=0.05)
    assert G.GHOST_FLAG in got["153+168"].flags
    assert not any(label.startswith(G.FENCE_PREFIX) for label in got)
    assert [g.action for g in ghosts if g.channel == "hexa"] == [G.FENCED]


def test_copy_the_trace_does_not_confirm_is_left_alone(params):
    # 156+157 has a real, broader peak of its own over the copy: the hexa trace does not follow
    # 204's shape, so nothing is subtracted and the analysts' whole peak stands.
    analytes, anchors = setup([analyte("204", "octa", 50.0, role="IS_204"),
                               analyte("156+157", "hexa", 50.04)])
    c = chroms(octa=trace(gauss(50.0, 20000)),
               hexa=trace(gauss(50.05, 12000, 0.05) + gauss(50.0, 3000), seed=4))
    raw = by_label(integrate_run(c, analytes, anchors, params))["156+157"]
    peaks, ghosts = G.integrate_with_ghosts(c, analytes, anchors, params, scale=SCALE)
    got = by_label(peaks)["156+157"]
    assert not [g for g in ghosts if g.channel == "hexa"]
    assert got.area == raw.area
    assert G.GHOST_FLAG not in got.flags


def test_nothing_is_subtracted_under_an_internal_standard(params):
    # PCB 206's HCl-loss copy lands on the octa channel 0.03 min from an analyte's expected
    # retention time. Under a congener it is subtracted; under the internal standard 204, which
    # the analysts integrate whole, it is left.
    c = chroms(nona=trace(gauss(50.03, 200000), seed=5),
               octa=trace(gauss(50.0, 1000) + gauss(50.03, 0.0585 * 200000), seed=6))
    for label, role in (("196", "native"), ("204", "IS_204")):
        analytes, anchors = setup([analyte("206", "nona", 50.03), analyte(label, "octa", 50.0, role=role),
                                   analyte("200", "octa", 48.0)])
        raw = by_label(integrate_run(c, analytes, anchors, params))[label]
        peaks, ghosts = G.integrate_with_ghosts(c, analytes, anchors, params, scale=SCALE)
        subtracted = [g for g in ghosts if g.channel == "octa" and g.action == G.SUBTRACTED]
        if role == "native":
            assert subtracted and subtracted[0].source == "206"
            assert by_label(peaks)[label].area < 0.5 * raw.area
        else:
            assert not subtracted
            assert by_label(peaks)[label].area == raw.area


@pytest.mark.parametrize("source, parent, exempt, other, lower, far", [
    ("95", "penta", "58", "57", "tetra", "52"),
    ("174", "hepta", "159", "162", "hexa", "153+168"),
])
def test_nothing_is_subtracted_under_58_or_159(params, source, parent, exempt, other, lower, far):
    # PCB 95's HCl-loss copy elutes just before PCB 58 and PCB 174's with PCB 159, and their
    # reported peaks usually hold it: left there. The same copy under any other congener goes.
    ratio = load_ghost_fractions()[source][1]
    c = chroms(**{parent: trace(gauss(50.0, 200000), seed=10),
                  lower: trace(gauss(50.0, ratio * SCALE.fractions[1] * 200000), seed=11)})
    for label in (exempt, other):
        analytes, anchors = setup([analyte(source, parent, 50.0), analyte(label, lower, 50.05),
                                   analyte(far, lower, 47.0)])
        raw = by_label(integrate_run(c, analytes, anchors, params))[label]
        peaks, ghosts = G.integrate_with_ghosts(c, analytes, anchors, params, scale=SCALE)
        subtracted = [g for g in ghosts if g.channel == lower and g.action == G.SUBTRACTED]
        assert raw.detected
        if label == exempt:
            assert not subtracted
            assert by_label(peaks)[label].area == raw.area
            assert G.GHOST_FLAG not in by_label(peaks)[label].flags
        else:
            assert subtracted and subtracted[0].source == source
            assert by_label(peaks)[label].area < 0.5 * raw.area


def test_a_fragment_the_masses_do_not_allow_is_never_looked_for(params):
    # A tetra-CB losing Cl2 cannot give the di-CB product ion, so a peak on the di channel at
    # PCB 52's retention time is a di-CB, whatever its shape.
    analytes, anchors = setup([analyte("52", "tetra", 50.0), analyte("15", "di", 50.02),
                               analyte("37", "tri", 50.02)])
    c = chroms(tetra=trace(gauss(50.0, 20000), seed=7), di=trace(gauss(50.0, 3000), seed=8),
               tri=trace(gauss(50.0, 0.0236 * 20000), seed=9))
    raw = by_label(integrate_run(c, analytes, anchors, params))
    peaks, ghosts = G.integrate_with_ghosts(c, analytes, anchors, params, scale=SCALE)
    assert by_label(peaks)["15"].area == raw["15"].area
    assert {g.channel for g in ghosts} == {"tri"}


def test_bound_keeps_the_first_result_when_the_second_moved_for_another_reason():
    assert G.keeps_corrected(10000.0, 9000.0, 1000.0)  # lost about what the copy brought
    assert not G.keeps_corrected(10000.0, 5000.0, 1000.0)  # lost far more than any copy
    assert not G.keeps_corrected(10000.0, 11000.0, 1000.0)  # a correction only takes away
    assert G.keeps_corrected(10000.0, 10200.0, 0.0)  # within the tolerance
