"""Calibration anchors at the end of a transition (``core.calibrate.build_anchors``).

Synthetic calibration traces, built so each test isolates one guard: a copy of a more chlorinated
congener's peak is never taken for an analyte's calibration peak (``copy_peaks``), and a transition
whose anchors disagree with the drift the other transitions show is assigned a second time from
that drift (``_drift_second_look``). Both were measured on runs acquired under a shorter oven
programme, where PCB 169 had been anchored on PCB 170's copy and PCB 189 on PCB 208's.
"""

from dataclasses import replace

import numpy as np

from quant209.core.calibrate import build_anchors, copy_peaks, detect_cal_peaks
from quant209.core.config import load_canonical_peaks
from quant209.core.mzml import Chromatogram
from quant209.core.peaks import IntegrationParams

TIMES = np.round(np.arange(44.0, 64.0, 0.005), 4)
SIGMA = 0.02


def trace(transition, peaks, seed=0):
    """A chromatogram of Gaussian peaks ``(apex, height)`` on a flat, lightly noisy baseline."""
    y = np.full(TIMES.shape, 50.0)
    for apex, height in peaks:
        y += height * np.exp(-0.5 * ((TIMES - apex) / SIGMA) ** 2)
    y += np.random.default_rng(seed).normal(0.0, 5.0, TIMES.size)
    return Chromatogram(q1=transition[0], q3=transition[1], times=TIMES.copy(), intensities=y)


def analytes(**rts):
    """Packaged analytes by label, each given the canonical retention time passed for it."""
    by_label = {a.label: a for a in load_canonical_peaks()}
    return [replace(by_label[label.replace("_", "+")], cal_rt=rt) for label, rt in rts.items()]


def key(analyte):
    return (analyte.q1, analyte.q3)


def at(chrom, index):
    return float(chrom.times[index])


def test_a_small_peak_on_a_more_chlorinated_congeners_apex_is_its_copy():
    params = IntegrationParams()
    hexa = analytes(**{"167": 48.74, "156_157": 50.49, "169": 55.37})
    hepta = analytes(**{"172": 51.60, "191": 53.12, "170": 54.50})
    chroms = {
        key(hepta[0]): trace(key(hepta[0]), [(51.60, 20000), (53.12, 22000), (54.50, 25000)], seed=1),
        # 54.50: PCB 170's copy (900 against 170's 25,000). 53.00: as small, but on no hepta apex.
        # 55.37: PCB 169's own peak.
        key(hexa[0]): trace(key(hexa[0]), [(48.74, 30000), (50.49, 40000), (53.00, 900), (54.50, 900), (55.37, 44000)]),
    }
    copies = copy_peaks(chroms, hexa + hepta, params)
    hexa_chrom = chroms[key(hexa[0])]
    flagged = sorted(round(at(hexa_chrom, i), 2) for i in copies.get(key(hexa[0]), set()))
    assert flagged == [54.50]
    detected = [round(at(hexa_chrom, i), 2) for i in detect_cal_peaks(hexa_chrom, params)]
    assert 53.00 in detected and 55.37 in detected  # seen as peaks, just not as copies


def test_a_large_peak_on_a_parents_apex_is_not_a_copy():
    params = IntegrationParams()
    hexa = analytes(**{"167": 48.74, "156_157": 50.49, "169": 54.50})
    hepta = analytes(**{"172": 51.60, "191": 53.12, "170": 54.50})
    chroms = {
        key(hepta[0]): trace(key(hepta[0]), [(51.60, 20000), (53.12, 22000), (54.50, 25000)], seed=1),
        key(hexa[0]): trace(key(hexa[0]), [(48.74, 30000), (50.49, 40000), (54.50, 30000)]),
    }
    assert not copy_peaks(chroms, hexa + hepta, params).get(key(hexa[0]))


def test_the_last_analyte_is_anchored_on_its_own_peak_not_on_a_copy_at_its_prediction():
    # PCB 169 is predicted at 54.52 min, where PCB 170's copy sits, and its own peak is 0.85 min
    # later. Matched to the copy, it would get the copy's area as its calibration peak.
    params = IntegrationParams()
    hexa = analytes(**{"167": 48.74, "156_157": 50.49, "169": 54.52})
    hepta = analytes(**{"172": 51.60, "191": 53.12, "170": 54.50})
    chroms = {
        key(hepta[0]): trace(key(hepta[0]), [(51.60, 20000), (53.12, 22000), (54.50, 25000)], seed=1),
        key(hexa[0]): trace(key(hexa[0]), [(48.74, 30000), (50.49, 40000), (54.50, 900), (55.37, 44000)]),
    }
    anchors = build_anchors(chroms, hexa + hepta, params)
    assert abs(anchors["169"] - 55.37) <= 0.006
    assert abs(anchors["170"] - 54.50) <= 0.006
    assert abs(anchors["156+157"] - 50.49) <= 0.006


def test_late_analytes_predicted_too_early_are_assigned_again_from_the_other_transitions_drift():
    # The chromatography runs 0.25 min late at 57 min and 0.5 min late from 58 min. The octa
    # transition's last two analytes are predicted at their canonical times, and a stray peak at
    # 57.45 min lets the first assignment hand each of them the peak before its own. The tetra and
    # di transitions' widely spaced peaks are assigned right, agree with each other, and show the
    # drift. (Two transitions are needed for that: with only one, its correct anchors would
    # disagree just as much with the octa transition's misplaced ones.)
    params = IntegrationParams()
    octa = analytes(**{"196": 55.0, "203": 56.5, "195": 57.5, "194": 58.0})
    tetra = analytes(**{"52": 55.5, "48": 57.0, "44_47_65": 58.2, "64": 59.0})
    di = analytes(**{"4": 55.5, "10": 57.0, "9": 58.2, "7": 59.0})
    late = [(55.5, 20000), (57.25, 20000), (58.7, 20000), (59.5, 20000)]
    chroms = {
        key(octa[0]): trace(key(octa[0]), [(55.0, 20000), (56.625, 20000), (57.45, 15000), (57.875, 20000),
                             (58.5, 20000)]),
        key(tetra[0]): trace(key(tetra[0]), late, seed=2),
        key(di[0]): trace(key(di[0]), late, seed=3),
    }
    anchors = build_anchors(chroms, octa + tetra + di, params)
    assert abs(anchors["195"] - 57.875) <= 0.006
    assert abs(anchors["194"] - 58.5) <= 0.006
    assert abs(anchors["203"] - 56.625) <= 0.006
    # The tetra and di transitions agreed with the drift, so they are left exactly as first assigned.
    assert [round(anchors[a.label], 2) for a in tetra + di] == [55.5, 57.25, 58.7, 59.5] * 2


def test_transitions_that_agree_are_left_exactly_as_first_assigned():
    params = IntegrationParams()
    octa = analytes(**{"196": 55.0, "203": 56.5, "195": 57.5, "194": 58.0})
    tetra = analytes(**{"52": 55.5, "48": 57.0, "44_47_65": 58.2, "64": 59.0})
    chroms = {
        key(octa[0]): trace(key(octa[0]), [(55.0, 20000), (56.5, 20000), (57.5, 20000), (58.0, 20000)]),
        key(tetra[0]): trace(key(tetra[0]), [(55.5, 20000), (57.0, 20000), (58.2, 20000), (59.0, 20000)], seed=2),
    }
    anchors = build_anchors(chroms, octa + tetra, params)
    assert {a.label: round(anchors[a.label], 2) for a in octa + tetra} == {
        a.label: a.cal_rt for a in octa + tetra}
