"""Internal-standard spike amounts and the per-homolog surrogate assignment.

Two internal standards are spiked into every sample just before injection: d-PCB30
("30D", ``IS_d30`` role) normalises the low-chlorination homologs and PCB204 ("204",
``IS_204`` role) normalises the high-chlorination homologs. Which one applies to a given
analyte is recorded in its ``is_norm`` field. Ten additional ``13C``-labelled surrogates are
spiked pre-extraction, one per homolog, to correct each homolog's congeners for
extraction/cleanup recovery. Spike amounts reverse-engineered from the lab workbook
``SamplesList_Standards`` sheet.
"""

from __future__ import annotations

from collections.abc import Sequence

from .config import Analyte

IS_SPIKE_NG = 10.4125  # each internal standard (d-PCB30, PCB204) spiked per sample. The exact
# value comes from the lab workbook's SamplesList_Standards sheet
#
# ONE constant for BOTH internal standards, deliberately collapsed: in this method d-PCB30 and
# PCB204 are spiked at the same amount, so a per-standard table would be two identical numbers
# and an invitation to typo one of them. It is a property of the method, not of the mass
# formula. ``Mass = RRF * Area * (IS_spike / IS_area)`` takes the spike of *this analyte's own*
# internal standard. A future method that spikes the two at different amounts must split this
# into a per-IS lookup (keyed like :func:`internal_standard_label`, off ``Analyte.is_norm``) and
# fix every call site, rather than pick one of the two values here.
SURROGATE_SPIKE_NG = 9.84  # each 13C surrogate spiked per sample

# Each native homolog's recovery is tracked by its 13C-labelled surrogate (by chlorination).
HOMOLOG_SURROGATE: dict[str, str] = {
    "mono": "13C 3",
    "di": "13C 15",
    "tri": "13C 28",
    "tetra": "13C 52",
    "penta": "13C 118",
    "hexa": "13C 153",
    "hepta": "13C 180",
    "octa": "13C 194",
    "nona": "13C 206",
    "deca": "13C 209",
}

_IS_ROLES = ("IS_d30", "IS_204")


def internal_standard_label(is_norm: str, analytes: Sequence[Analyte]) -> str:
    """Resolve the internal-standard analyte label ('30D' or '204') for an ``is_norm`` tag.

    Reads the labelled IS_d30/IS_204 rows straight from the canonical peak table rather than
    hardcoding the label strings, so a relabelled standard in the CSV stays correct here.
    """
    for analyte in analytes:
        if analyte.role in _IS_ROLES and analyte.is_norm == is_norm:
            return analyte.label
    raise ValueError(f"no internal-standard row found for is_norm={is_norm!r}")
