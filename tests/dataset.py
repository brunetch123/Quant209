"""Where the suite finds a local instrument-data archive, and what things are called there.

Most of this suite is self-contained: it builds its own chromatograms and workbooks and
checks the numerics against values worked out by hand. A smaller set of tests instead runs
against a real analysis run and the analyst workbooks acquired with it, because the property
under test is a relationship between a retention-time table and a real chromatogram, and
neither half survives being synthesised.

That archive is laboratory data. It is not distributed with the source and its sample
identifiers are not published, so no run or sample name appears in this repository. To run
those tests against an archive of your own, copy ``dataset.example.json`` to ``dataset.json``
in this directory and fill in the names. ``dataset.json`` is ignored by git. Without it, on
a fresh clone or in CI, every test that needs instrument data skips, and the rest of the
suite, which is most of it, runs normally.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
CONFIG = HERE / "dataset.json"


@lru_cache(maxsize=1)
def dataset() -> dict[str, str]:
    """The local name map, or an empty mapping when no archive is configured."""
    if not CONFIG.is_file():
        return {}
    return json.loads(CONFIG.read_text(encoding="utf-8"))


def name(key: str) -> str:
    """The local name for ``key``, skipping the test when the archive isn't configured."""
    value = dataset().get(key)
    if not value:
        pytest.skip(f"no local dataset entry for {key!r} (see tests/dataset.py)")
    return value


def run_dir(key: str = "run") -> Path:
    """The run folder ``key`` names, skipping the test when it isn't on disk."""
    path = REPO_ROOT / name(key)
    if not path.is_dir():
        pytest.skip("run data not present")
    return path


def workbook(key: str = "workbook") -> Path:
    """Path to an analyst workbook in the reference run, by dataset key."""
    return run_dir() / f"{name(key)}.xlsx"


class Names:
    """Lazy attribute access to the name map: ``NAMES.sample`` is ``name("sample")``.

    Resolved at use rather than at import so that an unconfigured archive skips the
    individual test that wanted a name, instead of failing collection of the module.
    """

    def __getattr__(self, key: str) -> str:
        if key.startswith("_"):
            raise AttributeError(key)
        return name(key)


#: Shared instance. Import as ``from dataset import NAMES``.
NAMES = Names()
