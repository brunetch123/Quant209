import quant209  # noqa: F401  # sets BLAS thread-limit env vars before numpy is ever imported, on this memory-constrained host

from pathlib import Path

import pytest

from dataset import run_dir


@pytest.fixture
def repo_root() -> Path:
    """Path to the repository root (parent of the tests/ directory)."""
    return Path(__file__).resolve().parent.parent


@pytest.fixture
def reference_run() -> Path:
    """The analysis-run folder the numerical tests compare against.

    Which folder that is comes from ``tests/dataset.json``. See ``tests/dataset.py``.
    Skips cleanly when no archive is configured or the folder isn't present, so a fresh
    clone runs green without instrument data.
    """
    return run_dir()
