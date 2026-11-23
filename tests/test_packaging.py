"""The packaging contract: what PyInstaller cannot see, the spec must declare.

These failures are invisible to every other test in this suite. Each module below is
reached only by *name* at runtime. matplotlib picks a backend from a string, and the
workflow pages are imported through ``import_module`` so a page's heavy dependencies load
only when it is first opened. Static analysis therefore ships none of them, the developer
install keeps working because the real packages are on the path, and the fault only appears
on a tester's machine as ModuleNotFoundError.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

SPEC = Path(__file__).resolve().parent.parent / "packaging" / "quant209.spec"

#: Imported by name at runtime, so the spec must name them for a frozen build to work.
REQUIRED_HIDDEN_IMPORTS = [
    # matplotlib selects its backend from a string. The validation report writes vector PDF.
    "matplotlib.backends.backend_pdf",
    "matplotlib.backends.backend_agg",
    # the six workflow pages, loaded through quant209.app.workflow.load_page_class
    "quant209.app.workflow.select_page",
    "quant209.app.workflow.assign_page",
    "quant209.app.workflow.convert_page",
    "quant209.app.workflow.process_page",
    "quant209.app.workflow.review_page",
    "quant209.app.workflow.export_page",
    # imported inside handlers to keep them off the startup path
    "quant209.core.pipeline",
    "quant209.report.excel",
    "quant209.report.validation",
]


@pytest.fixture(scope="module")
def spec_text() -> str:
    if not SPEC.exists():
        pytest.skip("packaging spec not present")
    return SPEC.read_text(encoding="utf-8")


@pytest.mark.parametrize("module", REQUIRED_HIDDEN_IMPORTS)
def test_spec_declares_every_runtime_import(spec_text: str, module: str) -> None:
    assert f'"{module}"' in spec_text, (
        f"{module} is imported by name at runtime, so a frozen build will not contain it "
        f"unless packaging/quant209.spec lists it in hiddenimports"
    )


def test_every_workflow_page_module_is_declared(spec_text: str) -> None:
    """The page registry is the source of truth: a seventh page must be declared too."""
    from quant209.app.workflow import PAGE_REGISTRY

    for module_name, _class_name in PAGE_REGISTRY.values():
        assert f'"{module_name}"' in spec_text, (
            f"{module_name} is in the page registry but not in the spec's hiddenimports"
        )


def test_packaged_data_files_are_declared(spec_text: str) -> None:
    """Every packaged table is read through importlib.resources, never imported, so the
    frozen build only carries the ones the spec names. A table left out works from source
    and fails in the installed application."""
    from importlib import resources

    tables = sorted(
        p.name for p in resources.files("quant209.config").iterdir() if p.name.endswith(".csv")
    )
    assert tables == ["aroclor_abundance.csv", "canonical_peaks.csv", "ghost_fractions.csv"]
    for table in tables:
        assert table in spec_text, f"{table} is packaged but not declared in the spec"
    assert "assets" in spec_text


def test_the_build_is_windowed_and_one_folder(spec_text: str) -> None:
    """console=True would open a terminal behind the app. A one-file build would unpack
    itself on every launch and break the cold-start budget."""
    assert re.search(r"console\s*=\s*False", spec_text), "the app must not open a console"
    assert "COLLECT(" in spec_text, "must be a one-folder build"
