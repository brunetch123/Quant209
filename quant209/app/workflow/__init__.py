"""The six workflow steps, imported one at a time and only when opened.

Why a registry of module names instead of six imports
------------------------------------------------------
The cold-start budget is 3 seconds from double-click to a window the analyst can use. Two
of these pages import the heavy world: Review reaches for pyqtgraph and pandas to draw
chromatograms, and Export reaches for openpyxl. Importing this package must not drag them
in. The pages are therefore named here as strings and resolved with
:func:`importlib.import_module` at first navigation. A session that goes Select -> Process
-> Review never pays for Export at all, and the Select page (the only one on screen at
startup) costs nothing but Qt.

The gate in :mod:`quant209.app.state` guarantees the ordering these lazy imports rely on:
a page cannot be opened before its step is unlocked, so no page is ever constructed against
a state it was not written for.
"""

from __future__ import annotations

from importlib import import_module

from quant209.app.state import (
    STEP_ASSIGN,
    STEP_CONVERT,
    STEP_EXPORT,
    STEP_PROCESS,
    STEP_REVIEW,
    STEP_SELECT,
)

#: step -> (module, class). Indexed by step number, in workflow order.
PAGE_REGISTRY: dict[int, tuple[str, str]] = {
    STEP_SELECT: ("quant209.app.workflow.select_page", "SelectPage"),
    STEP_ASSIGN: ("quant209.app.workflow.assign_page", "AssignPage"),
    STEP_CONVERT: ("quant209.app.workflow.convert_page", "ConvertPage"),
    STEP_PROCESS: ("quant209.app.workflow.process_page", "ProcessPage"),
    STEP_REVIEW: ("quant209.app.workflow.review_page", "ReviewPage"),
    STEP_EXPORT: ("quant209.app.workflow.export_page", "ExportPage"),
}


def load_page_class(step: int) -> type:
    """Import and return the page class for ``step``. The import happens here, once."""
    module_name, class_name = PAGE_REGISTRY[step]
    return getattr(import_module(module_name), class_name)
