"""The logo assets, and the places that load them.

Every one of these is a file the application or the build reads by name at runtime. A
missing or renamed asset fails silently, as a blank label or a default Windows icon, so
each is checked here rather than left to be noticed on a tester's machine.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from quant209.app import theme

ASSETS = Path(theme.__file__).resolve().parent / "assets"
IMAGES = Path(__file__).resolve().parent.parent / "docs" / "images"
SPEC = Path(__file__).resolve().parent.parent / "packaging" / "quant209.spec"

#: Sizes the Windows shell asks for. 16 is the taskbar and the title bar, 256 is the large
#: icon view and the installer.
ICO_SIZES = {16, 24, 32, 48, 64, 128, 256}


@pytest.mark.parametrize("name", ["mark.svg", "icon.svg", "icon-app.svg", "quant209.ico"])
def test_the_application_assets_are_present(name: str):
    path = ASSETS / name
    assert path.is_file(), f"{name} is missing from {ASSETS}"
    assert path.stat().st_size > 0


@pytest.mark.parametrize("name", ["logo.svg", "logo-dark.svg", "logo.png", "logo-dark.png"])
def test_the_documentation_images_are_present(name: str):
    assert (IMAGES / name).is_file(), f"{name} is missing from docs/images"


@pytest.mark.parametrize("name", ["mark.svg", "icon.svg", "icon-app.svg"])
def test_the_svg_assets_parse(name: str):
    """Qt renders nothing at all for malformed SVG, and reports no error while doing it."""
    root = ET.parse(ASSETS / name).getroot()
    assert root.tag.endswith("svg")
    assert root.get("viewBox"), "no viewBox, so the drawing will not scale"


def test_the_windows_icon_carries_every_size_the_shell_asks_for():
    pytest.importorskip("PIL", reason="Pillow reads the icon's directory")
    from PIL import Image

    with Image.open(ASSETS / "quant209.ico") as ico:
        sizes = {w for w, _ in ico.info["sizes"]}
    assert ICO_SIZES <= sizes, f"missing {sorted(ICO_SIZES - sizes)}"


def test_the_build_spec_names_the_icon():
    """PyInstaller embeds this into the executable, and the shortcuts inherit it."""
    text = SPEC.read_text(encoding="utf-8")
    assert "quant209.ico" in text
    assert "icon=None" not in text


def test_the_installer_points_both_shortcuts_at_the_icon():
    vbs = SPEC.parent / "Install Quant 209.vbs"
    text = vbs.read_text(encoding="utf-8", errors="ignore")
    assert text.count("IconLocation") == 2, "desktop and Start Menu shortcuts both need one"


def test_the_readme_shows_the_logo_in_both_themes():
    readme = (Path(__file__).resolve().parent.parent / "README.md").read_text(encoding="utf-8")
    assert "docs/images/logo.png" in readme
    assert "docs/images/logo-dark.png" in readme
    assert "prefers-color-scheme: dark" in readme


def test_the_sidebar_paints_the_mark(qtbot):
    """A missing or unreadable SVG yields a null pixmap and an empty label, not an error."""
    from quant209.app.sidebar import Sidebar

    sidebar = Sidebar()
    qtbot.addWidget(sidebar)
    pixmap = sidebar.mark.pixmap()
    assert not pixmap.isNull(), "the mark did not render"
    assert pixmap.width() == theme.MARK_WIDTH

    # The rendered shape must keep the drawing's own proportions.
    box = ET.parse(ASSETS / "mark.svg").getroot().get("viewBox").split()
    aspect = float(box[3]) / float(box[2])
    assert abs(pixmap.height() / pixmap.width() - aspect) < 0.02


def test_the_window_icon_is_set(qtbot):
    """Without this the frozen application shows the default Windows executable icon."""
    from PySide6.QtWidgets import QApplication

    from quant209.app.main import configure_app

    app = QApplication.instance()
    configure_app(app)
    assert not app.windowIcon().isNull()
    assert app.windowIcon().availableSizes(), "the icon carries no bitmaps"
