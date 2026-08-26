# PyInstaller build spec for the Quant 209 desktop application.
#
# One *folder* build, not one file. A one-file build unpacks itself into a temporary
# directory on every launch, which on a PySide6 application of this size costs several
# seconds each time -- and cold start under three seconds is a requirement this product was
# rebuilt to meet (see docs/validation/checkpoint-E.md). The folder is zipped for delivery
# and installed once by the installer beside it.
#
#   .venv\Scripts\pyinstaller.exe packaging\quant209.spec --noconfirm
#
# Excluded Qt modules below are the large ones this application never imports. They are
# named explicitly rather than trimmed by hand afterwards so the build is reproducible.

from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files

project = Path(SPECPATH).parent

datas = [
    # The canonical 185-analyte peak table: read at runtime through importlib.resources,
    # so PyInstaller cannot see it by static analysis.
    (str(project / "quant209" / "config" / "canonical_peaks.csv"), "quant209/config"),
    # Combo-box chevrons referenced from the stylesheet.
    (str(project / "quant209" / "app" / "assets"), "quant209/app/assets"),
]
datas += collect_data_files("pyqtgraph")

# Qt modules the application does not use. PySide6-Essentials is already a reduced install;
# these are what remains that is both large and unreferenced.
excluded_qt = [
    "PySide6.QtWebEngineCore", "PySide6.QtWebEngineWidgets", "PySide6.QtWebEngineQuick",
    "PySide6.QtQuick", "PySide6.QtQuick3D", "PySide6.QtQml", "PySide6.Qt3DCore",
    "PySide6.Qt3DRender", "PySide6.QtMultimedia", "PySide6.QtMultimediaWidgets",
    "PySide6.QtBluetooth", "PySide6.QtNfc", "PySide6.QtPositioning", "PySide6.QtSensors",
    "PySide6.QtSerialPort", "PySide6.QtWebSockets", "PySide6.QtWebChannel",
    "PySide6.QtCharts", "PySide6.QtDataVisualization", "PySide6.QtPdf", "PySide6.QtPdfWidgets",
    "PySide6.QtDesigner", "PySide6.QtHelp", "PySide6.QtTest", "PySide6.QtSql",
]

analysis = Analysis(
    [str(project / "quant209" / "app" / "main.py")],
    pathex=[str(project)],
    binaries=[],
    datas=datas,
    hiddenimports=[
        # The six workflow pages are imported by *name* at runtime (see
        # quant209/app/workflow/__init__.py) so that a page's heavy dependencies load only
        # when it is first opened -- which is what keeps cold start under three seconds.
        # Static analysis cannot see through import_module, so they are named here.
        "quant209.app.workflow.select_page",
        "quant209.app.workflow.assign_page",
        "quant209.app.workflow.convert_page",
        "quant209.app.workflow.process_page",
        "quant209.app.workflow.review_page",
        "quant209.app.workflow.export_page",
        "quant209.app.charts",
        # Reached only through lazy imports inside handlers, so static analysis misses them.
        "quant209.core.pipeline",
        "quant209.core.deconvolve",
        "quant209.report.excel",
        "quant209.report.truth",
        "quant209.report.validation",
        "quant209.convert.msconvert",
        "pyqtgraph",
        # matplotlib selects a backend by *name* at runtime, so nothing imports these
        # modules statically and PyInstaller ships none of them. The validation report is
        # written as vector PDF (report/validation.py calls savefig(format="pdf")) after
        # matplotlib.use("Agg"), so both are needed -- without them the frozen application
        # raises ModuleNotFoundError: matplotlib.backends.backend_pdf the first time an
        # analyst asks for a validation report, which is exactly how this was found.
        "matplotlib.backends.backend_pdf",
        "matplotlib.backends.backend_agg",
    ],
    hookspath=[],
    runtime_hooks=[],
    excludes=[*excluded_qt, "tkinter", "IPython", "jupyter", "notebook", "pytest", "sphinx"],
    noarchive=False,
)

pyz = PYZ(analysis.pure)

exe = EXE(
    pyz,
    analysis.scripts,
    [],
    exclude_binaries=True,
    name="Quant209",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    # console=False: launching must never open a terminal window. Any unhandled exception
    # is shown by the application's own error dialog rather than printed to a console the
    # user does not have.
    console=False,
    icon=str(project / "quant209" / "app" / "assets" / "quant209.ico"),
)

COLLECT(
    exe,
    analysis.binaries,
    analysis.datas,
    strip=False,
    upx=False,
    name="Quant209",
)
