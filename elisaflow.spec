"""Build a windowed application folder with protected assay configuration."""

from importlib.util import find_spec
from pathlib import Path
from shutil import copy2, copytree, ignore_patterns

from PyInstaller.utils.hooks import collect_data_files, collect_submodules, copy_metadata

from elisaflow.configuration.protection import make_read_only

project_dir = Path(SPECPATH)
documentation = project_dir / "build" / "docs"
if not (documentation / "index.html").is_file():
    raise SystemExit("Bouw eerst de handleiding met Sphinx of gebruik scripts/build_windows.ps1.")
qt_runtime = Path(find_spec("shiboken6").origin).parent / "concrt140.dll"
runtime_metadata = [
    entry
    for package in ("elisaflow", "numpy", "pandas", "scipy", "matplotlib", "reportlab")
    for entry in copy_metadata(package)
]

analysis = Analysis(
    [str(project_dir / "src" / "elisaflow" / "__main__.py")],
    pathex=[str(project_dir / "src")],
    binaries=[(str(qt_runtime), ".")],
    datas=[(str(project_dir / "src" / "elisaflow" / "assets"), "elisaflow/assets"), *runtime_metadata, *collect_data_files("tzdata")],
    hiddenimports=collect_submodules("tzdata"),
    hooksconfig={"matplotlib": {"backends": ["Agg"]}},
    # Optional imports in pandas/SciPy can otherwise pull in unrelated packages
    # from a developer's environment. ElisaFlow uses CPU arrays and CSV/XML input.
    excludes=[
        "pymupdf", "fitz", "sphinx", "ruff", "black", "pytest",
        "cupy", "cupyx", "numba", "llvmlite",
        "IPython", "ipykernel", "ipywidgets", "jupyter_client", "jupyter_core",
        "openpyxl", "xlsxwriter", "xlrd", "lxml",
        "tkinter", "_tkinter",
        "PySide6.QtNetwork", "PySide6.QtQml", "PySide6.QtQuick",
        "PySide6.QtOpenGL", "PySide6.QtOpenGLWidgets", "PySide6.QtPrintSupport",
    ],
)
# The desktop uses raster Qt Widgets and reports use Agg. These optional
# resources are not used by either workflow; retain all TrueType math fonts.
analysis.binaries = [
    entry for entry in analysis.binaries
    if Path(entry[0]).name.lower() != "opengl32sw.dll"
]
unused_data = (
    "matplotlib/mpl-data/sample_data/",
    "matplotlib/mpl-data/fonts/afm/",
    "matplotlib/mpl-data/fonts/pdfcorefonts/",
    "PySide6/translations/",
)
analysis.datas = [
    entry for entry in analysis.datas
    if not entry[0].replace("\\", "/").startswith(unused_data)
    or entry[0].replace("\\", "/") == "PySide6/translations/qtbase_nl.qm"
]
pyz = PYZ(analysis.pure)
exe = EXE(
    pyz,
    analysis.scripts,
    [],
    exclude_binaries=True,
    name="ElisaFlow",
    icon=str(project_dir / "src" / "elisaflow" / "assets" / "elisaflow.ico"),
    console=False,
    upx=False,
)
bundle = COLLECT(
    exe,
    analysis.binaries,
    analysis.datas,
    name="ElisaFlow",
    upx=False,
)

# Keep application paths, assay settings and templates beside the executable.
bundle_dir = Path(DISTPATH) / "ElisaFlow"
copy2(project_dir / "app.toml", bundle_dir / "app.toml")
copy2(project_dir / "LICENSE", bundle_dir / "LICENSE")
copytree(documentation, bundle_dir / "docs", ignore=ignore_patterns(".doctrees", ".buildinfo"))
(bundle_dir / "resultaten").mkdir(exist_ok=True)
for directory, pattern in (("configs", "*.json"), ("protocols", "*.skax")):
    target = bundle_dir / directory
    target.mkdir(exist_ok=True)
    for source in sorted((project_dir / directory).glob(pattern)):
        destination = target / source.name
        copy2(source, destination)
        if directory == "configs":
            make_read_only(destination)
