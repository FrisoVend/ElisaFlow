"""Shared report and desktop app branding."""

from pathlib import Path

APP_NAME = "ElisaFlow"
APP_USER_MODEL_ID = "ReinierDiagnostiek.ElisaFlow"

INK = "#263442"
MUTED = "#65717D"
BRAND = "#087C8E"
ACCENT = "#12A59F"
BORDER = "#D3DFE0"
HEADER = "#EAF4F3"
STRIPE = "#F7FAFA"
LOGO_PATH = Path(__file__).with_name("assets") / "reinier_diagnostiek_logo_transparent.png"
ICON_PATH = LOGO_PATH.with_name("elisaflow.ico")
