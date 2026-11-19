"""Create services and launch the shared desktop window."""

import ctypes
import logging
import sys

from elisaflow.assays.access_control import AccessControl
from elisaflow.assays.management import AssayManagement
from elisaflow.branding import APP_NAME, APP_USER_MODEL_ID
from elisaflow.configuration.application import ApplicationSettings


def launch(config: ApplicationSettings, *, show_settings: bool = False) -> int:
    # Load the data services before Qt, keeping pandas imports out of Qt's import hooks.
    access = AccessControl(config.password_file)
    assays = AssayManagement(config.configs_dir, access, protocols_dir=config.protocols_dir)
    # MainWindow loads the calculation modules before importing Qt.
    try:
        from elisaflow.ui.main_window import MainWindow
    except ModuleNotFoundError as error:
        if error.name and error.name.split(".")[0] in {"PySide6", "shiboken6"}:
            raise RuntimeError('Installeer eerst de UI-afhankelijkheid: python -m pip install -e ".[ui]"') from error
        raise

    from PySide6.QtWidgets import QApplication

    if sys.platform == "win32":
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(APP_USER_MODEL_ID)
    application = QApplication.instance() or QApplication([sys.argv[0]])
    application.setApplicationName(APP_NAME)
    application.setApplicationDisplayName(APP_NAME)
    window = MainWindow(config, assays, access)
    if show_settings:
        window.tabs.setCurrentWidget(window.assay_settings)
    window.show()
    logging.getLogger(__name__).info("Desktop ready")
    return application.exec()
