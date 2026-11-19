"""Qt prompts for selecting an assay and entering the analyst name."""

import sys
from pathlib import Path

from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QInputDialog,
    QMessageBox,
    QWidget,
)

from elisaflow.branding import APP_NAME
from elisaflow.formatting import format_protocol_name
from elisaflow.runs.analyst import validate_analyst_name
from elisaflow.ui.desktop import ProtocolSelectionCancelled
from elisaflow.ui.theme import PageHeader, apply_theme, form_layout, page_layout, surface

# Keep an application created by the command-line workflow alive between calls.
_application: QApplication | None = None


def prompt_run_analyst(parent: QWidget | None = None) -> str | None:
    """Require a nonempty name, allowing the operator to cancel before a run."""
    global _application
    _application = QApplication.instance() or QApplication([sys.argv[0]])
    while True:
        value, accepted = QInputDialog.getText(parent, "Naam analist", "Vul uw naam in:")
        if not accepted:
            return None
        try:
            return validate_analyst_name(value)
        except ValueError as error:
            QMessageBox.warning(parent, "Naam ongeldig", str(error))


class RunSelectionDialog(QDialog):
    def __init__(self, protocols: tuple[Path, ...], parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle(f"{APP_NAME} | Start meting")
        self.resize(520, 300)
        apply_theme(self)

        layout = page_layout(self)
        layout.addWidget(PageHeader("Selecteer een assay", "Het bijbehorende protocol wordt geopend in SkanIt."))
        card = surface()
        form = form_layout(card, inset=True)
        self.protocol = QComboBox()
        self.protocol.setPlaceholderText("Selecteer een assay")
        for path in protocols:
            self.protocol.addItem(format_protocol_name(path.name), path)
        self.protocol.setCurrentIndex(-1)
        form.addRow("Assay:", self.protocol)
        layout.addWidget(card)
        layout.addStretch()

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
        )
        self.start_button = buttons.button(QDialogButtonBox.StandardButton.Ok)
        self.start_button.setText("Start meting")
        self.start_button.setProperty("primary", True)
        self.start_button.setEnabled(False)

        self.annuleer_button = buttons.button(QDialogButtonBox.StandardButton.Cancel)
        self.annuleer_button.setText("Annuleer")

        self.protocol.currentIndexChanged.connect(
            lambda index: self.start_button.setEnabled(index >= 0),
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def selected_protocol(self) -> Path:
        path = self.protocol.currentData()
        if path is None:
            raise RuntimeError("Geen SkanIt protocol geselecteerd.")
        return path


def select_run_protocol(protocols: tuple[Path, ...]) -> Path:
    """Run a modal selector, reusing an existing Qt application when available."""
    global _application
    _application = QApplication.instance() or QApplication([sys.argv[0]])
    dialog = RunSelectionDialog(protocols)
    try:
        if dialog.exec() != QDialog.DialogCode.Accepted:
            raise ProtocolSelectionCancelled("Geen SkanIt protocol geselecteerd.")
        return dialog.selected_protocol()
    finally:
        dialog.deleteLater()
