"""One desktop window for measurement runs and assay management."""

from pathlib import Path

from elisaflow.ui.run_worker import RunWorker

# Load the worker's calculation dependencies before importing Qt.
# isort: split
from PySide6.QtCore import Qt, Slot
from PySide6.QtGui import QCloseEvent
from PySide6.QtWidgets import (
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from elisaflow.assays.access_control import AccessControl
from elisaflow.assays.management import AssayManagement
from elisaflow.branding import APP_NAME
from elisaflow.configuration.application import ApplicationSettings
from elisaflow.runs.experiments import load_experiment
from elisaflow.ui.assay_settings import AssaySettingsWidget
from elisaflow.ui.desktop import ProtocolSelectionCancelled, open_document, open_documentation, prompt_analyst_name, select_protocol
from elisaflow.ui.theme import PageHeader, apply_theme, helper_text, page_layout, scroll_page, surface


class MainWindow(QMainWindow):
    def __init__(self, config: ApplicationSettings, assays: AssayManagement, access: AccessControl) -> None:
        super().__init__()
        self._config = config
        self._worker: RunWorker | None = None
        self._report: Path | None = None
        self.setWindowTitle(APP_NAME)
        self.resize(920, 680)
        self.setMinimumSize(760, 560)
        apply_theme(self)

        content = QWidget(self)
        content.setObjectName("appContent")
        content_layout = QVBoxLayout(content)
        content_layout.setContentsMargins(0, 0, 0, 0)
        content_layout.setSpacing(0)
        self.setCentralWidget(content)
        self.tabs = QTabWidget(content)
        self.tabs.setObjectName("mainTabs")
        content_layout.addWidget(self.tabs, 1)
        footer = QHBoxLayout()
        footer.setContentsMargins(24, 0, 24, 20)
        footer.addStretch()
        self.documentation_button = QPushButton("Open documentatie", self)
        self.documentation_button.setProperty("utility", True)
        self.documentation_button.clicked.connect(self._open_documentation)
        footer.addWidget(self.documentation_button)
        content_layout.addLayout(footer)
        run_page = QWidget(self)
        run_page.setObjectName("appContent")
        layout = page_layout(run_page)
        layout.addWidget(PageHeader("Start meting"))
        run_card = surface()
        run_layout = QVBoxLayout(run_card)
        run_layout.setContentsMargins(24, 20, 24, 20)
        run_layout.setSpacing(16)
        title = QLabel("Nieuwe meting")
        title.setObjectName("sectionTitle")
        run_layout.addWidget(title)
        run_layout.addWidget(helper_text(
            "Selecteer een assay en voer de meting uit in SkanIt. "
            "ElisaFlow analyseert de resultaten en opent het PDF-rapport."
        ))
        actions = QHBoxLayout()
        self.start_button = QPushButton("Start meting")
        self.start_button.setProperty("primary", True)
        self.start_button.clicked.connect(self._start_run)
        actions.addWidget(self.start_button)
        actions.addStretch()
        run_layout.addLayout(actions)
        layout.addWidget(run_card)

        saved_card = surface()
        saved_layout = QVBoxLayout(saved_card)
        saved_layout.setContentsMargins(24, 20, 24, 20)
        saved_layout.setSpacing(16)
        saved_title = QLabel("Opgeslagen meting")
        saved_title.setObjectName("sectionTitle")
        saved_layout.addWidget(saved_title)
        saved_layout.addWidget(helper_text("Heranalyse verwerkt een opgeslagen meting met de toen gebruikte instellingen."))
        saved_actions = QHBoxLayout()
        self.reanalyze_button = QPushButton("Heranalyse")
        self.reanalyze_button.clicked.connect(self._start_reanalysis)
        saved_actions.addWidget(self.reanalyze_button)
        saved_actions.addStretch()
        saved_layout.addLayout(saved_actions)
        layout.addWidget(saved_card)

        status_card = surface()
        status_layout = QVBoxLayout(status_card)
        status_layout.setContentsMargins(24, 20, 24, 20)
        status_layout.setSpacing(12)
        self.status_title = QLabel("Klaar voor de volgende meting")
        self.status_title.setObjectName("statusTitle")
        status_layout.addWidget(self.status_title)
        self.status = QLabel("Klaar om een meting te starten.")
        self.status.setWordWrap(True)
        self.status.setTextFormat(Qt.TextFormat.PlainText)
        self.status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        status_layout.addWidget(self.status)
        self.report_path = QLineEdit()
        self.report_path.setReadOnly(True)
        self.report_path.setAccessibleName("Locatie van het rapport")
        self.report_path.hide()
        status_layout.addWidget(self.report_path)
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setTextVisible(False)
        self.progress.hide()
        status_layout.addWidget(self.progress)
        status_actions = QHBoxLayout()
        self.cancel_button = QPushButton("Verwerking annuleren")
        self.cancel_button.setToolTip(f"Stopt het wachten en verwerken in {APP_NAME}; een meting in SkanIt wordt niet gestopt.")
        self.cancel_button.clicked.connect(self._cancel_run)
        self.cancel_button.hide()
        status_actions.addWidget(self.cancel_button)
        self.open_report_button = QPushButton("Open rapport")
        self.open_report_button.clicked.connect(self._open_report)
        self.open_report_button.hide()
        status_actions.addWidget(self.open_report_button)
        status_actions.addStretch()
        status_layout.addLayout(status_actions)
        layout.addWidget(status_card)
        layout.addStretch()

        self.assay_settings = AssaySettingsWidget(assays, access, self)
        self.tabs.addTab(scroll_page(run_page), "Start meting")
        self.tabs.addTab(self.assay_settings, "Assay beheer")

    def _set_status(self, title: str, message: str, state: str = "") -> None:
        self.status_title.setText(title)
        self.status_title.setProperty("state", state)
        self.status_title.style().unpolish(self.status_title)
        self.status_title.style().polish(self.status_title)
        self.status.setText(message)

    @Slot()
    def _start_reanalysis(self) -> None:
        self._start_run(reanalysis=True)

    @Slot()
    def _start_run(self, *, reanalysis: bool = False) -> None:
        if self._worker is not None:
            return
        try:
            protocol, experiment_dir = None, None
            if reanalysis:
                if not self._config.results_dir.is_dir():
                    raise FileNotFoundError(f"De resultatenmap bestaat niet: {self._config.results_dir}")
                folder = QFileDialog.getExistingDirectory(
                    self, "Kies een experimentmap voor heranalyse", str(self._config.results_dir),
                )
                if not folder:
                    return
                experiment_dir = Path(folder).resolve()
                load_experiment(experiment_dir)
            else:
                protocol = select_protocol(self._config.protocols_dir)
            analyst = prompt_analyst_name(parent=self)
            if analyst is None:
                return
        except ProtocolSelectionCancelled:
            return
        except (OSError, ValueError, RuntimeError) as error:
            QMessageBox.warning(self, "Heranalyse kon niet gestart worden" if reanalysis else "Meting kon niet gestart worden", str(error))
            return

        self._report = None
        self.report_path.hide()
        self.open_report_button.hide()
        self.start_button.setEnabled(False)
        self.reanalyze_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.cancel_button.show()
        self.progress.show()
        if reanalysis:
            self._set_status("Heranalyse wordt verwerkt", f"Opgeslagen meting wordt opnieuw verwerkt: {experiment_dir.name}.")
        else:
            self._set_status("Meting wordt verwerkt", f"Meting gestart: {protocol.stem}.")
        worker = RunWorker(self._config, analyst, self, protocol_template=protocol, experiment_dir=experiment_dir)
        self._worker = worker
        worker.succeeded.connect(self._run_succeeded)
        worker.failed.connect(self._run_failed)
        worker.cancelled.connect(self._run_cancelled)
        worker.finished.connect(self._run_finished)
        worker.finished.connect(worker.deleteLater)
        worker.start()

    @Slot()
    def _cancel_run(self) -> None:
        if self._worker is None:
            return
        self._worker.cancel()
        self.cancel_button.setEnabled(False)
        self._set_status("Verwerking wordt geannuleerd", "Annuleren aangevraagd. Een lopende verwerking wordt eerst afgerond.")

    @Slot(object)
    def _run_succeeded(self, report: Path) -> None:
        self._report = report
        self._set_status("Rapport gereed", "Meting afgerond. Het PDF-rapport is opgeslagen.", "success")
        self.report_path.setText(str(report))
        self.report_path.setCursorPosition(0)
        self.report_path.setToolTip(str(report))
        self.report_path.show()
        self.open_report_button.show()
        self._open_report()

    @Slot(str)
    def _run_failed(self, message: str) -> None:
        self._set_status("Meting mislukt", f"Meting mislukt: {message}", "error")
        QMessageBox.warning(self, "Meting mislukt", message)

    @Slot()
    def _run_cancelled(self) -> None:
        self._set_status("Verwerking geannuleerd", "Verwerking geannuleerd. Reeds opgeslagen bestanden blijven bewaard.")

    @Slot()
    def _run_finished(self) -> None:
        self._worker = None
        self.start_button.setEnabled(True)
        self.reanalyze_button.setEnabled(True)
        self.cancel_button.hide()
        self.progress.hide()

    @Slot()
    def _open_report(self) -> None:
        if self._report is not None:
            try:
                open_document(self._report)
            except (OSError, ValueError) as error:
                QMessageBox.warning(self, "Rapport kon niet geopend worden", str(error))

    @Slot()
    def _open_documentation(self) -> None:
        try:
            open_documentation()
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, "Documentatie kon niet geopend worden", str(error))

    def closeEvent(self, event: QCloseEvent) -> None:
        if self._worker is not None:
            event.ignore()
            QMessageBox.information(self, "Meting actief", f"Rond de meting af of annuleer deze voordat je {APP_NAME} sluit.")
            return
        super().closeEvent(event)
