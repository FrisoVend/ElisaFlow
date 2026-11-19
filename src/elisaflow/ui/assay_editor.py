"""Compose focused assay editor pages and persist a validated configuration."""

from copy import deepcopy
from pathlib import Path

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QMessageBox,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from elisaflow.assays.access_control import AccessControl
from elisaflow.assays.management import AssayManagement
from elisaflow.configuration.assay_settings import AssaySettings, MeasurementSettings, ReportSettings
from elisaflow.formatting import format_protocol_name
from elisaflow.models import CalculationMethod
from elisaflow.ui.editors.calculation import FormulaCalculationEditor
from elisaflow.ui.editors.calibration import CalibrationCalculationEditor
from elisaflow.ui.editors.measurement import MeasurementEditor
from elisaflow.ui.editors.report import ReportEditor
from elisaflow.ui.editors.samples import PatientEditor, SampleCollection, SampleEditor
from elisaflow.ui.editors.validation import ValidationEditor
from elisaflow.ui.theme import PageHeader, apply_theme, editor_tabs, form_layout, helper_text, page_layout, scroll_page, surface


class AssayEditor(QDialog):
    def __init__(
        self,
        settings: AssaySettings | None,
        assays: AssayManagement,
        access: AccessControl,
        parent: QWidget | None = None,
        *,
        protocol_path: Path | None = None,
        read_only: bool = False,
    ) -> None:
        super().__init__(parent)
        self._settings = settings
        self._assays = assays
        self._read_only = read_only
        self._editable = not read_only and access.can_edit_assays()
        self._protocol_path = None
        if settings is None:
            if protocol_path is None:
                raise ValueError("Selecteer een SKAX-protocol om een assay aan te maken.")
            self._protocol_path = Path(protocol_path).expanduser().resolve()
            protocol = assays.prepare_new_assay(self._protocol_path)
            self._protocol_name = self._protocol_path.name
            self._method = CalculationMethod.CALIBRATION if protocol.standards else CalculationMethod.FORMULA
            config = {
                "lot_number": "",
                "measurement": {"mode": "single"},
                "calculation": {
                    "method": self._method.value,
                    **({"fit_model": "4pl"} if self._method == CalculationMethod.CALIBRATION else {"expression": "sample"}),
                },
            }
        else:
            protocol = assays.load_protocol_samples(settings.protocol_name)
            protocol.validate_samples(settings)
            protocol.validate_measurements(settings)
            self._protocol_name = settings.protocol_name
            self._method = settings.calculation.method
            config = settings.to_dict()
        self._initial_config = deepcopy(config)
        self._control_names = protocol.controls
        self._has_blank = protocol.has_blank
        self._wavelengths = protocol.wavelengths_nm
        title = "Nieuw assay" if settings is None else "Bekijk assay" if read_only else "Verander assay"
        self.setWindowTitle(f"{title}: {format_protocol_name(self._protocol_name)}")
        self.resize(1100, 820)
        self.setMinimumSize(860, 640)
        apply_theme(self)
        layout = page_layout(self)
        layout.addWidget(PageHeader(format_protocol_name(self._protocol_name), "Alleen bekijken" if not self._editable else title))
        details = surface()
        form = form_layout(details, inset=True)
        method_labels = {CalculationMethod.FORMULA: "Berekening", CalculationMethod.CALIBRATION: "Standaardcurve"}
        if settings is None:
            self.calculation_method = QComboBox()
            for kind in CalculationMethod:
                self.calculation_method.addItem(method_labels[kind], kind.value)
            self.calculation_method.setCurrentIndex(self.calculation_method.findData(self._method.value))
            self.calculation_method.setEnabled(self._editable)
            form.addRow("Resultaatmethode:", self.calculation_method)
        else:
            form.addRow("Resultaatmethode:", QLabel(method_labels[self._method]))
        self.lot_number = QLineEdit(config["lot_number"])
        self.lot_number.setReadOnly(not self._editable)
        form.addRow("Lotnummer kit:", self.lot_number)
        layout.addWidget(details)
        body = QHBoxLayout()
        body.setSpacing(20)
        self.navigation = QListWidget()
        self.navigation.setObjectName("editorNavigation")
        self.navigation.setFixedWidth(196)
        self.navigation.setAccessibleName("Assay-instellingen")
        panel = surface()
        panel_layout = QVBoxLayout(panel)
        panel_layout.setContentsMargins(20, 12, 20, 20)
        panel_layout.setSpacing(16)
        self.page_title = QLabel()
        self.page_title.setObjectName("sectionTitle")
        panel_layout.addWidget(self.page_title)
        self.pages = QStackedWidget()
        panel_layout.addWidget(self.pages, 1)
        body.addWidget(self.navigation)
        body.addWidget(panel, 1)
        layout.addLayout(body, 1)

        measurement = MeasurementSettings.from_dict(config.get("measurement", {}))
        self.measurement = MeasurementEditor(measurement, self._editable, wavelengths_nm=protocol.wavelengths_nm)
        self._add_page("Meting", self.measurement)
        self._calculation_pages: dict[CalculationMethod, FormulaCalculationEditor | CalibrationCalculationEditor] = {}
        self._calculation_stack = QStackedWidget()
        self.sample_tabs = editor_tabs()
        self.patients = PatientEditor(config.get("patients", {}), self._editable, has_blank=protocol.has_blank)
        self.controls = SampleCollection(
            config.get("controls", {}),
            self._editable,
            names={name: name for name in protocol.controls},
            has_blank=protocol.has_blank,
        )
        self.blank = SampleEditor(config.get("blank", {}), self._editable)
        self.standards = SampleCollection(
            config.get("standards", {}),
            self._editable,
            standard=True,
            names={value: ", ".join(names) for value, names in protocol.standards.items()},
            has_blank=protocol.has_blank,
        )
        for label, editor in (("Patiënten", self.patients), ("Controles", self.controls), ("Blanco", self.blank), ("Standaarden", self.standards)):
            index = self.sample_tabs.addTab(scroll_page(editor), label)
            if editor is self.standards:
                self.sample_tabs.setTabVisible(index, bool(protocol.standards))
        self._add_page("Monsterinstellingen", self.sample_tabs, scroll=False)
        self._add_page("Berekening", self._calculation_stack)
        self.validation = ValidationEditor(config.get("validation", {}), protocol, self._editable)
        self._add_page("QC-criteria", self.validation, scroll=False)
        self.report = ReportEditor(ReportSettings.from_dict(config.get("report", {})), self._editable)
        self._add_page("Rapport", self.report)
        self.navigation.currentRowChanged.connect(self._show_page)
        self.navigation.setCurrentRow(0)
        self._show_calculation(self._method.value)
        self.measurement.changed.connect(self._measurement_changed)
        if settings is None:
            self.calculation_method.currentIndexChanged.connect(lambda index: self._show_calculation(self.calculation_method.itemData(index)))
        actions = QDialogButtonBox.StandardButton.Close
        if not read_only:
            actions |= QDialogButtonBox.StandardButton.Save
        buttons = QDialogButtonBox(actions)
        buttons.button(QDialogButtonBox.StandardButton.Close).setText("Sluit")
        save = buttons.button(QDialogButtonBox.StandardButton.Save)
        if save is not None:
            save.setText("Opslaan")
            save.setEnabled(self._editable)
            save.setProperty("primary", True)
        buttons.accepted.connect(self._save)
        buttons.rejected.connect(self.reject)
        footer = QHBoxLayout()
        footer.addWidget(helper_text("Alleen bekijken · wijzigingen zijn vergrendeld" if not self._editable else "Wijzigingen worden bewaard na Opslaan."), 1)
        footer.addWidget(buttons)
        layout.addLayout(footer)

    def _show_page(self, index: int) -> None:
        self.pages.setCurrentIndex(index)
        self.page_title.setText(self.navigation.item(index).text())

    def _add_page(self, title: str, widget: QWidget, *, scroll: bool = True) -> None:
        self.pages.addWidget(scroll_page(widget) if scroll else widget)
        self.navigation.addItem(title)

    def _show_calculation(self, type_name: str) -> None:
        self._method = CalculationMethod(type_name)
        if self._method not in self._calculation_pages:
            editor = FormulaCalculationEditor if self._method == CalculationMethod.FORMULA else CalibrationCalculationEditor
            options = {"has_blank": self._has_blank} if self._method == CalculationMethod.FORMULA else {}
            page = editor(self._initial_config, self._control_names, self._editable, **options)
            self._calculation_pages[self._method] = page
            self._calculation_stack.addWidget(page)
        self.calculation = self._calculation_pages[self._method]
        self._calculation_stack.setCurrentWidget(self.calculation)
        self.report.set_calibration(self._method == CalculationMethod.CALIBRATION)
        self.measurement.set_method(self._method)
        self._measurement_changed()

    def _measurement_changed(self) -> None:
        try:
            selection = self.measurement.selection()
        except ValueError:
            return
        wavelengths = tuple(value for value in self._wavelengths if value != selection.reference_wavelength_nm)
        for page in self._calculation_pages.values():
            page.set_measurement(selection, wavelengths)
        self.validation.set_measurement(selection, wavelengths)

    def _save(self) -> None:
        if self._read_only:
            return
        config = deepcopy(self._initial_config)
        config["lot_number"] = self.lot_number.text().strip()
        try:
            self._measurement_changed()
            config["measurement"] = self.measurement.to_config()
            config["patients"] = self.patients.to_config()
            config["controls"] = self.controls.to_config()
            config["blank"] = self.blank.to_config()
            config["standards"] = self.standards.to_config()
            config["validation"] = self.validation.to_config()
            config["report"] = self.report.to_config()
            config.update(self.calculation.to_config())
            settings = self._assays.validate_assay(self._protocol_name, config, protocol_path=self._protocol_path)
            if self._settings is None:
                self._assays.create_assay(settings, protocol_path=self._protocol_path)
            else:
                self._assays.save_assay(settings, expected_settings=self._settings)
        except (ValueError, PermissionError, NotImplementedError, OSError) as error:
            QMessageBox.warning(self, "Assay kon niet worden opgeslagen", str(error))
            return
        self.accept()
