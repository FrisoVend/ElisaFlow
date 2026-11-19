"""Select processing options from the immutable SKAX wavelength list."""

from PySide6.QtCore import QSignalBlocker, Signal
from PySide6.QtWidgets import QComboBox, QLabel, QLineEdit, QVBoxLayout, QWidget

from elisaflow.configuration.assay_settings import MeasurementSettings, SaturationSettings
from elisaflow.models import CalculationMethod
from elisaflow.ui.editors.fields import _number
from elisaflow.ui.theme import FormSection, helper_text


class SaturationEditor(FormSection):
    changed = Signal()

    def __init__(self, wavelength: int, settings: SaturationSettings, editable: bool, wavelengths: tuple[int, ...]):
        super().__init__(f"Verzadiging bij {wavelength} nm")
        self._wavelength = wavelength
        self._wavelengths = wavelengths
        form = self.form
        self.mode = QComboBox()
        self.mode.addItem("Uitsluiten", "exclude")
        self.mode.addItem("Berekenen uit een andere golflengte", "calculate")
        self.mode.setCurrentIndex(self.mode.findData(settings.mode))
        form.addRow("Verzadigde metingen:", self.mode)
        self.source_wavelength_nm = QComboBox()
        for source in wavelengths:
            if source != wavelength:
                self.source_wavelength_nm.addItem(f"{source} nm", source)
        if settings.source_wavelength_nm is not None:
            self.source_wavelength_nm.setCurrentIndex(self.source_wavelength_nm.findData(settings.source_wavelength_nm))
        self.slope = QLineEdit("" if settings.slope is None else str(settings.slope))
        self.intercept = QLineEdit("" if settings.intercept is None else str(settings.intercept))
        form.addRow("Brongolflengte:", self.source_wavelength_nm)
        form.addRow("Helling:", self.slope)
        form.addRow("Intercept:", self.intercept)
        self.equation = helper_text("Ext doel = intercept + helling × Ext bron (na eventuele referentieaftrek, vóór blancoaftrek).")
        form.addRow(self.equation)
        self.setEnabled(editable)
        self.mode.currentIndexChanged.connect(self._refresh)
        self.source_wavelength_nm.currentIndexChanged.connect(lambda: self.changed.emit())
        self.slope.textChanged.connect(lambda: self.changed.emit())
        self.intercept.textChanged.connect(lambda: self.changed.emit())
        self._refresh()

    def set_reference(self, reference: int | None) -> None:
        self.setVisible(self._wavelength != reference)
        sources = tuple(value for value in self._wavelengths if value not in (self._wavelength, reference))
        current = tuple(self.source_wavelength_nm.itemData(index) for index in range(self.source_wavelength_nm.count()))
        if current != sources:
            selected = self.source_wavelength_nm.currentData()
            with QSignalBlocker(self.source_wavelength_nm):
                self.source_wavelength_nm.clear()
                for source in sources:
                    self.source_wavelength_nm.addItem(f"{source} nm", source)
                self.source_wavelength_nm.setCurrentIndex(max(0, self.source_wavelength_nm.findData(selected)))

    def _refresh(self) -> None:
        calculate = self.mode.currentData() == "calculate"
        for widget in (self.source_wavelength_nm, self.slope, self.intercept, self.equation):
            self.layout().setRowVisible(widget, calculate)
        self.changed.emit()

    def selection(self) -> SaturationSettings:
        if self.mode.currentData() == "exclude":
            return SaturationSettings()
        return SaturationSettings(
            "calculate", self.source_wavelength_nm.currentData(),
            _number(self.slope, "Helling"), _number(self.intercept, "Intercept"),
        )


class MeasurementEditor(QWidget):
    changed = Signal()

    def __init__(self, settings: MeasurementSettings, editable: bool, *, wavelengths_nm: tuple[int, ...]):
        super().__init__()
        self._editable = editable
        self._method = CalculationMethod.CALIBRATION
        self._wavelengths = wavelengths_nm
        self._preferred = settings.preferred_wavelength_nm or 450
        self._specified_saturation = set(settings.saturation)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(16)
        general = FormSection("Meetinstellingen")
        layout.addWidget(general)
        form = general.form
        self.mode = QComboBox()
        self.mode.addItem("Eén meetresultaat", "single")
        self.mode.addItem("Twee golflengten; omschakelen op concentratie", "switch_on_concentration")
        self.mode.setCurrentIndex(self.mode.findData(settings.mode))

        self.wavelengths_nm = QLabel(", ".join(map(str, wavelengths_nm)))
        form.addRow("Golflengten uit SKAX (nm):", self.wavelengths_nm)
        form.addRow("Meetmethode:", self.mode)
        self.reference_wavelength_nm = QComboBox()
        self.reference_wavelength_nm.addItem("Geen", None)
        for wavelength in wavelengths_nm:
            self.reference_wavelength_nm.addItem(f"{wavelength} nm", wavelength)
        self.reference_wavelength_nm.setCurrentIndex(max(0, self.reference_wavelength_nm.findData(settings.reference_wavelength_nm)))
        self.reference_wavelength_nm.setEnabled(editable)
        form.addRow("Referentiegolflengte aftrekken:", self.reference_wavelength_nm)
        self.preferred_wavelength_nm = QComboBox()
        self.switching_settings = FormSection("Omschakeling tussen golflengten")
        layout.addWidget(self.switching_settings)
        switching_form = self.switching_settings.form
        switching_form.addRow("Voorkeursgolflengte:", self.preferred_wavelength_nm)
        self.cutoff_concentration = QLineEdit("" if settings.cutoff_concentration is None else str(settings.cutoff_concentration))
        self.cutoff_concentration.setPlaceholderText("Verplicht bij curvewisseling")
        self.cutoff_concentration.setReadOnly(not editable)
        switching_form.addRow("Omschakelconcentratie (vóór verdunningscorrectie):", self.cutoff_concentration)
        self.saturation = {}
        for wavelength in wavelengths_nm:
            editor = SaturationEditor(wavelength, settings.saturation.get(wavelength, SaturationSettings()), editable, wavelengths_nm)
            self.saturation[wavelength] = editor
            layout.addWidget(editor)
        layout.addStretch()
        self.mode.currentIndexChanged.connect(lambda: self._refresh(notify=True))
        self.preferred_wavelength_nm.currentIndexChanged.connect(lambda: self.changed.emit())
        self.reference_wavelength_nm.currentIndexChanged.connect(lambda: self._refresh(notify=True))
        self.cutoff_concentration.textChanged.connect(lambda: self._refresh(notify=True))
        for editor in self.saturation.values():
            editor.changed.connect(lambda: self._refresh(notify=True))
        self._refresh()

    def set_method(self, method: CalculationMethod) -> None:
        self._method = method
        self._refresh()

    def _refresh(self, *, notify: bool = False) -> None:
        switching = self._method == CalculationMethod.CALIBRATION and self.mode.currentData() == "switch_on_concentration"
        self.mode.setEnabled(self._editable and self._method == CalculationMethod.CALIBRATION)
        self.switching_settings.setVisible(switching)
        reference = self.reference_wavelength_nm.currentData()
        analysis = tuple(value for value in self._wavelengths if value != reference)
        for editor in self.saturation.values():
            editor.set_reference(reference)
        preferred = self.preferred_wavelength_nm.currentData() or self._preferred
        if analysis and preferred not in analysis:
            preferred = analysis[0]
        current = tuple(self.preferred_wavelength_nm.itemData(index) for index in range(self.preferred_wavelength_nm.count()))
        if current != analysis:
            with QSignalBlocker(self.preferred_wavelength_nm):
                self.preferred_wavelength_nm.clear()
                for wavelength in analysis:
                    self.preferred_wavelength_nm.addItem(f"{wavelength} nm", wavelength)
                self.preferred_wavelength_nm.setCurrentIndex(self.preferred_wavelength_nm.findData(preferred))
        self._preferred = preferred
        self.preferred_wavelength_nm.setEnabled(self._editable)
        if notify:
            self.changed.emit()

    def selection(self) -> MeasurementSettings:
        switching = self._method == CalculationMethod.CALIBRATION and self.mode.currentData() == "switch_on_concentration"
        reference = self.reference_wavelength_nm.currentData()
        saturation = {wavelength: editor.selection() for wavelength, editor in self.saturation.items() if wavelength != reference}
        return MeasurementSettings(
            mode="switch_on_concentration" if switching else "single",
            preferred_wavelength_nm=self.preferred_wavelength_nm.currentData() if switching else None,
            reference_wavelength_nm=reference,
            cutoff_concentration=_number(self.cutoff_concentration, "Omschakelconcentratie") if switching else None,
            saturation={wavelength: rule for wavelength, rule in saturation.items()
                        if rule.mode != "exclude" or wavelength in self._specified_saturation},
        )

    def to_config(self) -> dict:
        settings = self.selection()
        settings.analysis_wavelengths(self._wavelengths)
        return settings.to_dict()
