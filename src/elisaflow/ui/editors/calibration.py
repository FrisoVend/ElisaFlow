"""Own calibration curve drafts and their wavelength-specific presentation."""

from PySide6.QtCore import QSignalBlocker
from PySide6.QtWidgets import QLineEdit, QStackedWidget, QVBoxLayout, QWidget

from elisaflow.configuration.assay_settings import MeasurementSettings
from elisaflow.ui.editors.curve import CurveEditor
from elisaflow.ui.editors.fields import _number, _range_fields, _result_config
from elisaflow.ui.theme import FormSection, editor_tabs


class CalibrationCalculationEditor(QWidget):
    """Keep shared result settings and the existing draft for each wavelength."""

    def __init__(self, config: dict, control_names: tuple[str, ...], editable: bool) -> None:
        super().__init__()
        self._editable = editable
        self._control_names = control_names
        self._switching = False
        self._single_wavelength: int | None = None
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(16)
        result = config.get("result", {}) if config.get("calculation", {}).get("method") == "calibration" else {}
        self.result_label = QLineEdit(result.get("label", "Conc"))
        self.result_unit = QLineEdit(result.get("unit", ""))
        self.shared_settings = FormSection("Algemene instellingen")
        layout.addWidget(self.shared_settings)
        shared = self.shared_settings.form
        shared.addRow("Resultaatlabel:", self.result_label)
        shared.addRow("Resultaateenheid:", self.result_unit)
        curve = config.get("calculation", {})
        if curve.get("method") != "calibration":
            curve = {"fit_model": "4pl"}
        self._initial_fits = {int(key): fit for key, fit in curve.get("wavelengths", {}).items()}
        preferred = config.get("measurement", {}).get("preferred_wavelength_nm", 450)
        curve_section = FormSection("Standaardcurve")
        layout.addWidget(curve_section)
        self.curve = CurveEditor(self._initial_fits.get(preferred, curve), control_names, editable=editable)
        self._single_curve_stack = QStackedWidget()
        self._single_curve_stack.addWidget(self.curve)
        curve_section.form.addRow(self._single_curve_stack)
        self.wavelength_tabs = editor_tabs()
        self.wavelength_curves: dict[int, CurveEditor] = {}
        self._active_wavelengths = ()
        curve_section.form.addRow(self.wavelength_tabs)
        layout.addStretch()
        bounds, self.concentration_minimum, self.concentration_maximum = _range_fields(*curve.get("concentration_range", (None, None)))
        shared.addRow("Concentratiebereik:", bounds)
        self._show_curve_mode()
        for field in (self.result_label, self.result_unit, self.concentration_minimum, self.concentration_maximum):
            field.setReadOnly(not editable)
        self.curve.layout_changed.connect(self._curve_layout_changed)
        self.wavelength_tabs.currentChanged.connect(self._curve_layout_changed)

    def set_measurement(self, settings: MeasurementSettings, wavelengths: tuple[int, ...]) -> None:
        self._switching = settings.switching
        if self._switching:
            for wavelength in wavelengths:
                self._curve_for_wavelength(wavelength).layout().setContentsMargins(0, 16, 0, 0)
            labels = tuple(settings.wavelength_label(wavelength) for wavelength in wavelengths)
            if wavelengths != self._active_wavelengths:
                current = self.wavelength_tabs.currentWidget()
                with QSignalBlocker(self.wavelength_tabs):
                    self.wavelength_tabs.clear()
                    for wavelength, label in zip(wavelengths, labels, strict=True):
                        self.wavelength_tabs.addTab(self.wavelength_curves[wavelength], label)
                    if self.wavelength_tabs.indexOf(current) >= 0:
                        self.wavelength_tabs.setCurrentWidget(current)
            for index, label in enumerate(labels):
                self.wavelength_tabs.setTabText(index, label)
            self._active_wavelengths = wavelengths
        else:
            if len(wavelengths) == 1:
                wavelength = wavelengths[0]
                if self._single_wavelength is None and not self.wavelength_curves and not self._initial_fits:
                    self.wavelength_curves[wavelength] = self.curve
                self.curve = self._curve_for_wavelength(wavelength)
                self._single_wavelength = wavelength
            self._single_curve_stack.addWidget(self.curve)
            self._single_curve_stack.setCurrentWidget(self.curve)
            self.curve.layout().setContentsMargins(0, 0, 0, 0)
            self._active_wavelengths = ()
        self._show_curve_mode()

    def _curve_for_wavelength(self, wavelength: int) -> CurveEditor:
        if wavelength not in self.wavelength_curves:
            # Retain incomplete drafts only for the wavelength they describe.
            draft = self.curve.draft()
            if wavelength != self._single_wavelength:
                draft.pop("parameters")
                draft.pop("response_scale")
                draft.pop("preprocessing")
                draft["fixed_parameters"] = False
            config = self._initial_fits.get(wavelength, draft)
            curve = CurveEditor(config, self._control_names, editable=self._editable)
            curve.layout_changed.connect(self._curve_layout_changed)
            self.wavelength_curves[wavelength] = curve
        return self.wavelength_curves[wavelength]

    def _show_curve_mode(self) -> None:
        self._single_curve_stack.setVisible(not self._switching)
        self.wavelength_tabs.setVisible(self._switching)
        self._curve_layout_changed()

    def _curve_layout_changed(self) -> None:
        page = self.wavelength_tabs.currentWidget()
        if page is not None:
            height = page.heightForWidth(self.width()) if page.hasHeightForWidth() else page.sizeHint().height()
            self.wavelength_tabs.setFixedHeight(height + self.wavelength_tabs.tabBar().sizeHint().height() + 1)
        self.layout().invalidate()
        self.updateGeometry()

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        if event.size().width() != event.oldSize().width():
            self._curve_layout_changed()

    def to_config(self) -> dict:
        curve = {"method": "calibration"}
        if self._switching:
            fits = {}
            for wavelength in self._active_wavelengths:
                try:
                    fits[str(wavelength)] = self.wavelength_curves[wavelength].to_config()
                except ValueError as error:
                    self.wavelength_tabs.setCurrentWidget(self.wavelength_curves[wavelength])
                    raise ValueError(f"Kalibratie bij {wavelength} nm: {error}") from error
            curve["wavelengths"] = fits
        else:
            curve.update(self.curve.to_config())
        if self.concentration_minimum.text().strip() or self.concentration_maximum.text().strip():
            curve["concentration_range"] = [
                _number(self.concentration_minimum, "Minimale concentratie"),
                _number(self.concentration_maximum, "Maximale concentratie"),
            ]
        return {"calculation": curve, "result": _result_config(self.result_label, self.result_unit)}
