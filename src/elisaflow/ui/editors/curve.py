"""One calibration curve editor, shared by single- and multiple-wavelength assays."""

from dataclasses import fields

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QCheckBox, QComboBox, QLineEdit, QWidget

from elisaflow.calibration.curves import CURVE_MODELS
from elisaflow.models import ResponseScale
from elisaflow.ui.editors.fields import _number, _response_scale_selector
from elisaflow.ui.editors.formula import FormulaEditor
from elisaflow.ui.theme import form_layout


class CurveEditor(QWidget):
    layout_changed = Signal()

    def __init__(self, config: dict, control_names: tuple[str, ...], *, editable: bool) -> None:
        super().__init__()
        self._form = form_layout(self)
        self.fit_model = QComboBox()
        for label, model in (
            ("4pl", "4pl"), ("5pl", "5pl"), ("Kwadratisch", "quadratic"), ("Kubisch", "cubic"),
            ("Machtsfunctie met offset", "power_law"),
        ):
            self.fit_model.addItem(label, model)
        self.fit_model.setCurrentIndex(self.fit_model.findData(config.get("fit_model", "4pl")))
        self._form.addRow("Curvemodel:", self.fit_model)
        self.fixed_parameters = QCheckBox("Vaste parameters gebruiken")
        self.fixed_parameters.setChecked(config.get("fixed_parameters", "parameters" in config))
        self._form.addRow(self.fixed_parameters)
        default_scale = ResponseScale.BLANK_CORRECTED if self.fixed_parameters.isChecked() else ResponseScale.RAW
        self.response_scale = _response_scale_selector(config.get("response_scale", default_scale))
        self.response_scale.setToolTip(
            f"{self.response_scale.toolTip()} De keuze geldt voor standaarden, controles en patiënten, vóór voorbewerking."
        )
        self._form.addRow("Responsschaal:", self.response_scale)
        self.parameters = {}
        names = dict.fromkeys(parameter.name for model in CURVE_MODELS.values() for parameter in fields(model.parameter_type))
        for name in names:
            value = config.get("parameters", {}).get(name, 1.0 if name == "E" else None)
            field = QLineEdit("" if value is None else str(value))
            field.setAccessibleName(f"Parameter {name}")
            self.parameters[name] = field
            self._form.addRow(f"Parameter {name}:", field)
        self.equation = QLineEdit()
        self.equation.setReadOnly(True)
        self._form.addRow("Formule:", self.equation)
        self.weighting_exponent = QLineEdit(str(config.get("weighting_exponent", 0.0)))
        self.weighting_exponent.setAlignment(Qt.AlignmentFlag.AlignRight)
        self.weighting_exponent.setToolTip(
            "Weging = 1/y^p op de gekozen responsschaal; p=0 is ongewogen. "
            "Een exponent ongelijk aan nul vereist positieve invoerwaarden. Decimale exponenten zijn toegestaan."
        )
        self._form.addRow("Wegingsexponent (1/y^p):", self.weighting_exponent)
        self.preprocessing = FormulaEditor(
            config.get("preprocessing"), control_names, sample_name="Monster", optional=True, title="Voorbewerking kalibratie",
        )
        self._form.addRow("Voorbewerking:", self.preprocessing)
        self.fit_model.currentIndexChanged.connect(self._refresh)
        self.fixed_parameters.toggled.connect(self._refresh)
        self._refresh()
        self.setEnabled(editable)

    def _refresh(self) -> None:
        fixed = self.fixed_parameters.isChecked()
        self._form.setRowVisible(self.weighting_exponent, not fixed)
        self._form.setRowVisible(self.preprocessing, fixed)
        if not fixed:
            self.preprocessing.clear()
        model = self.fit_model.currentData()
        selected = {parameter.name for parameter in fields(CURVE_MODELS[model].parameter_type)}
        for name, field in self.parameters.items():
            self._form.setRowVisible(field, fixed and name in selected)
        equations = {"quadratic": "y = A + Bx + Cx²", "cubic": "y = A + Bx + Cx² + Dx³", "power_law": "y = C + A·x^B"}
        self.equation.setText(equations.get(model, ""))
        self._form.setRowVisible(self.equation, fixed and model in equations)
        self._form.invalidate()
        self.updateGeometry()
        self.layout_changed.emit()

    def to_config(self) -> dict:
        config = {"fit_model": self.fit_model.currentData(), "response_scale": self.response_scale.currentData()}
        if self.fixed_parameters.isChecked():
            config["parameters"] = {
                parameter.name: _number(self.parameters[parameter.name], f"Parameter {parameter.name}")
                for parameter in fields(CURVE_MODELS[config["fit_model"]].parameter_type)
            }
            preprocessing = self.preprocessing.to_config()
            if preprocessing is not None:
                config["preprocessing"] = preprocessing
        else:
            config["weighting_exponent"] = _number(self.weighting_exponent, "Wegingsexponent")
        return config

    def draft(self) -> dict:
        """Copy all entered text, including incomplete or currently hidden fields."""
        return {
            "fit_model": self.fit_model.currentData(),
            "fixed_parameters": self.fixed_parameters.isChecked(),
            "response_scale": self.response_scale.currentData(),
            "weighting_exponent": self.weighting_exponent.text(),
            "parameters": {name: field.text() for name, field in self.parameters.items()},
            "preprocessing": self.preprocessing.to_config(),
        }


