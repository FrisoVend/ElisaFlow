"""Edit QC, calculation and report settings for a SKAX protocol."""

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QLineEdit,
    QSizePolicy,
    QWidget,
)

from elisaflow.configuration.assay_settings import ReplicateFilterSettings
from elisaflow.models import ResponseScale
from elisaflow.ui.theme import FormSection, form_layout


def _response_scale_selector(value: ResponseScale | str, *, raw_only: bool = False) -> QComboBox:
    selector = QComboBox()
    for label, scale in (("Vóór blancoaftrek", ResponseScale.RAW), ("Na blancoaftrek", ResponseScale.BLANK_CORRECTED)):
        if not raw_only or scale == ResponseScale.RAW:
            selector.addItem(label, scale.value)
    selector.setCurrentIndex(selector.findData(ResponseScale.parse(value).value))
    selector.setMinimumContentsLength(18)
    selector.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
    selector.setToolTip(
        "Beide schalen zijn na eventuele referentiecorrectie en verzadigingsomrekening. "
        "Zonder blanco zijn beide schalen gelijk."
    )
    return selector


def _number(field: QLineEdit, label: str) -> float:
    try:
        return float(field.text().strip())
    except ValueError as error:
        raise ValueError(f"{label} moet een getal zijn.") from error


def _result_config(label: QLineEdit, unit: QLineEdit) -> dict[str, str]:
    result = {"label": label.text().strip()}
    if value := unit.text().strip():
        result["unit"] = value
    return result


def _range_fields(minimum: object = None, maximum: object = None) -> tuple:
    widget = QWidget()
    widget.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    layout = QHBoxLayout(widget)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(12)
    fields = []
    for value, label in ((minimum, "Minimum"), (maximum, "Maximum")):
        field = QLineEdit("" if value is None else str(value))
        field.setPlaceholderText(label)
        layout.addWidget(field)
        fields.append(field)
    return widget, fields[0], fields[1]


class _ReplicateFilterEditor(FormSection):
    def __init__(self, config: dict | None, *, has_blank: bool) -> None:
        super().__init__("Replicaatfilter (extinctie)")
        self.enabled = QCheckBox("Replicaatfilter inschakelen")
        self.enabled.setChecked(config is not None)
        self.form.addRow(self.enabled)
        parameters = QWidget()
        parameters.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        form = form_layout(parameters)
        values = config or {}
        self.response_scale = _response_scale_selector(values.get("response_scale", ResponseScale.default_for_blank(has_blank)))
        self.response_scale.setToolTip(
            f"{self.response_scale.toolTip()} De keuze geldt voor het replicaatgemiddelde en de procentuele afwijking, vóór formules en kalibratievoorbewerking."
        )
        form.addRow("Responsschaal:", self.response_scale)
        self.fields = {}
        for name, label, default in (
            ("max_deviation_percent", "Maximale extinctieafwijking (%):", None),
            ("max_exclusions", "Maximaal aantal uitsluitingen:", 1),
        ):
            value = values.get(name, default)
            field = QLineEdit("" if value is None else str(value))
            self.fields[name] = field
            form.addRow(label, field)
        parameters.setEnabled(self.enabled.isChecked())
        self.enabled.toggled.connect(parameters.setEnabled)
        self.form.addRow(parameters)

    def to_config(self) -> dict | None:
        if not self.enabled.isChecked():
            return None
        values = {
            "max_deviation_percent": _number(self.fields["max_deviation_percent"], "Extinctieafwijking"),
            "response_scale": self.response_scale.currentData(),
        }
        try:
            values["max_exclusions"] = int(self.fields["max_exclusions"].text().strip())
        except ValueError as error:
            raise ValueError("Het maximale aantal uitsluitingen moet een geheel getal zijn.") from error
        return ReplicateFilterSettings.from_dict(values).to_dict()
