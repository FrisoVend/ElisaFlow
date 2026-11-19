"""Typed, immutable assay configuration parsed at the JSON boundary."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import asdict, dataclass, field, fields, replace
from math import isfinite
from types import MappingProxyType
from typing import TYPE_CHECKING, ClassVar

from elisaflow.calibration.curves import CURVE_MODELS, CalibrationModelDefinition, CalibrationParameters
from elisaflow.configuration.validation import finite_number, validate_config
from elisaflow.formulas.expressions import Formula
from elisaflow.models import CalculationMethod, ResponseScale, SampleType, ValidationMetric
from elisaflow.sample_selection import select_samples

if TYPE_CHECKING:
    from elisaflow.models import Sample


def _nonempty_string(value: object, field_name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"'{field_name}' moet een niet-lege tekst zijn.")


def _validate_formula(formula: Formula | None) -> None:
    if formula is not None and not isinstance(formula, Formula):
        raise ValueError("Voorbewerking moet een Formula-object zijn.")


def _standard_concentration(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (str, int, float)):
        raise ValueError("Standaardsleutels moeten eindige, niet-negatieve concentraties zijn.")
    try:
        concentration = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError("Standaardsleutels moeten eindige, niet-negatieve concentraties zijn.") from error
    if not isfinite(concentration) or concentration < 0:
        raise ValueError("Standaardsleutels moeten eindige, niet-negatieve concentraties zijn.")
    return concentration


def _freeze_standards[T](values: Mapping, prepare_value: Callable[[object], T]) -> Mapping[float, T]:
    standards = {}
    for key, value in values.items():
        concentration = _standard_concentration(key)
        if concentration in standards:
            raise ValueError(f"Dubbele standaardconcentratie: {concentration:g}.")
        standards[concentration] = prepare_value(value)
    return MappingProxyType(standards)


def _model_definition(name: object) -> CalibrationModelDefinition:
    if not isinstance(name, str) or name not in CURVE_MODELS:
        available = ", ".join(sorted(CURVE_MODELS))
        raise ValueError(f"Niet-ondersteund kalibratiemodel {name!r}. Beschikbare modellen: {available}.")
    return CURVE_MODELS[name]


def _freeze_wavelengths(values: Mapping, value_type: type, label: str) -> Mapping[int, object]:
    if not isinstance(values, Mapping):
        raise ValueError(f"{label}: golflengten moeten een object met golflengten als sleutels zijn.")
    wavelengths = {}
    for key, value in values.items():
        if isinstance(key, bool) or not isinstance(key, (str, int)):
            raise ValueError(f"{label}: golflengten moeten positieve gehele getallen zijn.")
        try:
            wavelength = int(key)
        except ValueError as error:
            raise ValueError(f"{label}: golflengten moeten positieve gehele getallen zijn.") from error
        if wavelength <= 0 or wavelength in wavelengths:
            raise ValueError(f"{label}: golflengten moeten unieke positieve gehele getallen zijn.")
        if not isinstance(value, value_type):
            raise ValueError(f"{label}: golflengten moeten {value_type.__name__}-objecten bevatten.")
        wavelengths[wavelength] = value
    return MappingProxyType(wavelengths)


def _validate_concentration_range(settings: CalibrationSettings | WavelengthCalibrationSettings) -> None:
    if settings.concentration_range is not None:
        value = settings.concentration_range
        if not isinstance(value, (list, tuple)) or len(value) != 2:
            raise ValueError("Het concentratiebereik moet een ondergrens en bovengrens bevatten.")
        lower = finite_number(value[0], "Ondergrens concentratiebereik")
        upper = finite_number(value[1], "Bovengrens concentratiebereik")
        if lower < 0 or upper <= lower:
            raise ValueError("Het concentratiebereik moet niet-negatief zijn en oplopen van ondergrens naar bovengrens.")
        object.__setattr__(settings, "concentration_range", (lower, upper))


@dataclass(frozen=True)
class ValidationLimits:
    minimum: float | None = None
    maximum: float | None = None
    response_scale: ResponseScale | None = None

    def __post_init__(self) -> None:
        if self.response_scale is not None:
            object.__setattr__(self, "response_scale", ResponseScale.parse(self.response_scale))
        for name in ("minimum", "maximum"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, finite_number(value, name))
        if self.minimum is not None and self.maximum is not None and self.minimum > self.maximum:
            raise ValueError("Het validatieminimum mag niet groter zijn dan het maximum.")

    @classmethod
    def from_dict(cls, config: object) -> ValidationLimits:
        values = validate_config(config, {"min", "max", "response_scale"}, "Validatiegrenzen")
        if "response_scale" in values:
            ResponseScale.parse(values["response_scale"])
        return cls(values.get("min"), values.get("max"), values.get("response_scale"))

    def to_dict(self) -> dict:
        return {
            name: value for name, value in (
                ("min", self.minimum), ("max", self.maximum), ("response_scale", self.response_scale.value if self.response_scale is not None else None),
            )
            if value is not None
        }


@dataclass(frozen=True)
class ValidationCriterion:
    """One valid acceptance criterion, including its sample identity and response scale."""

    sample_type: SampleType
    sample_name: str
    metric: ValidationMetric
    limits: ValidationLimits
    nominal_concentration: float | None = None

    def __post_init__(self) -> None:
        if self.sample_type not in (SampleType.BLANK, SampleType.CONTROL, SampleType.STANDARD):
            raise ValueError("Validatiecriteria horen bij een blanco, controle of standaard.")
        _nonempty_string(self.sample_name, "Naam validatiemonster")
        if not isinstance(self.metric, ValidationMetric) or not isinstance(self.limits, ValidationLimits):
            raise ValueError("Validatiecriteria vereisen een ValidationMetric en ValidationLimits.")
        if self.sample_type == SampleType.STANDARD:
            object.__setattr__(self, "nominal_concentration", _standard_concentration(self.nominal_concentration))
        elif self.nominal_concentration is not None:
            raise ValueError("Alleen een standaardcriterium heeft een nominale concentratie.")
        scale = self.limits.response_scale
        if self.metric == ValidationMetric.EXTINCTION:
            if scale is None:
                raise ValueError("Extinctiecriteria vereisen een expliciete responsschaal (response_scale).")
            if self.sample_type == SampleType.BLANK and scale != ResponseScale.RAW:
                raise ValueError("Blanco-extinctiecriteria vereisen response_scale 'raw', vóór blancoaftrek.")
        elif scale is not None:
            raise ValueError("response_scale is alleen van toepassing op extinctiecriteria.")

    @property
    def uses_blank_value(self) -> bool:
        return self.sample_type == SampleType.BLANK and self.metric == ValidationMetric.EXTINCTION

    def select(self, samples: Iterable[Sample]) -> tuple[Sample, ...]:
        return select_samples(
            samples, self.sample_type,
            name=self.sample_name if self.sample_type == SampleType.CONTROL else None,
            concentration=self.nominal_concentration,
        )


@dataclass(frozen=True)
class ReplicateFilterSettings:
    """Filter extinctions on the selected scale against their original mean."""

    max_deviation_percent: float
    max_exclusions: int = 1
    response_scale: ResponseScale = field(kw_only=True)

    def __post_init__(self) -> None:
        object.__setattr__(self, "response_scale", ResponseScale.parse(self.response_scale))
        maximum = finite_number(self.max_deviation_percent, "max_deviation_percent")
        if maximum < 0:
            raise ValueError("Het maximale afwijkingspercentage (max_deviation_percent) mag niet negatief zijn.")
        object.__setattr__(self, "max_deviation_percent", maximum)
        if isinstance(self.max_exclusions, bool) or not isinstance(self.max_exclusions, int) or self.max_exclusions < 0:
            raise ValueError("Het maximale aantal uitsluitingen (max_exclusions) moet een niet-negatief geheel getal zijn.")

    @classmethod
    def from_dict(cls, config: object) -> ReplicateFilterSettings:
        values = validate_config(config, {field.name for field in fields(cls)}, "Replicaatfilterinstellingen")
        if "max_deviation_percent" not in values:
            raise ValueError("Replicaatfilterinstellingen vereisen een maximaal afwijkingspercentage (max_deviation_percent).")
        ResponseScale.parse(values.get("response_scale"))
        return cls(**values)

    def to_dict(self) -> dict:
        return {
            "max_deviation_percent": self.max_deviation_percent,
            "max_exclusions": self.max_exclusions,
            "response_scale": self.response_scale.value,
        }


@dataclass(frozen=True)
class PatientSettings:
    replicate_filter: ReplicateFilterSettings | None = None

    def __post_init__(self) -> None:
        if self.replicate_filter is not None and not isinstance(self.replicate_filter, ReplicateFilterSettings):
            raise ValueError("Het replicaatfilter moet een ReplicateFilterSettings-object zijn.")

    @classmethod
    def from_dict(cls, config: object) -> PatientSettings:
        values = validate_config(config, {"replicate_filter"}, "Patiëntinstellingen")
        return cls(
            replicate_filter=ReplicateFilterSettings.from_dict(values["replicate_filter"]) if "replicate_filter" in values else None,
        )

    def to_dict(self) -> dict:
        return {"replicate_filter": self.replicate_filter.to_dict()} if self.replicate_filter is not None else {}


@dataclass(frozen=True)
class SampleSettings:
    lot_number: str | None = None
    replicate_filter: ReplicateFilterSettings | None = None

    def __post_init__(self) -> None:
        if self.lot_number is not None:
            _nonempty_string(self.lot_number, "lot_number")
        if self.replicate_filter is not None and not isinstance(self.replicate_filter, ReplicateFilterSettings):
            raise ValueError("Het replicaatfilter moet een ReplicateFilterSettings-object zijn.")

    @classmethod
    def from_dict(cls, config: object) -> SampleSettings:
        values = validate_config(config, {"lot_number", "replicate_filter"}, "Monsterinstellingen")
        return cls(
            lot_number=values.get("lot_number"),
            replicate_filter=ReplicateFilterSettings.from_dict(values["replicate_filter"]) if "replicate_filter" in values else None,
        )

    def to_dict(self) -> dict:
        config = {}
        if self.lot_number is not None:
            config["lot_number"] = self.lot_number
        if self.replicate_filter is not None:
            config["replicate_filter"] = self.replicate_filter.to_dict()
        return config


def _sample_settings(value: object, label: str) -> SampleSettings:
    if not isinstance(value, SampleSettings):
        raise ValueError(f"{label} moeten SampleSettings-objecten bevatten.")
    return value


@dataclass(frozen=True)
class FormulaOverrideSettings:
    """An independent control input scale and optional replacement expression."""

    expression: str | None = None
    response_scale: ResponseScale = field(kw_only=True)

    def __post_init__(self) -> None:
        if self.expression is not None:
            _nonempty_string(self.expression, "expression")
        object.__setattr__(self, "response_scale", ResponseScale.parse(self.response_scale))

    @classmethod
    def from_dict(cls, config: object) -> FormulaOverrideSettings:
        values = validate_config(config, {"expression", "response_scale"}, "Controleformule")
        return cls(values.get("expression"), response_scale=values.get("response_scale"))

    def to_dict(self) -> dict:
        config = {"response_scale": self.response_scale.value}
        if self.expression is not None:
            config["expression"] = self.expression
        return config


@dataclass(frozen=True)
class FormulaSettings:
    formula: Formula
    control_overrides: Mapping[str, FormulaOverrideSettings] = field(default_factory=dict)
    control_results: Mapping[str, ResultSettings] = field(default_factory=dict)
    response_scale: ResponseScale = field(kw_only=True)
    method: ClassVar[CalculationMethod] = CalculationMethod.FORMULA
    _overrides: Mapping[str, Formula] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "response_scale", ResponseScale.parse(self.response_scale))
        if not isinstance(self.formula, Formula) or not isinstance(self.control_overrides, Mapping):
            raise ValueError("Formuleberekening vereist een Formula-object en een koppeling van afwijkende controleformules.")
        overrides = {}
        for name, settings in self.control_overrides.items():
            _nonempty_string(name, "Naam afwijkende controleformule")
            if not isinstance(settings, FormulaOverrideSettings):
                raise ValueError("Afwijkende controleformules moeten FormulaOverrideSettings-objecten bevatten.")
            overrides[name] = Formula(settings.expression if settings.expression is not None else self.formula.expression, self.formula.references)
        object.__setattr__(self, "control_overrides", MappingProxyType(dict(self.control_overrides)))
        object.__setattr__(self, "_overrides", MappingProxyType(overrides))
        if not isinstance(self.control_results, Mapping):
            raise ValueError("Controleresultaten moeten controlenamen aan resultaatinstellingen koppelen.")
        for name, result in self.control_results.items():
            _nonempty_string(name, "Controlenaam")
            if not isinstance(result, ResultSettings):
                raise ValueError("Controleresultaten moeten ResultSettings-objecten bevatten.")
        object.__setattr__(self, "control_results", MappingProxyType(dict(self.control_results)))

    @property
    def overrides(self) -> Mapping[str, Formula]:
        return self._overrides

    def response_scale_for_control(self, name: str) -> ResponseScale:
        return self.control_overrides[name].response_scale if name in self.control_overrides else self.response_scale

    @classmethod
    def from_dict(cls, config: object) -> FormulaSettings:
        values = validate_config(
            config, {"method", "expression", "references", "control_overrides", "control_results", "response_scale"}, "Formuleberekening",
        )
        if values.get("method") != "formula":
            raise ValueError("De formuleberekeningsmethode moet 'formula' zijn.")
        results = values.get("control_results", {})
        if not isinstance(results, dict):
            raise ValueError("Controleresultaten moeten een object met controlenamen als sleutels zijn.")
        overrides = values.get("control_overrides", {})
        if not isinstance(overrides, dict):
            raise ValueError("Afwijkende controleformules moeten een object met controlenamen als sleutels zijn.")
        return cls(
            Formula(values.get("expression"), values.get("references", {})),
            {name: FormulaOverrideSettings.from_dict(override) for name, override in overrides.items()},
            {name: ResultSettings.from_dict(result) for name, result in results.items()},
            response_scale=values.get("response_scale"),
        )

    def to_dict(self) -> dict:
        config = {"method": "formula", **self.formula.to_dict(), "response_scale": self.response_scale.value}
        if self.control_overrides:
            config["control_overrides"] = {name: override.to_dict() for name, override in self.control_overrides.items()}
        if self.control_results:
            config["control_results"] = {name: result.to_dict() for name, result in self.control_results.items()}
        return config


@dataclass(frozen=True)
class CalibrationSettings:
    """Validated calibration settings with typed parameters; errors are fractions."""

    method: ClassVar[CalculationMethod] = CalculationMethod.CALIBRATION

    fit_model: str
    parameters: CalibrationParameters | None = None
    concentration_range: tuple[float, float] | None = None
    weighting_exponent: float = 0.0
    preprocessing: Formula | None = None
    response_scale: ResponseScale = field(kw_only=True)

    def __post_init__(self) -> None:
        _validate_formula(self.preprocessing)
        if self.preprocessing is not None and self.parameters is None:
            raise ValueError("Voorbewerking is alleen beschikbaar voor curven met vaste parameters.")
        _validate_concentration_range(self)
        definition = _model_definition(self.fit_model)
        object.__setattr__(self, "response_scale", ResponseScale.parse(self.response_scale))
        object.__setattr__(self, "weighting_exponent", finite_number(self.weighting_exponent, "weighting_exponent"))
        if self.parameters is not None:
            if not isinstance(self.parameters, definition.parameter_type):
                raise ValueError(f"De parameters moeten een {definition.parameter_type.__name__}-object zijn.")
            if self.weighting_exponent != 0:
                raise ValueError("De wegingsexponent is alleen van toepassing bij curvepassing.")
            if self.fit_model in ("quadratic", "cubic") and self.concentration_range is None:
                raise ValueError("Vaste polynomiale kalibratie vereist een concentratiebereik om een monotone tak te selecteren.")

    @classmethod
    def from_dict(cls, config: object) -> CalibrationSettings:
        values = dict(validate_config(config, {field.name for field in fields(cls)} | {"method"}, "Curve-instellingen"))
        if values.pop("method", None) != "calibration":
            raise ValueError("De kalibratiemethode moet 'calibration' zijn.")
        definition = _model_definition(values.get("fit_model"))
        ResponseScale.parse(values.get("response_scale"))
        if values.get("parameters") is not None:
            values["parameters"] = definition.parse_parameters(values["parameters"])
        if "preprocessing" in values:
            values["preprocessing"] = Formula.from_dict(values["preprocessing"])
        return cls(**values)

    def to_dict(self) -> dict:
        config = {
            "method": "calibration", "fit_model": self.fit_model,
            "weighting_exponent": self.weighting_exponent, "response_scale": self.response_scale.value,
        }
        if self.parameters is not None:
            config["parameters"] = asdict(self.parameters)
        if self.concentration_range is not None:
            config["concentration_range"] = list(self.concentration_range)
        if self.preprocessing:
            config["preprocessing"] = self.preprocessing.to_dict()
        return config


@dataclass(frozen=True)
class WavelengthCalibrationSettings:
    """Independent curves with shared limits and optional fixed-curve preprocessing."""

    wavelengths: Mapping[int, CalibrationSettings]
    concentration_range: tuple[float, float] | None = None
    method: ClassVar[CalculationMethod] = CalculationMethod.CALIBRATION

    def __post_init__(self) -> None:
        _validate_concentration_range(self)
        wavelengths = _freeze_wavelengths(self.wavelengths, CalibrationSettings, "Kalibratie")
        if len(wavelengths) != 2:
            raise ValueError("Golflengtewisseling vereist kalibratie-instellingen voor precies twee golflengten.")
        for curve in wavelengths.values():
            if curve.concentration_range not in (None, self.concentration_range):
                raise ValueError("Het concentratiebereik moet gezamenlijk in calculation worden ingesteld.")
        object.__setattr__(self, "wavelengths", MappingProxyType({
            key: replace(curve, concentration_range=self.concentration_range) for key, curve in wavelengths.items()
        }))

    def for_wavelength(self, wavelength_nm: int | None) -> CalibrationSettings:
        if wavelength_nm not in self.wavelengths:
            raise ValueError(f"Geen kalibratie-instellingen beschikbaar voor {wavelength_nm} nm.")
        return self.wavelengths[wavelength_nm]

    @classmethod
    def from_dict(cls, config: object) -> WavelengthCalibrationSettings:
        values = dict(validate_config(config, {field.name for field in fields(cls)} | {"method"}, "Kalibratie-instellingen per golflengte"))
        if values.pop("method", None) != "calibration":
            raise ValueError("De kalibratiemethode moet 'calibration' zijn.")
        curves = values.get("wavelengths")
        if not isinstance(curves, dict):
            raise ValueError("Kalibratiegolflengten moeten een object met golflengten als sleutels zijn.")
        values["wavelengths"] = {}
        for key, curve in curves.items():
            fit = validate_config(curve, {"fit_model", "weighting_exponent", "parameters", "response_scale", "preprocessing"}, f"Kalibratie bij {key} nm")
            values["wavelengths"][key] = CalibrationSettings.from_dict({
                "method": "calibration", **fit, "concentration_range": values.get("concentration_range"),
            })
        return cls(**values)

    def to_dict(self) -> dict:
        config = {
            "method": "calibration",
            "wavelengths": {
                str(wavelength): {key: value for key, value in curve.to_dict().items() if key not in {"method", "concentration_range"}}
                for wavelength, curve in self.wavelengths.items()
            },
        }
        if self.concentration_range is not None:
            config["concentration_range"] = list(self.concentration_range)
        return config


@dataclass(frozen=True)
class ResultSettings:
    label: str = "Resultaat"
    unit: str | None = None

    def __post_init__(self) -> None:
        _nonempty_string(self.label, "Resultaatlabel")
        object.__setattr__(self, "label", self.label.strip())
        if self.unit is not None:
            _nonempty_string(self.unit, "Resultaateenheid")
            object.__setattr__(self, "unit", self.unit.strip())

    @classmethod
    def from_dict(cls, config: object) -> ResultSettings:
        return cls(**validate_config(config, {"label", "unit"}, "Resultaatinstellingen"))

    def to_dict(self) -> dict:
        return {"label": self.label, **({"unit": self.unit} if self.unit else {})}


@dataclass(frozen=True)
class ReportSettings:
    """Display precision and warning thresholds, these do not affect analysis."""

    MAX_DECIMALS: ClassVar[int] = 6

    result_decimals: int = 1
    extinction_decimals: int = 3
    calibration_decimals: int = 1
    percentage_decimals: int = 1
    cv_warning_percent: float = 20
    recovery_warning_percent: float = 20

    show_result_cv: bool = True

    def __post_init__(self) -> None:
        if not isinstance(self.show_result_cv, bool):
            raise ValueError("CV van resultaten tonen moet waar of onwaar zijn.")
        for name in ("result_decimals", "extinction_decimals", "calibration_decimals", "percentage_decimals"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= self.MAX_DECIMALS:
                raise ValueError(f"'{name}' moet een geheel getal tussen 0 en {self.MAX_DECIMALS} zijn.")
        for name in ("cv_warning_percent", "recovery_warning_percent"):
            value = finite_number(getattr(self, name), name)
            if value < 0:
                raise ValueError(f"'{name}' mag niet negatief zijn.")
            object.__setattr__(self, name, value)

    @classmethod
    def from_dict(cls, config: object) -> ReportSettings:
        values = validate_config(config, {item.name for item in fields(cls)}, "Rapportinstellingen")
        return cls(**values)

    def to_dict(self) -> dict:
        return asdict(self)


ValidationRules = Mapping[ValidationMetric, ValidationLimits]


def _freeze_rules(rules: ValidationRules) -> ValidationRules:
    if not isinstance(rules, Mapping) or any(
        not isinstance(metric, ValidationMetric) or not isinstance(limits, ValidationLimits) for metric, limits in rules.items()
    ):
        raise ValueError("Validatiecriteria moeten ValidationMetric aan ValidationLimits koppelen.")
    return MappingProxyType(dict(rules))


def _parse_rules(config: object) -> ValidationRules:
    values = validate_config(config, {metric.value for metric in ValidationMetric}, "Validatiecriteria")
    return {ValidationMetric(metric): ValidationLimits.from_dict(limits) for metric, limits in values.items()}


def _rules_config(rules: ValidationRules) -> dict:
    return {metric.value: limits.to_dict() for metric, limits in rules.items()}


@dataclass(frozen=True)
class ValidationScope:
    """Acceptance criteria, independent of sample processing metadata."""

    controls: Mapping[str, ValidationRules] = field(default_factory=dict)
    blank: ValidationRules = field(default_factory=dict)
    standards: Mapping[float, ValidationRules] = field(default_factory=dict)
    criteria: tuple[ValidationCriterion, ...] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.controls, Mapping) or not isinstance(self.standards, Mapping):
            raise ValueError("Validatiecontroles en -standaarden moeten koppelingen zijn.")
        for name in self.controls:
            _nonempty_string(name, "Controlenaam")
        object.__setattr__(self, "controls", MappingProxyType({name: _freeze_rules(rules) for name, rules in self.controls.items()}))
        object.__setattr__(self, "blank", _freeze_rules(self.blank))
        object.__setattr__(self, "standards", _freeze_standards(self.standards, _freeze_rules))
        criteria = [
            ValidationCriterion(SampleType.CONTROL, name, metric, limits)
            for name, rules in self.controls.items() for metric, limits in rules.items()
        ]
        criteria.extend(ValidationCriterion(SampleType.BLANK, "Blanco", metric, limits) for metric, limits in self.blank.items())
        criteria.extend(
            ValidationCriterion(SampleType.STANDARD, f"Standaard {concentration:g}", metric, limits, concentration)
            for concentration, rules in self.standards.items() for metric, limits in rules.items()
        )
        object.__setattr__(self, "criteria", tuple(criteria))

    def has_metric(self, metric: ValidationMetric) -> bool:
        return any(criterion.metric == metric for criterion in self.criteria)

    @classmethod
    def from_dict(cls, config: object) -> ValidationScope:
        values = validate_config(config, {"controls", "blank", "standards"}, "Validatiebereik")
        for name in ("controls", "standards"):
            if not isinstance(values.get(name, {}), dict):
                raise ValueError(f"Validatie '{name}' moet een object zijn.")
        return cls(
            controls={name: _parse_rules(rules) for name, rules in values.get("controls", {}).items()},
            blank=_parse_rules(values.get("blank", {})),
            standards={key: _parse_rules(rules) for key, rules in values.get("standards", {}).items()},
        )

    def to_dict(self) -> dict:
        config = {}
        if self.controls:
            config["controls"] = {name: _rules_config(rules) for name, rules in self.controls.items()}
        if self.blank:
            config["blank"] = _rules_config(self.blank)
        if self.standards:
            config["standards"] = {str(value): _rules_config(rules) for value, rules in self.standards.items()}
        return config


@dataclass(frozen=True)
class SaturationSettings:
    """Convert reference-corrected, unblanked Ext values using fixed coefficients."""

    mode: str = "exclude"
    source_wavelength_nm: int | None = None
    slope: float | None = None
    intercept: float | None = None

    def __post_init__(self) -> None:
        if self.mode not in ("exclude", "calculate"):
            raise ValueError("Verzadiging vereist modus 'exclude' of 'calculate'.")
        if self.mode == "exclude":
            if any(value is not None for value in (self.source_wavelength_nm, self.slope, self.intercept)):
                raise ValueError("Uitsluiten van verzadiging gebruikt geen brongolflengte, helling of intercept.")
            return
        source = self.source_wavelength_nm
        if isinstance(source, bool) or not isinstance(source, int) or source <= 0:
            raise ValueError("source_wavelength_nm moet een positief geheel getal zijn.")
        for name in ("slope", "intercept"):
            object.__setattr__(self, name, finite_number(getattr(self, name), name))
        if self.slope <= 0:
            raise ValueError("De helling (slope) van de verzadigingsconversie moet positief zijn.")

    @classmethod
    def from_dict(cls, config: object) -> SaturationSettings:
        return cls(**validate_config(config, {item.name for item in fields(cls)}, "Verzadigingsinstellingen"))

    def to_dict(self) -> dict:
        return asdict(self) if self.mode == "calculate" else {"mode": "exclude"}


@dataclass(frozen=True)
class MeasurementSettings:
    """Select the measurement input independently of the calibration model."""

    mode: str = "single"
    preferred_wavelength_nm: int | None = None
    reference_wavelength_nm: int | None = None
    cutoff_concentration: float | None = None
    saturation: Mapping[int, SaturationSettings] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.mode not in ("single", "switch_on_concentration"):
            raise ValueError("De meetmethode moet 'single' of 'switch_on_concentration' zijn.")
        wavelength = self.preferred_wavelength_nm
        for name in ("preferred_wavelength_nm", "reference_wavelength_nm"):
            value = getattr(self, name)
            if value is not None and (isinstance(value, bool) or not isinstance(value, int) or value <= 0):
                raise ValueError(f"{name} moet een positief geheel getal zijn.")
        if self.switching != (wavelength is not None):
            raise ValueError("Stel alleen bij meetmethode 'switch_on_concentration' een voorkeursgolflengte (preferred_wavelength_nm) in.")
        if self.switching:
            cutoff = finite_number(self.cutoff_concentration, "cutoff_concentration")
            if cutoff <= 0:
                raise ValueError("De omschakelconcentratie (cutoff_concentration) moet positief zijn.")
            object.__setattr__(self, "cutoff_concentration", cutoff)
        elif self.cutoff_concentration is not None:
            raise ValueError("cutoff_concentration vereist meetmethode 'switch_on_concentration'.")
        object.__setattr__(self, "saturation", _freeze_wavelengths(self.saturation, SaturationSettings, "Verzadiging"))
        for target, rule in self.saturation.items():
            if rule.source_wavelength_nm == target:
                raise ValueError("De brongolflengte voor verzadigingsconversie moet verschillen van de doelgolflengte.")
            if rule.mode == "calculate" and self.reference_wavelength_nm in (target, rule.source_wavelength_nm):
                raise ValueError("De referentiegolflengte mag geen bron of doel zijn voor verzadigingsconversie.")
        if wavelength is not None and wavelength == self.reference_wavelength_nm:
            raise ValueError("De voorkeurs- en referentiegolflengte moeten verschillen.")

    def analysis_wavelengths(self, wavelengths: tuple[int, ...]) -> tuple[int, ...]:
        """Validate selections against the immutable raw wavelengths from SKAX."""
        reference = self.reference_wavelength_nm
        if reference is not None and reference not in wavelengths:
            raise ValueError(f"Referentiegolflengte {reference} nm ontbreekt in het SKAX-protocol.")
        for target, rule in self.saturation.items():
            if target not in wavelengths or rule.mode == "calculate" and rule.source_wavelength_nm not in wavelengths:
                raise ValueError("De doel- en brongolflengten voor verzadiging moeten in het SKAX-protocol voorkomen.")
        analysis = tuple(value for value in wavelengths if value != reference)
        expected = 2 if self.switching else 1
        if len(analysis) != expected:
            raise ValueError(f"Meetmethode {self.mode!r} vereist precies {expected} analysegolflengte(n), plus een eventuele referentie, in SKAX.")
        if self.switching and self.preferred_wavelength_nm not in analysis:
            raise ValueError(f"Voorkeursgolflengte {self.preferred_wavelength_nm} nm ontbreekt in de SKAX-analysegolflengten.")
        return analysis

    def wavelength_label(self, wavelength_nm: int) -> str:
        reference = self.reference_wavelength_nm
        return f"{wavelength_nm} − {reference} nm" if reference is not None else f"{wavelength_nm} nm"

    @property
    def switching(self) -> bool:
        return self.mode == "switch_on_concentration"

    @classmethod
    def from_dict(cls, config: object) -> MeasurementSettings:
        values = dict(validate_config(config, {item.name for item in fields(cls)}, "Meetinstellingen"))
        if "saturation" in values:
            rules = values["saturation"]
            if not isinstance(rules, dict):
                raise ValueError("Verzadiging (saturation) moet een object met golflengten als sleutels zijn.")
            values["saturation"] = {key: SaturationSettings.from_dict(rule) for key, rule in rules.items()}
        return cls(**values)

    def to_dict(self) -> dict:
        config = {"mode": self.mode}
        if self.preferred_wavelength_nm is not None:
            config["preferred_wavelength_nm"] = self.preferred_wavelength_nm
        if self.reference_wavelength_nm is not None:
            config["reference_wavelength_nm"] = self.reference_wavelength_nm
        if self.cutoff_concentration is not None:
            config["cutoff_concentration"] = self.cutoff_concentration
        if self.saturation:
            config["saturation"] = {str(wavelength): rule.to_dict() for wavelength, rule in self.saturation.items()}
        return config


@dataclass(frozen=True)
class ValidationPlan:
    result: ValidationScope = field(default_factory=ValidationScope)
    wavelengths: Mapping[int, ValidationScope] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.result, ValidationScope):
            raise ValueError("Een validatieplan vereist criteria voor het resultaat en de golflengten.")
        object.__setattr__(self, "wavelengths", _freeze_wavelengths(self.wavelengths, ValidationScope, "Validatie"))
        for scope in self.wavelengths.values():
            if scope.has_metric(ValidationMetric.RESULT):
                raise ValueError("Resultaatcriteria moeten in validation.result worden ingesteld, niet in validation.wavelengths.")

    @classmethod
    def from_dict(cls, config: object) -> ValidationPlan:
        values = validate_config(config, {"result", "wavelengths"}, "Validatieplan")
        wavelengths = values.get("wavelengths", {})
        if not isinstance(wavelengths, dict):
            raise ValueError("Validatiegolflengten moeten een object met golflengten als sleutels zijn.")
        return cls(
            ValidationScope.from_dict(values.get("result", {})), {key: ValidationScope.from_dict(scope) for key, scope in wavelengths.items()}
        )

    def to_dict(self) -> dict:
        config = {}
        if result := self.result.to_dict():
            config["result"] = result
        if self.wavelengths:
            config["wavelengths"] = {str(wavelength): scope.to_dict() for wavelength, scope in self.wavelengths.items()}
        return config


@dataclass(frozen=True)
class AssaySettings:
    protocol_name: str
    lot_number: str
    calculation: FormulaSettings | CalibrationSettings | WavelengthCalibrationSettings
    result: ResultSettings = field(default_factory=ResultSettings)
    controls: Mapping[str, SampleSettings] = field(default_factory=dict)
    blank: SampleSettings = field(default_factory=SampleSettings)
    standards: Mapping[float, SampleSettings] = field(default_factory=dict)
    report: ReportSettings = field(default_factory=ReportSettings)
    validation: ValidationPlan = field(default_factory=ValidationPlan)
    measurement: MeasurementSettings = field(default_factory=MeasurementSettings)
    patients: PatientSettings = field(default_factory=PatientSettings)

    def to_dict(self) -> dict:
        config = {
            "lot_number": self.lot_number,
            "measurement": self.measurement.to_dict(),
            "calculation": self.calculation.to_dict(),
            "result": self.result.to_dict(),
            "report": self.report.to_dict(),
        }
        for name, samples in (("controls", self.controls), ("standards", self.standards)):
            if samples:
                config[name] = {str(key): sample.to_dict() for key, sample in samples.items()}
        for name, section in (("patients", self.patients), ("blank", self.blank), ("validation", self.validation)):
            if contents := section.to_dict():
                config[name] = contents
        return config

    def __post_init__(self) -> None:
        _nonempty_string(self.protocol_name, "protocol_name")
        _nonempty_string(self.lot_number, "lot_number")
        for value, kind, name in (
            (self.calculation, (FormulaSettings, CalibrationSettings, WavelengthCalibrationSettings), "calculation"),
            (self.result, ResultSettings, "result"),
            (self.measurement, MeasurementSettings, "measurement"),
            (self.report, ReportSettings, "report"),
            (self.validation, ValidationPlan, "validation"),
            (self.blank, SampleSettings, "blank"),
            (self.patients, PatientSettings, "patients"),
        ):
            if not isinstance(value, kind):
                raise ValueError(f"Ongeldig instellingenobject voor {name}.")
        if self.blank.replicate_filter is not None:
            raise ValueError("Replicaatfiltering van blanco’s wordt niet ondersteund; filtering vindt plaats na blanco-aftrek.")
        if not isinstance(self.controls, Mapping) or not isinstance(self.standards, Mapping):
            raise ValueError("Controles en standaarden moeten koppelingen zijn.")
        for name, sample in self.controls.items():
            _nonempty_string(name, "Controlenaam")
            _sample_settings(sample, "Controles")
        object.__setattr__(self, "controls", MappingProxyType(dict(self.controls)))
        object.__setattr__(self, "standards", _freeze_standards(self.standards, lambda sample: _sample_settings(sample, "Standaarden")))
        if self.measurement.switching:
            if isinstance(self.calculation, FormulaSettings):
                raise ValueError("Golflengtewisseling is alleen van toepassing op kalibratieassays.")
            if not isinstance(self.calculation, WavelengthCalibrationSettings):
                raise ValueError("Golflengtewisseling vereist calculation.wavelengths met curve-instellingen voor beide golflengten.")
            if self.measurement.preferred_wavelength_nm not in self.calculation.wavelengths:
                raise ValueError("De voorkeursgolflengte moet instellingen hebben in calculation.wavelengths.")
            if self.validation.result.has_metric(ValidationMetric.EXTINCTION):
                raise ValueError("Stel bij golflengtewisseling de extinctiecriteria in onder validation.wavelengths.")
        elif isinstance(self.calculation, WavelengthCalibrationSettings):
            raise ValueError("calculation.wavelengths vereist meetmethode 'switch_on_concentration'.")
        elif self.validation.wavelengths:
            raise ValueError("Criteria voor een enkele meting horen in validation.result; criteria per golflengte vereisen golflengtewisseling.")

    def calculation_for_wavelength(self, wavelength_nm: int) -> FormulaSettings | CalibrationSettings:
        if isinstance(self.calculation, WavelengthCalibrationSettings):
            return self.calculation.for_wavelength(wavelength_nm)
        return self.calculation

    @classmethod
    def from_dict(cls, protocol_name: str, config: object) -> AssaySettings:
        values = validate_config(
            config, {"lot_number", "calculation", "result", "patients", "controls", "blank", "standards", "report", "validation", "measurement"}, "Assay-instellingen"
        )
        calculation = values.get("calculation")
        if not isinstance(calculation, dict) or calculation.get("method") not in ("formula", "calibration"):
            raise ValueError("calculation.method moet 'formula' of 'calibration' zijn.")
        calculation_type = (
            FormulaSettings if calculation["method"] == "formula"
            else WavelengthCalibrationSettings if "wavelengths" in calculation
            else CalibrationSettings
        )
        controls, standards = values.get("controls", {}), values.get("standards", {})
        if not isinstance(controls, dict) or not isinstance(standards, dict):
            raise ValueError("Controles en standaarden moeten objecten zijn.")
        return cls(
            protocol_name=protocol_name,
            lot_number=values.get("lot_number"),
            calculation=calculation_type.from_dict(calculation),
            result=ResultSettings.from_dict(values.get("result", {})),
            patients=PatientSettings.from_dict(values.get("patients", {})),
            controls={name: SampleSettings.from_dict(value) for name, value in controls.items()},
            blank=SampleSettings.from_dict(values.get("blank", {})),
            standards={key: SampleSettings.from_dict(value) for key, value in standards.items()},
            report=ReportSettings.from_dict(values.get("report", {})),
            validation=ValidationPlan.from_dict(values.get("validation", {})),
            measurement=MeasurementSettings.from_dict(values.get("measurement", {})),
        )
