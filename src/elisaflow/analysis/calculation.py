"""Calculate retained wells and record the exact input scale with their results."""

from dataclasses import dataclass, replace
from math import isfinite, isinf

import numpy as np

from elisaflow.analysis.sample_updates import clear_result, map_wells
from elisaflow.calibration.curves import CalibratedValue, CalibrationCurve, create_calibration_curve
from elisaflow.calibration.errors import CalibrationError
from elisaflow.configuration.assay_settings import CalibrationSettings, FormulaSettings
from elisaflow.formulas.expressions import Formula
from elisaflow.models import CalculationDetails, ExtinctionScale, ResponseScale, ResultStatus, Sample, SampleType, WellMeasurement
from elisaflow.sample_selection import select_samples

NAN = float("nan")


@dataclass(frozen=True)
class CalculationOutcome:
    """Executed calculation: current wells, details and their actual input scale."""

    samples: tuple[Sample, ...]
    details: CalculationDetails
    input_scale: ExtinctionScale | None


def calculate_samples(
    samples: tuple[Sample, ...],
    calculation: FormulaSettings | CalibrationSettings,
    *,
    wavelength_nm: int,
    reference_wavelength_nm: int | None,
    blank_value: float | None,
    allow_fit_failure: bool,
) -> CalculationOutcome:
    """Calculate after blank correction and replicate filtering, before QC.

    Control references use retained means before overrides or preprocessing.
    Only recognized fitted-curve failures may become unavailable outcomes;
    invalid configuration and unexpected implementation errors propagate.
    """
    return _SampleCalculation(
        calculation, wavelength_nm=wavelength_nm, reference_wavelength_nm=reference_wavelength_nm,
        blank_value=blank_value, allow_fit_failure=allow_fit_failure,
    ).calculate(samples)


class _SampleCalculation:
    def __init__(
        self, calculation: FormulaSettings | CalibrationSettings, *, wavelength_nm: int,
        reference_wavelength_nm: int | None, blank_value: float | None, allow_fit_failure: bool,
    ):
        self.calculation = calculation
        self.allow_fit_failure = allow_fit_failure
        self.blank_value = blank_value
        self.input_scale = ExtinctionScale.for_calculation(
            calculation,
            wavelength_nm=wavelength_nm,
            reference_wavelength_nm=reference_wavelength_nm,
            blank_value=blank_value,
        )

    def calculate(self, samples: tuple[Sample, ...]) -> CalculationOutcome:
        outcome = self._calculate_calibration(samples) if isinstance(self.calculation, CalibrationSettings) else self._calculate_formula(samples)
        samples = tuple(
            replace(sample, measurements=tuple(map(clear_result, sample.measurements)))
            if sample.sample_type == SampleType.UNKNOWN and any(well.reference_invalid and not well.disabled for well in sample.measurements)
            else sample
            for sample in outcome.samples
        )
        return replace(outcome, samples=samples)

    def _input(self, well: WellMeasurement, response_scale: ResponseScale) -> float:
        return response_scale.value_for(well) if well.included else NAN

    def _evaluate(
        self, samples: tuple[Sample, ...], values: np.ndarray, formula: Formula | None, response_scale: ResponseScale,
    ) -> np.ndarray:
        if formula is None:
            return values
        return np.asarray(formula.evaluate(values, self._control_means(samples, formula, response_scale)))

    def _control_means(self, samples: tuple[Sample, ...], formula: Formula, response_scale: ResponseScale) -> dict[str, float]:
        controls = {}
        for name in formula.control_names:
            readings = [
                self._input(well, response_scale)
                for sample in samples
                if sample.sample_type == SampleType.CONTROL and sample.name == name
                for well in sample.retained
            ]
            if not readings:
                raise ValueError(f"Geen metingen gevonden voor controle {name!r}.")
            if not np.isfinite(readings).all():
                raise ValueError(f"Metingen voor controle {name!r} moeten eindig zijn.")
            controls[name] = float(np.mean(readings))
        return controls

    def _calculate_formula(self, samples: tuple[Sample, ...]) -> CalculationOutcome:
        calculated = []
        for sample in samples:
            scale = (
                ExtinctionScale.for_calculation(
                    self.calculation, wavelength_nm=self.input_scale.wavelength_nm,
                    reference_wavelength_nm=self.input_scale.reference_wavelength_nm,
                    blank_value=self.blank_value, control_name=sample.name,
                ) if sample.sample_type == SampleType.CONTROL else self.input_scale
            )
            response_scale = scale.response_scale
            formula = (
                self.calculation.overrides.get(sample.name, self.calculation.formula)
                if sample.sample_type == SampleType.CONTROL
                else self.calculation.formula
            )
            values = self._evaluate(samples, np.array([self._input(well, response_scale) for well in sample.measurements]), formula, response_scale)
            calculated.append(
                replace(
                    sample,
                    calculation_input_scale=scale,
                    measurements=tuple(
                        replace(well, result=float(value), result_status=ResultStatus.IN_RANGE if isfinite(value) else ResultStatus.UNAVAILABLE)
                        for well, value in zip(sample.measurements, values, strict=True)
                    ),
                )
            )
        return CalculationOutcome(
            tuple(calculated), CalculationDetails(formula_settings=self.calculation), self.input_scale,
        )

    def _calculate_calibration(self, samples: tuple[Sample, ...]) -> CalculationOutcome:
        settings = self.calculation
        readings = [(sample, well) for sample in samples if sample.sample_type != SampleType.BLANK for well in sample.measurements]
        inputs = self._evaluate(
            samples, np.array([self._input(well, self.input_scale.response_scale) for _, well in readings]),
            settings.preprocessing, self.input_scale.response_scale,
        )
        standards = [index for index, (sample, well) in enumerate(readings) if sample.sample_type == SampleType.STANDARD and well.included]
        if settings.parameters is None:
            if any(readings[index][0].concentration is None for index in standards):
                raise ValueError("Curvepassing vereist een nominale concentratie voor elke standaard.")
            try:
                curve = create_calibration_curve(
                    settings,
                    concentrations=[readings[index][0].concentration for index in standards],
                    extinction_values=inputs[standards],
                )
            except CalibrationError as error:
                if not self.allow_fit_failure:
                    raise
                return CalculationOutcome(samples, CalculationDetails(calibration_settings=settings, calibration_error=str(error)), None)
        else:
            curve = create_calibration_curve(settings)
        values = curve.calculate_results(inputs)
        saturation = self._saturation_evidence(samples)
        back_calculated = dict(zip(standards, curve.diagnostics.back_calculated, strict=True)) if curve.was_fitted else {}
        calculated = {
            well.key: self._calculated_well(
                sample, well, self._saturated_result(sample, well, curve, saturation) or value, float(response), back_calculated.get(index),
            )
            for index, ((sample, well), value, response) in enumerate(zip(readings, values, inputs, strict=True))
        }
        samples = map_wells(samples, lambda sample, well: well if sample.sample_type == SampleType.BLANK else calculated[well.key])
        samples = tuple(
            replace(sample, calculation_input_scale=self.input_scale) if sample.sample_type != SampleType.BLANK else sample for sample in samples
        )
        details = CalculationDetails(
            curve=curve,
            calibration_settings=settings,
            curve_r_squared=self._calibration_r_squared(samples, curve) if curve.was_fitted else None,
        )
        return CalculationOutcome(samples, details, self.input_scale)

    def _calibration_r_squared(self, samples: tuple[Sample, ...], curve: CalibrationCurve) -> float | None:
        """Compare standard means on the recorded curve-input scale."""
        pairs = [
            (sample.concentration, self.input_scale.statistics(sample).mean)
            for sample in select_samples(samples, SampleType.STANDARD)
            if sample.concentration is not None
        ]
        pairs = [(nominal, observed) for nominal, observed in pairs if isfinite(nominal) and nominal >= 0 and isfinite(observed)]
        if not pairs:
            return None
        nominal, observed = np.asarray(pairs).T
        predictions = np.asarray(curve.calculate_extinction(nominal))
        total = np.sum((observed - observed.mean()) ** 2)
        if total == 0:
            return None
        value = float(1 - np.sum((observed - predictions) ** 2) / total)
        return value if isfinite(value) else None

    def _saturation_evidence(self, samples: tuple[Sample, ...]) -> tuple[float, float, float] | None:
        """Saturation exceeds reliable raw standard signals at this wavelength."""
        reliable = [
            well.raw_extinction if well.raw_extinction is not None else well.source_extinction
            for sample in samples if sample.sample_type == SampleType.STANDARD
            for well in sample.retained if not well.saturated
        ]
        reliable = [value for value in reliable if isfinite(value)]
        if not reliable:
            return None
        formula = self.calculation.preprocessing
        coefficients = (1.0, 0.0) if formula is None else formula.affine_coefficients(self._control_means(samples, formula, self.calculation.response_scale))
        if coefficients is None or coefficients[0] == 0:
            return None
        return max(reliable), *coefficients

    def _saturated_result(
        self, sample: Sample, well: WellMeasurement, curve: CalibrationCurve, evidence: tuple[float, float, float] | None,
    ) -> CalibratedValue | None:
        if sample.sample_type != SampleType.UNKNOWN or not well.saturation_unresolved or well.disabled or well.reference_invalid or evidence is None:
            return None
        raw_bound, slope, intercept = evidence
        response = raw_bound - (well.reference_extinction or 0)
        response -= self.input_scale.blank_subtracted or 0
        return curve.calculate_bound(slope * response + intercept, above=slope > 0)

    @staticmethod
    def _calculated_well(
        sample: Sample,
        well: WellMeasurement,
        calculated: CalibratedValue,
        response: float,
        back_calculated: float | None,
    ) -> WellMeasurement:
        value = calculated.concentration if well.included else NAN
        qualified = well.included or (
            sample.sample_type == SampleType.UNKNOWN and well.saturation_unresolved and not well.disabled and not well.reference_invalid
            and calculated.status in (ResultStatus.BELOW_RANGE, ResultStatus.ABOVE_RANGE)
        )
        status = calculated.status if qualified else ResultStatus.UNAVAILABLE
        limit = calculated.limit if qualified else NAN
        finite_before_dilution = isfinite(value)
        if sample.sample_type == SampleType.UNKNOWN:
            value *= sample.dilution
            limit *= sample.dilution
        if finite_before_dilution and not isfinite(value):
            raise ValueError("Voor verdunning gecorrigeerde concentraties moeten eindig zijn.")
        if isinf(limit):
            raise ValueError("Voor verdunning gecorrigeerde concentratiegrenzen moeten eindig zijn.")
        return replace(
            well,
            result=value,
            result_status=status,
            result_limit=limit,
            curve_extinction=response,
            back_calculated_concentration=back_calculated,
        )
