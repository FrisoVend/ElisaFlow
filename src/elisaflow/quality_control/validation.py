"""Evaluate acceptance limits against immutable samples and their retained wells."""

from elisaflow.configuration.assay_settings import ValidationLimits, ValidationScope
from elisaflow.models import ResultStatus, Sample, ValidationResult, ValidationStatus
from elisaflow.quality_control.summary import QualityControlSummary, summarize_quality_control


class AssayValidator:
    def __init__(self, scope: ValidationScope):
        self.scope = scope

    def validate(
        self,
        samples: tuple[Sample, ...],
        blank_value: float | None = None,
        *,
        wavelength_nm: int | None = None,
    ) -> tuple[ValidationResult, ...]:
        results = []
        for criterion in self.scope.criteria:
            limits = criterion.limits
            summary = summarize_quality_control(samples, criterion, blank_value=blank_value)
            qualification = summary.result_qualification
            value, status, message = self._check_limits(summary, limits)
            results.append(
                ValidationResult(
                    sample_type=criterion.sample_type,
                    sample_name=criterion.sample_name,
                    metric=criterion.metric,
                    value=value,
                    minimum=limits.minimum,
                    maximum=limits.maximum,
                    status=status,
                    message=message,
                    result_status=qualification.status if qualification is not None else None,
                    result_limit=qualification.limit if qualification is not None else None,
                    nominal_concentration=criterion.nominal_concentration,
                    wavelength_nm=wavelength_nm,
                    cv_percent=summary.cv_percent if value is not None else None,
                    response_scale=limits.response_scale,
                )
            )
        return tuple(results)

    @staticmethod
    def _check_limits(
        summary: QualityControlSummary, limits: ValidationLimits,
    ) -> tuple[float | None, ValidationStatus, str]:
        if limits.minimum is None and limits.maximum is None:
            return None, ValidationStatus.NOT_EVALUATED, "Geen validatiegrenzen ingesteld."
        if summary.mean is None:
            qualification = summary.result_qualification
            result_status = qualification.status if qualification is not None else None
            result_limit = qualification.limit if qualification is not None else None
            if result_limit is not None:
                if result_status == ResultStatus.BELOW_RANGE and limits.minimum is not None and result_limit <= limits.minimum:
                    return (
                        None, ValidationStatus.FAILED,
                        f"Het resultaat is kleiner dan {result_limit:g} en ligt onder het minimum {limits.minimum:g}.",
                    )
                if result_status == ResultStatus.ABOVE_RANGE and limits.maximum is not None and result_limit >= limits.maximum:
                    return (
                        None, ValidationStatus.FAILED,
                        f"Het resultaat is groter dan {result_limit:g} en ligt boven het maximum {limits.maximum:g}.",
                    )
            return None, ValidationStatus.NOT_EVALUATED, summary.message
        value = summary.mean
        if limits.minimum is not None and value < limits.minimum:
            return value, ValidationStatus.FAILED, f"Het gemiddelde ligt onder het minimum {limits.minimum:g}."
        if limits.maximum is not None and value > limits.maximum:
            return value, ValidationStatus.FAILED, f"Het gemiddelde ligt boven het maximum {limits.maximum:g}."
        return value, ValidationStatus.PASSED, "Het gemiddelde ligt binnen de ingestelde grenzen."
