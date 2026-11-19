"""Expected calibration failures, distinct from unexpected implementation errors."""


class CalibrationError(Exception):
    """A curve cannot provide results for the supplied standards or limits."""


class CalibrationDataError(CalibrationError, ValueError):
    """Calibration inputs do not support the requested model or range."""


class CalibrationFitError(CalibrationError, RuntimeError):
    """Numerical fitting or configured fit acceptance failed."""
