import unittest

import numpy as np
import pandas as pd
from fixtures import reference_parameters, samples_with_statistics

from elisaflow.calibration.logistic import FivePLModel, FivePLParameters, FourPLModel, FourPLParameters
from elisaflow.calibration.plot import plot_calibration
from elisaflow.calibration.polynomial import CubicModel, CubicParameters, QuadraticModel, QuadraticParameters
from elisaflow.calibration.power_law import PowerLawModel, PowerLawParameters
from elisaflow.formulas.expressions import Formula
from elisaflow.models import Statistics


class CalibrationPlotTests(unittest.TestCase):
    def test_preprocessing_uses_plain_formula_notation_inside_the_legend(self):
        formula = Formula("sample * (0.86 / ref)", {"ref": "QC {1} $x$"})
        figure = plot_calibration(
            samples_with_statistics(self.samples.to_dict("records")), self.model, self.parameters,
            input_statistics=self.input_statistics, model_name="4pl", preprocessing=formula, figsize=(3.8, 2.44),
        )
        figure.canvas.draw()
        axes = figure.axes[0]
        self.assertEqual(axes.get_title(), "")
        legend = axes.get_legend()
        self.assertIn("Voorbewerking: Monster * (0.86 / 'QC {1} $x$')", legend.get_title().get_text())
        self.assertFalse(legend.get_title().get_parse_math())
        self.assertTrue(legend.get_frame_on())
        bounds = legend.get_window_extent(figure.canvas.get_renderer())
        self.assertGreaterEqual(bounds.x0, axes.bbox.x0)
        self.assertLessEqual(bounds.x1, axes.bbox.x1)
        self.assertGreaterEqual(bounds.y0, axes.bbox.y0)
        self.assertLessEqual(bounds.y1, axes.bbox.y1)

    def test_fixed_curve_can_be_plotted_without_standards_and_without_fit_statistics(self):
        figure = plot_calibration(
            (), self.model, self.parameters, input_statistics=(), model_name="4pl", curve_was_fitted=False,
            concentration_range=(0, 500), preprocessing=Formula("sample * 2"),
        )
        figure.canvas.draw()
        axes = figure.axes[0]
        line = axes.lines[0]
        self.assertEqual(line.get_label(), "Vaste curve")
        self.assertEqual(line.get_xdata()[-1], 500)
        self.assertGreater(line.get_xdata()[0], 0)
        self.assertEqual(axes.get_xlabel(), "Concentratie")
        self.assertEqual(len(axes.collections), 0)
        np.testing.assert_allclose(line.get_ydata(), self.model.calculate_extinctions(line.get_xdata(), self.parameters))
        self.assertNotIn("R²", axes.get_legend().get_title().get_text())

    def test_fixed_curve_uses_concentrations_before_patient_dilution_for_its_display_range(self):
        samples = samples_with_statistics([
            {"sample_type": "Unknown", "sample_dilution": 10, "result_mean": value} for value in (100, 1000)
        ])
        figure = plot_calibration(samples, self.model, self.parameters, input_statistics=(), curve_was_fitted=False)
        np.testing.assert_allclose(figure.axes[0].lines[0].get_xdata()[[0, -1]], [10, 100])

    def test_fixed_models_have_a_display_domain_even_without_measurements(self):
        for model, parameters, limits in (
            (FourPLModel(), self.parameters, None),
            (FivePLModel(), FivePLParameters(2.5, -1.3, 2, 0.1, 0.4), None),
            (QuadraticModel(), QuadraticParameters(0.1, 0.02, 0.0001), (0, 100)),
            (CubicModel(), CubicParameters(0.1, 0.02, 0.0001, 0.000001), (0, 100)),
            (PowerLawModel(), PowerLawParameters(0.3, -0.7, 0.1), None),
        ):
            with self.subTest(model=type(model).__name__):
                figure = plot_calibration((), model, parameters, input_statistics=(), curve_was_fitted=False, concentration_range=limits)
                figure.canvas.draw()
                values = figure.axes[0].lines[0].get_xdata()
                self.assertGreater(values[0], 0)
                self.assertGreater(values[-1], values[0])
                self.assertTrue(np.isfinite(figure.axes[0].lines[0].get_ydata()).all())

    def test_explicit_missing_inputs_do_not_fall_back_to_finite_sample_values(self):
        samples = samples_with_statistics(self.samples.to_dict("records"))
        missing = (Statistics(),) * 4
        with self.assertRaisesRegex(ValueError, "positieve concentraties"):
            plot_calibration(samples, self.model, self.parameters, input_statistics=missing)

    def test_model_and_all_coefficients_fit_inside_the_graph_at_report_size(self) -> None:
        for name, model, parameters, expected in (
            ("4pl", FourPLModel(), FourPLParameters(4.5, -1.2, 30, 0.08), "4PL\nA=4.5; B=-1.2\nC=30; D=0.08"),
            ("5pl", FivePLModel(), FivePLParameters(2.5, -1.3, 2, 0.1, 0.4), "5PL\nA=2.5; B=-1.3; C=2\nD=0.1; E=0.4"),
            ("quadratic", QuadraticModel(), QuadraticParameters(0.1, 0.02, 0.0001), "QUADRATIC\nA=0.1; B=0.02\nC=0.0001"),
            ("cubic", CubicModel(), CubicParameters(0.1, 0.02, 0.0001, 0.000001), "CUBIC\nA=0.1; B=0.02\nC=0.0001; D=1e-06"),
            ("power_law", PowerLawModel(), PowerLawParameters(0.3, 1.7, 0.12), "POWER_LAW\nA=0.3; B=1.7\nC=0.12"),
        ):
            with self.subTest(model=name):
                figure = plot_calibration(
                    samples_with_statistics(self.samples.to_dict("records")), model, parameters,
                    input_statistics=self.input_statistics, model_name=name, figsize=(3.8, 2.44),
                )
                figure.canvas.draw()
                axes = figure.axes[0]
                legend = axes.get_legend()
                self.assertEqual(legend.get_title().get_text(), expected)
                np.testing.assert_allclose(figure.get_size_inches(), (3.8, 2.44))
                bounds = legend.get_window_extent(figure.canvas.get_renderer())
                self.assertGreaterEqual(bounds.x0, axes.bbox.x0)
                self.assertLessEqual(bounds.x1, axes.bbox.x1)
                self.assertGreaterEqual(bounds.y0, axes.bbox.y0)
                self.assertLessEqual(bounds.y1, axes.bbox.y1)

    def setUp(self) -> None:
        self.model = FourPLModel()
        self.parameters = reference_parameters()
        self.samples = pd.DataFrame(
            {
                "sample_type": ["Standard"] * 4 + ["Unknown"],
                "sample_concentration": [0, 1, 5, 25, 5],
                "extinction_mean": [100] * 5,
                "curve_extinction_mean": [0.2, 0.7, 1.2, 1.8, 999],
                "curve_extinction_sd": [0.01, np.nan, 0.05, np.inf, 10],
            }
        )
        self.input_statistics = (Statistics(0.2, 0.01), Statistics(0.7), Statistics(1.2, 0.05), Statistics(1.8))

    def test_plot_uses_sample_means_and_finite_sds_on_the_fit_scale(self) -> None:
        original = self.samples.copy(deep=True)
        figure = plot_calibration(
            samples_with_statistics(self.samples.to_dict("records")), self.model, self.parameters,
            input_statistics=self.input_statistics, concentration_unit="U/ml",
        )
        axes = figure.axes[0]
        self.assertEqual(axes.get_xscale(), "log")
        self.assertEqual(axes.get_title(), "")
        self.assertEqual(axes.get_yscale(), "linear")
        self.assertEqual(axes.get_xlabel(), "Standaard concentratie [U/ml]")
        measurements = axes.collections[0].get_offsets()
        np.testing.assert_allclose(measurements[:, 0], [1, 5, 25])
        np.testing.assert_allclose(measurements[:, 1], [0.7, 1.2, 1.8])
        fit = axes.lines[0]
        self.assertEqual(fit.get_xdata()[0], 1)
        self.assertEqual(fit.get_xdata()[-1], 25)
        self.assertTrue((fit.get_xdata() > 0).all())
        np.testing.assert_allclose(np.diff(np.log(fit.get_xdata())), np.log(25) / (len(fit.get_xdata()) - 1))
        np.testing.assert_allclose(fit.get_ydata(), self.model.calculate_extinctions(fit.get_xdata(), self.parameters))
        self.assertEqual(len(axes.containers), 1)
        bars = axes.containers[0].lines[2][0].get_segments()
        np.testing.assert_allclose(
            bars,
            [
                [[5, 1.15], [5, 1.25]],
            ],
        )
        pd.testing.assert_frame_equal(self.samples, original)

    def test_standard_means_remain_visible_without_an_sd_column(self) -> None:
        samples = self.samples.drop(columns=["curve_extinction_mean", "curve_extinction_sd"])
        samples["extinction_mean"] = [0.2, 0.7, 1.2, 1.8, 999]
        figure = plot_calibration(
            samples_with_statistics(samples.to_dict("records")), self.model, self.parameters,
            input_statistics=tuple(Statistics(value) for value in (0.2, 0.7, 1.2, 1.8)),
        )
        axes = figure.axes[0]
        self.assertEqual(axes.get_xlabel(), "Standaard concentratie")
        np.testing.assert_allclose(axes.collections[0].get_offsets()[:, 1], [0.7, 1.2, 1.8])
        self.assertEqual(len(axes.containers), 0)

    def test_plot_requires_a_finite_mean_at_a_positive_concentration(self) -> None:
        for concentrations, means in (([0, -1, np.nan], [0.2, 0.7, 1.2]), ([0, 1, 5], [0.2, np.nan, np.inf])):
            with self.subTest(concentrations=concentrations):
                samples = pd.DataFrame({
                    "sample_type": ["Standard"] * 3,
                    "sample_concentration": concentrations,
                    "extinction_mean": means,
                })
                with self.assertRaisesRegex(ValueError, 'positieve concentraties'):
                    plot_calibration(
                        samples_with_statistics(samples.to_dict("records")), self.model, self.parameters,
                        input_statistics=tuple(Statistics(value) for value in means),
                    )

    def test_plot_shows_unavailable_r_squared(self) -> None:
        figure = plot_calibration(
            samples_with_statistics(self.samples.to_dict("records")), self.model, self.parameters,
            input_statistics=self.input_statistics, r_squared=None,
        )
        self.assertEqual(figure.axes[0].lines[0].get_label(), "Fit (R² = n/a)")

    def test_plot_renders_the_supplied_r_squared_without_recalculating(self) -> None:
        figure = plot_calibration(
            samples_with_statistics(self.samples.to_dict("records")), self.model, self.parameters,
            input_statistics=self.input_statistics, r_squared=0.12345,
        )
        self.assertEqual(figure.axes[0].lines[0].get_label(), "Fit (R² = 0.1235)")


if __name__ == "__main__":
    unittest.main()
