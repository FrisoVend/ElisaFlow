import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from elisaflow.configuration.application import ApplicationSettings, default_config_path


class ApplicationSettingsTests(unittest.TestCase):
    def test_password_file_defaults_to_the_configuration_directory(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary).resolve()
            config = directory / "app.toml"
            config.write_text("[paths]\n", encoding="utf-8")
            self.assertEqual(ApplicationSettings.from_file(config).password_file, directory / ".secrets" / "assay-editor.hash")

    def test_password_file_accepts_relative_and_absolute_paths(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary).resolve()
            config = directory / "app.toml"
            for value in ("private/editor.hash", (directory / "shared" / "editor.hash").as_posix()):
                with self.subTest(value=value):
                    config.write_text(f"[access_control]\npassword_file = '{value}'\n", encoding="utf-8")
                    self.assertEqual(ApplicationSettings.from_file(config).password_file, (directory / value).resolve())

    def test_invalid_access_control_configuration_is_rejected(self) -> None:
        with TemporaryDirectory() as temporary:
            config = Path(temporary) / "app.toml"
            for contents in ('access_control = "invalid"', '[access_control]\npassword_file = ""', '[access_control]\npassword_file = 123'):
                with self.subTest(contents=contents):
                    config.write_text(contents, encoding="utf-8")
                    with self.assertRaises(ValueError):
                        ApplicationSettings.from_file(config)

    def test_frozen_configuration_is_found_beside_executable(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "app.toml").write_text('[paths]\nresults = "exports"\n', encoding="utf-8")
            with (
                patch("elisaflow.configuration.application.sys.frozen", True, create=True),
                patch("elisaflow.configuration.application.sys.executable", str(directory / "ElisaFlow.exe")),
                patch("elisaflow.configuration.application.Path.cwd", return_value=directory / "unrelated"),
            ):
                config_path = default_config_path()
                settings = ApplicationSettings.from_file(config_path)
            self.assertEqual(config_path, directory / "app.toml")
            self.assertEqual(settings.results_dir, directory / "exports")

    def test_frozen_configuration_ignores_local_overrides(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            local = directory / "app.local.toml"
            local.write_text('not valid TOML = [', encoding="utf-8")
            (directory / "app.toml").write_text("[paths]\n", encoding="utf-8")
            with (
                patch("elisaflow.configuration.application.sys.frozen", True, create=True),
                patch("elisaflow.configuration.application.sys.executable", str(directory / "ElisaFlow.exe")),
            ):
                config = default_config_path()
                self.assertEqual(config, directory / "app.toml")
                self.assertEqual(ApplicationSettings.from_file(config).results_dir, directory / "resultaten")

    def test_source_configuration_uses_working_directory(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            with (
                patch("elisaflow.configuration.application.sys.frozen", False, create=True),
                patch("elisaflow.configuration.application.Path.cwd", return_value=directory),
            ):
                self.assertEqual(default_config_path(), directory / "app.toml")
                local = directory / "app.local.toml"
                local.write_text('not valid TOML = [', encoding="utf-8")
                (directory / "app.toml").write_text("[paths]\n", encoding="utf-8")
                self.assertEqual(default_config_path(), directory / "app.toml")
                self.assertEqual(ApplicationSettings.from_file(default_config_path()).results_dir, directory / "resultaten")

    def test_missing_app_toml_does_not_fall_back_to_a_local_file(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "app.local.toml").write_text("[paths]\n", encoding="utf-8")
            for frozen in (False, True):
                with (
                    self.subTest(frozen=frozen),
                    patch("elisaflow.configuration.application.sys.frozen", frozen, create=True),
                    patch("elisaflow.configuration.application.sys.executable", str(directory / "ElisaFlow.exe")),
                    patch("elisaflow.configuration.application.Path.cwd", return_value=directory),
                ):
                    self.assertEqual(default_config_path(), directory / "app.toml")
                    with self.assertRaises(FileNotFoundError):
                        ApplicationSettings.from_file(default_config_path())

    def test_shipped_configuration_uses_the_standard_results_folder(self) -> None:
        directory = Path(__file__).resolve().parents[1]
        self.assertEqual(ApplicationSettings.from_file(directory / "app.toml").results_dir, directory / "resultaten")


if __name__ == "__main__":
    unittest.main()
