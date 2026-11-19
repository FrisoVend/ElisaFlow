import io
import unittest
from unittest.mock import patch

from elisaflow.runs.analyst import validate_analyst_name
from elisaflow.ui.desktop import prompt_analyst_name


class UserIdentityTests(unittest.TestCase):
    def test_entered_names_are_trimmed_and_preserved(self):
        for name in ("Friso van Dijk", "Zoë de Vries", "Anne-Marie O'Neil", r"LAB\operator"):
            with self.subTest(name=name):
                self.assertEqual(validate_analyst_name(f"  {name}  "), name)

    def test_names_must_be_nonempty_and_without_control_characters(self):
        for name in (None, 42, "", "   ", "Name\nOther", "Name\rOther", "Name\tOther", "Name\x00", "Name\x7f"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                validate_analyst_name(name)

    def test_console_prompt_retries_invalid_names(self):
        with (
            patch("elisaflow.ui.desktop.sys.stdin", io.StringIO()),
            patch("builtins.input", side_effect=["   ", "Name\nOther", "  Friso van Dijk  "]) as enter,
            patch("builtins.print") as message,
        ):
            self.assertEqual(prompt_analyst_name(), "Friso van Dijk")
        self.assertEqual(enter.call_count, 3)
        self.assertEqual(message.call_count, 2)

    def test_console_prompt_can_be_cancelled(self):
        for error in (EOFError(), KeyboardInterrupt()):
            with self.subTest(error=type(error)), patch("elisaflow.ui.desktop.sys.stdin", io.StringIO()), patch("builtins.input", side_effect=error):
                self.assertIsNone(prompt_analyst_name())


if __name__ == "__main__":
    unittest.main()
