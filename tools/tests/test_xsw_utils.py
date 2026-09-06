import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import replacePar
from xsw_utils import (
    XswParameterError,
    get_parameter,
    get_parameter_values,
    replace_parameter,
    replace_parameter_in_text,
    validate_parameter_value,
)


class XswUtilsTests(unittest.TestCase):
    def test_reads_single_and_double_quoted_values(self):
        text = (
            '<parameter title="Zmax" value="0.1"/>\n'
            "<parameter title='Zmax' value='0.2'/>\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.xsw"
            path.write_text(text, encoding="utf-8")
            self.assertEqual(get_parameter_values(path, "Zmax"), ["0.1", "0.2"])

    def test_get_parameter_rejects_missing_and_duplicate_parameters(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.xsw"
            path.write_text(
                '<parameter title="Zmax" value="0.1"/>\n'
                '<parameter title="Zmax" value="0.2"/>\n',
                encoding="utf-8",
            )
            with self.assertRaisesRegex(XswParameterError, "occurs 2 times"):
                get_parameter(path, "Zmax")
            with self.assertRaisesRegex(XswParameterError, "not found"):
                get_parameter(path, "Zmin")

    def test_replaces_only_the_exact_parameter_name(self):
        text = (
            '<parameter title="Z" value="1"/>\n'
            '<parameter title="Zmax" value="2"/>\n'
        )
        updated, count = replace_parameter_in_text(text, "Z", "3")
        self.assertEqual(count, 1)
        self.assertIn('title="Z" value="3"', updated)
        self.assertIn('title="Zmax" value="2"', updated)

    def test_replace_parameter_writes_only_for_one_match(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.xsw"
            original = (
                '<parameter title="Zmax" value="0.1"/>\n'
                '<parameter title="Zmax" value="0.2"/>\n'
            )
            path.write_text(original, encoding="utf-8")
            self.assertEqual(replace_parameter(path, "Zmax", "0.3"), 2)
            self.assertEqual(path.read_text(encoding="utf-8"), original)

    def test_rejects_unsupported_xml_attribute_characters(self):
        for character in ('&', '<', '"', "'"):
            with self.subTest(character=character):
                with self.assertRaises(XswParameterError):
                    validate_parameter_value("left" + character + "right")

    def test_replace_par_rejects_invalid_value_without_modifying_file(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.xsw"
            original = '<parameter title="Zmax" value="0.1"/>\n'
            path.write_text(original, encoding="utf-8")
            errors = io.StringIO()
            with contextlib.redirect_stderr(errors):
                result = replacePar.main(["replacePar", "Zmax", "1&2", str(path)])
            self.assertEqual(result, 1)
            self.assertIn("unsupported XML character", errors.getvalue())
            self.assertEqual(path.read_text(encoding="utf-8"), original)


if __name__ == "__main__":
    unittest.main()
