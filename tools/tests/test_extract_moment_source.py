import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import add
import extract_moment_source
from moment_source_utils import generated_xsw_paths, time_interval
from tests.test_support import write_spectrum, write_xsw


class ExtractMomentSourceTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary_directory.name)
        self.original = self.directory / "source.xsw"
        write_xsw(self.original, Zmax="0.1")
        self.plus_xsw, self.minus_xsw = generated_xsw_paths(self.original)
        write_xsw(self.plus_xsw, MomentSource="0", Zmax="0.101")
        write_xsw(self.minus_xsw, MomentSource="0", Zmax="0.099")
        self.plus_spectrum = self.directory / "results" / "plus" / "all"
        self.minus_spectrum = self.directory / "results" / "minus" / "all"
        self.output_spectrum = self.directory / "results" / "extracted" / "all"

    def tearDown(self):
        self.temporary_directory.cleanup()

    def arguments(self, zero_negative=False, force=False):
        return SimpleNamespace(
            xsw=str(self.original), zero_negative=zero_negative, force=force
        )

    def run_extract(self, args=None):
        if args is None:
            args = self.arguments()
        with mock.patch.object(
            extract_moment_source,
            "result_spectrum_path",
            side_effect=[self.plus_spectrum, self.minus_spectrum],
        ), mock.patch.object(
            extract_moment_source,
            "extracted_spectrum_path",
            return_value=self.output_spectrum,
        ), contextlib.redirect_stdout(io.StringIO()):
            extract_moment_source.extract(args)

    def write_normalized_pair(self):
        delta_t = time_interval(self.original, 0.099, 0.101)["years"]
        minus = [(1, 10 * delta_t, 20 * delta_t), (10, 30 * delta_t, 40 * delta_t)]
        plus = [(1, 12 * delta_t, 17 * delta_t), (10, 35 * delta_t, 39 * delta_t)]
        write_spectrum(self.plus_spectrum, plus, comments=("# plus spectrum",))
        write_spectrum(self.minus_spectrum, minus)
        return [(2, -3), (5, -1)]

    def test_extracts_normalized_difference_and_preserves_header(self):
        expected = self.write_normalized_pair()
        self.run_extract()

        text = self.output_spectrum.read_text(encoding="utf-8")
        self.assertTrue(text.startswith("# plus spectrum\n"))
        self.assertIn("(F_plus-F_minus)/deltaT_years", text)
        data = add.readTable(self.output_spectrum)
        for row, expected_row in enumerate(expected):
            self.assertAlmostEqual(data[1][row], expected_row[0], places=12)
            self.assertAlmostEqual(data[2][row], expected_row[1], places=12)

    def test_zero_negative_option_clamps_negative_values(self):
        expected = self.write_normalized_pair()
        self.run_extract(self.arguments(zero_negative=True))
        data = add.readTable(self.output_spectrum)
        for row, expected_row in enumerate(expected):
            self.assertAlmostEqual(data[1][row], expected_row[0], places=12)
            self.assertEqual(data[2][row], 0)

    def test_rejects_different_energy_grids(self):
        write_spectrum(self.plus_spectrum, [(1, 1), (10, 2)])
        write_spectrum(self.minus_spectrum, [(1, 1), (20, 2)])
        with self.assertRaisesRegex(ValueError, "different energy grids"):
            self.run_extract()

    def test_reports_missing_result_spectrum(self):
        write_spectrum(self.minus_spectrum, [(1, 1), (10, 2)])
        with self.assertRaisesRegex(ValueError, "spectrum not found"):
            self.run_extract()

    def test_existing_output_requires_force(self):
        self.write_normalized_pair()
        write_spectrum(self.output_spectrum, [(1, 1), (10, 2)])
        with self.assertRaisesRegex(ValueError, "use --force"):
            self.run_extract()
        self.run_extract(self.arguments(force=True))
        self.assertIn(
            "(F_plus-F_minus)/deltaT_years",
            self.output_spectrum.read_text(encoding="utf-8"),
        )

    def test_rejects_asymmetric_generated_redshifts(self):
        write_xsw(self.plus_xsw, MomentSource="0", Zmax="0.102")
        with self.assertRaisesRegex(ValueError, "not symmetric"):
            self.run_extract()


if __name__ == "__main__":
    unittest.main()
