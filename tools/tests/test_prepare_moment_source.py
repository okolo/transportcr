import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import prepare_moment_source
from moment_source_utils import generated_xsw_paths, time_interval
from tests.test_support import write_xsw
from xsw_utils import get_parameter


class PrepareMomentSourceTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary_directory.name)

    def tearDown(self):
        self.temporary_directory.cleanup()

    def arguments(self, path, epsilon_fraction=0.01, force=False):
        return SimpleNamespace(
            xsw=str(path), epsilon_fraction=epsilon_fraction, force=force
        )

    def test_creates_symmetric_continuous_source_configurations(self):
        original = self.directory / "source.xsw"
        write_xsw(original, Zmax="0.1", microStep="0.001")
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            prepare_moment_source.prepare(self.arguments(original))

        plus, minus = generated_xsw_paths(original)
        self.assertAlmostEqual(float(get_parameter(plus, "Zmax")), 0.101)
        self.assertAlmostEqual(float(get_parameter(minus, "Zmax")), 0.099)
        self.assertEqual(get_parameter(plus, "MomentSource"), "0")
        self.assertEqual(get_parameter(minus, "MomentSource"), "0")
        self.assertEqual(get_parameter(plus, "microStep"), "0.001")
        self.assertNotIn("microStep reduced", output.getvalue())

    def test_reduces_microstep_to_half_the_interval_light_distance(self):
        original = self.directory / "source.xsw"
        write_xsw(original, Zmax="0.1", microStep="100")
        expected = 0.5 * time_interval(original, 0.099, 0.101)["light_mpc"]
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            prepare_moment_source.prepare(self.arguments(original))

        plus, minus = generated_xsw_paths(original)
        self.assertAlmostEqual(float(get_parameter(plus, "microStep")), expected)
        self.assertAlmostEqual(float(get_parameter(minus, "microStep")), expected)
        self.assertIn("microStep reduced", output.getvalue())

    def test_rejects_wrong_source_mode_and_disabled_redshift_losses(self):
        for overrides, message in (
            ({"MomentSource": "0"}, "MomentSource must be Single"),
            ({"CELredshift": "false"}, "CELredshift must be true"),
        ):
            with self.subTest(overrides=overrides):
                original = self.directory / ("source-{0}.xsw".format(len(message)))
                write_xsw(original, **overrides)
                with self.assertRaisesRegex(ValueError, message):
                    prepare_moment_source.prepare(self.arguments(original))

    def test_rejects_invalid_epsilon_and_zmin_above_zminus(self):
        original = self.directory / "source.xsw"
        write_xsw(original)
        for epsilon in (0, 1, -0.1, 1.1):
            with self.subTest(epsilon=epsilon):
                with self.assertRaisesRegex(ValueError, "epsilon-fraction"):
                    prepare_moment_source.prepare(
                        self.arguments(original, epsilon_fraction=epsilon)
                    )

        write_xsw(original, Zmin="0.0995")
        with self.assertRaisesRegex(ValueError, "above Zminus"):
            prepare_moment_source.prepare(self.arguments(original))

    def test_existing_outputs_require_force(self):
        original = self.directory / "source.xsw"
        write_xsw(original)
        plus, minus = generated_xsw_paths(original)
        plus.write_text("sentinel\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "use --force"):
            prepare_moment_source.prepare(self.arguments(original))

        with contextlib.redirect_stdout(io.StringIO()):
            prepare_moment_source.prepare(self.arguments(original, force=True))
        self.assertNotEqual(plus.read_text(encoding="utf-8"), "sentinel\n")
        self.assertTrue(minus.is_file())


if __name__ == "__main__":
    unittest.main()
