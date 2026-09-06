import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import addPar


class AddParTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary_directory.name)

    def tearDown(self):
        self.temporary_directory.cleanup()

    def write(self, name, text):
        path = self.directory / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_inserts_sample_tag_before_the_following_parameter(self):
        sample = self.write(
            "sample.xsw",
            '<config>\n  <par title="Wanted" value="7"/>\n'
            '  <par title="After" value="8"/>\n</config>\n',
        )
        target = self.write(
            "target.xsw",
            '<config>\n  <par title="Before" value="1"/>\n'
            '  <par title="After" value="8"/>\n</config>\n',
        )
        result = addPar.main(["addPar", "Wanted", str(sample), str(target)])
        self.assertEqual(result, 0)
        updated = target.read_text(encoding="utf-8")
        self.assertLess(updated.index('title="Wanted"'), updated.index('title="After"'))
        self.assertEqual(updated.count('title="Wanted"'), 1)

    def test_missing_sample_parameter_leaves_target_unchanged(self):
        sample = self.write(
            "sample.xsw", '<config>\n  <par title="Other" value="7"/>\n</config>\n'
        )
        target = self.write("target.xsw", "<config/>\n")
        original = target.read_text(encoding="utf-8")
        with contextlib.redirect_stdout(io.StringIO()):
            result = addPar.main(["addPar", "Wanted", str(sample), str(target)])
        self.assertEqual(result, 1)
        self.assertEqual(target.read_text(encoding="utf-8"), original)

    def test_missing_neighbour_leaves_target_unchanged(self):
        sample = self.write(
            "sample.xsw",
            '<config>\n  <par title="Wanted" value="7"/>\n'
            '  <par title="After" value="8"/>\n</config>\n',
        )
        target = self.write(
            "target.xsw", '<config>\n  <par title="Before" value="1"/>\n</config>\n'
        )
        original = target.read_text(encoding="utf-8")
        with contextlib.redirect_stdout(io.StringIO()):
            result = addPar.main(["addPar", "Wanted", str(sample), str(target)])
        self.assertEqual(result, 0)
        self.assertEqual(target.read_text(encoding="utf-8"), original)
        self.assertFalse(Path(str(target) + ".tmp").exists())


if __name__ == "__main__":
    unittest.main()
