import io
import contextlib
import tempfile
import unittest
from pathlib import Path

import add


class AddTests(unittest.TestCase):
    def write_table(self, directory, name, text):
        path = Path(directory) / name
        path.write_text(text, encoding="utf-8")
        return path

    def parse_output(self, output):
        rows = [
            [float(value) for value in line.split()]
            for line in output.getvalue().splitlines()
            if line.strip()
        ]
        return [list(column) for column in zip(*rows)]

    def test_read_table_ignores_comments_and_blank_lines(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_table(
                directory,
                "spectrum",
                "# header\n\n  1 2  \n   \n10 20\n",
            )
            self.assertEqual(add.readTable(str(path)), [[1.0, 10.0], [2.0, 20.0]])

    def test_read_table_rejects_variable_column_count(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_table(directory, "spectrum", "1 2\n10 20 30\n")
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                add.readTable(str(path))

    def test_add_rejects_components_with_different_column_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            first = self.write_table(directory, "first", "1 2\n10 20\n")
            second = self.write_table(directory, "second", "1 2 3\n10 20 30\n")
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                add.add([1, 1], [first, second], output=io.StringIO())

    def test_negative_values_can_be_preserved_or_clamped(self):
        with tempfile.TemporaryDirectory() as directory:
            plus = self.write_table(directory, "plus", "1 1\n10 1\n")
            minus = self.write_table(directory, "minus", "1 2\n10 2\n")

            preserved = io.StringIO()
            add.add([1, -1], [plus, minus], output=preserved, zero_negative=False)
            self.assertEqual(self.parse_output(preserved)[1], [-1, -1])

            clamped = io.StringIO()
            add.add([1, -1], [plus, minus], output=clamped, zero_negative=True)
            self.assertEqual(self.parse_output(clamped)[1], [0, 0])

    def test_output_grid_spacing_follows_last_component(self):
        with tempfile.TemporaryDirectory() as directory:
            coarse = self.write_table(directory, "coarse", "1 1\n10 1\n100 1\n")
            fine = self.write_table(directory, "fine", "1 1\n2 1\n4 1\n8 1\n")

            fine_last = io.StringIO()
            add.add([1, 0], [coarse, fine], output=fine_last)
            coarse_last = io.StringIO()
            add.add([0, 1], [fine, coarse], output=coarse_last)

            fine_grid = self.parse_output(fine_last)[0]
            coarse_grid = self.parse_output(coarse_last)[0]
            self.assertEqual(fine_grid, [1, 2, 4, 8, 16, 32, 64])
            self.assertEqual(coarse_grid, [1, 10, 100])


if __name__ == "__main__":
    unittest.main()
