import importlib
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import helpPar


TOOLS_DIRECTORY = Path(__file__).resolve().parents[1]
PYPROJECT = TOOLS_DIRECTORY / "pyproject.toml"
EXPECTED_ENTRY_POINTS = {
    "add": "add:main",
    "addPar": "addPar:main",
    "extract_moment_source": "extract_moment_source:main",
    "getPar": "getPar:main",
    "helpPar": "helpPar:main",
    "prepare_moment_source": "prepare_moment_source:main",
    "replacePar": "replacePar:main",
}


def declared_entry_points():
    text = PYPROJECT.read_text(encoding="utf-8")
    section = re.search(
        r"^\[project\.scripts\]\s*$([\s\S]*?)(?=^\[|\Z)", text, re.MULTILINE
    )
    if section is None:
        return {}
    return dict(
        re.findall(r'^([A-Za-z0-9_]+)\s*=\s*"([^"]+)"\s*$', section.group(1), re.MULTILINE)
    )


class EntryPointTests(unittest.TestCase):
    def test_help_par_default_points_to_repository_switches(self):
        self.assertEqual(
            helpPar.default_xsw_path(),
            TOOLS_DIRECTORY.parent / "bin" / "switches.xsw",
        )

    def test_all_expected_entry_points_resolve_to_main_functions(self):
        self.assertEqual(declared_entry_points(), EXPECTED_ENTRY_POINTS)
        for target in EXPECTED_ENTRY_POINTS.values():
            module_name, function_name = target.split(":", 1)
            with self.subTest(target=target):
                function = getattr(importlib.import_module(module_name), function_name)
                self.assertTrue(callable(function))

    @unittest.skipUnless(
        os.environ.get("TRANSPORTCR_TEST_EDITABLE_INSTALL") == "1",
        "set TRANSPORTCR_TEST_EDITABLE_INSTALL=1 to test an editable installation",
    )
    def test_editable_install_creates_console_commands(self):
        with tempfile.TemporaryDirectory() as directory:
            environment = Path(directory) / "venv"
            subprocess.run([sys.executable, "-m", "venv", str(environment)], check=True)
            python = environment / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
            subprocess.run(
                [
                    str(python),
                    "-m",
                    "pip",
                    "install",
                    "--editable",
                    str(TOOLS_DIRECTORY),
                ],
                check=True,
            )
            scripts = environment / ("Scripts" if os.name == "nt" else "bin")
            for command in EXPECTED_ENTRY_POINTS:
                with self.subTest(command=command):
                    executable = scripts / (command + (".exe" if os.name == "nt" else ""))
                    self.assertTrue(executable.is_file())
                    arguments = ["--help"] if command in {
                        "extract_moment_source",
                        "helpPar",
                        "prepare_moment_source",
                    } else []
                    result = subprocess.run(
                        [str(executable)] + arguments,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                    )
                    expected_status = 0 if command in {
                        "extract_moment_source",
                        "prepare_moment_source",
                    } else 1
                    self.assertEqual(result.returncode, expected_status)


if __name__ == "__main__":
    unittest.main()
