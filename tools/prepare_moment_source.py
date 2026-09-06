#!/usr/bin/env python3
"""Prepare continuous-source XSW files for moment-source extraction."""

import argparse
import sys
from pathlib import Path

from moment_source_utils import (
    DEFAULT_EPSILON_FRACTION,
    generated_xsw_paths,
    parse_float_parameter,
    result_spectrum_path,
    runtime_directory,
    time_interval,
    validate_original_xsw,
)
from xsw_utils import replace_parameter_in_text


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Create Zmax+epsilon and Zmax-epsilon continuous-source XSW files "
            "from a single-source configuration."
        )
    )
    parser.add_argument("xsw", help="source XSW file with MomentSource=Single")
    parser.add_argument(
        "--epsilon-fraction",
        type=float,
        default=DEFAULT_EPSILON_FRACTION,
        help="epsilon/Zmax (default: %(default)s)",
    )
    parser.add_argument(
        "--force", action="store_true", help="overwrite generated XSW files"
    )
    return parser.parse_args()


def replace_one(text, parameter_name, value):
    updated, count = replace_parameter_in_text(text, parameter_name, value)
    if count != 1:
        raise ValueError(
            "parameter {0} must occur exactly once, found {1}".format(
                parameter_name, count
            )
        )
    return updated


def prepare(args):
    original_path, zmax = validate_original_xsw(args.xsw)
    if not 0 < args.epsilon_fraction < 1:
        raise ValueError("--epsilon-fraction must be between 0 and 1")

    epsilon = args.epsilon_fraction * zmax
    zplus = zmax + epsilon
    zminus = zmax - epsilon
    if zminus <= 0:
        raise ValueError("Zmax-epsilon must be positive")

    zmin = parse_float_parameter(original_path, "Zmin", 0.0)
    if zmin > zminus:
        raise ValueError(
            "Zmin={0:.17g} is above Zminus={1:.17g}; the lower run would not "
            "contain the complete finite-difference shell".format(zmin, zminus)
        )

    plus_path, minus_path = generated_xsw_paths(original_path)
    existing = [path for path in (plus_path, minus_path) if path.exists()]
    if existing and not args.force:
        raise ValueError(
            "generated file already exists (use --force to overwrite): {0}".format(
                ", ".join(str(path) for path in existing)
            )
        )

    interval = time_interval(original_path, zminus, zplus)
    original_microstep = parse_float_parameter(original_path, "microStep")
    if original_microstep <= 0:
        raise ValueError("microStep must be positive")
    maximum_microstep = 0.5 * interval["light_mpc"]
    generated_microstep = min(original_microstep, maximum_microstep)

    with open(original_path, "r") as input_file:
        original_text = input_file.read()
    base_text = replace_one(original_text, "MomentSource", "0")
    if generated_microstep < original_microstep:
        base_text = replace_one(
            base_text, "microStep", format(generated_microstep, ".17g")
        )
        print(
            "microStep reduced from {0:.17g} Mpc to {1:.17g} Mpc "
            "(0.5 * deltaT light distance).".format(
                original_microstep, generated_microstep
            )
        )

    plus_text = replace_one(base_text, "Zmax", format(zplus, ".17g"))
    minus_text = replace_one(base_text, "Zmax", format(zminus, ".17g"))

    for path, text in ((plus_path, plus_text), (minus_path, minus_text)):
        with open(path, "w") as output_file:
            output_file.write(text)

    plus_result = result_spectrum_path(original_path, plus_path)
    minus_result = result_spectrum_path(original_path, minus_path)
    print("Created {0}".format(plus_path))
    print("Created {0}".format(minus_path))
    print("Zminus={0:.17g}, Zplus={1:.17g}".format(zminus, zplus))
    print(
        "deltaT={0:.17g} years ({1:.17g} light-travel Mpc)".format(
            interval["years"], interval["light_mpc"]
        )
    )
    runtime_path = runtime_directory()
    print("Run from {0}:".format(runtime_path))
    print("  ./propagation {0}".format(plus_path))
    print("  ./propagation {0}".format(minus_path))
    print("Expected spectra:")
    print("  {0}".format(plus_result))
    print("  {0}".format(minus_result))
    print("Then run:")
    extractor = Path(__file__).resolve().with_name("extract_moment_source.py")
    print("  python3 {0} {1}".format(extractor, original_path))


def main():
    try:
        prepare(parse_arguments())
    except (OSError, ValueError) as error:
        print("error: {0}".format(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
