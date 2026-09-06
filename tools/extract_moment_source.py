#!/usr/bin/env python3
"""Extract a normalized moment-source spectrum from two continuous runs."""

import argparse
import io
import sys

from add import add, readTable
from moment_source_utils import (
    extracted_spectrum_path,
    generated_xsw_paths,
    parse_bool_parameter,
    parse_float_parameter,
    result_spectrum_path,
    time_interval,
    validate_original_xsw,
)
from xsw_utils import get_parameter


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Subtract the prepared continuous-source spectra and divide by "
            "their cosmological time interval."
        )
    )
    parser.add_argument("xsw", help="original XSW file with MomentSource=Single")
    parser.add_argument(
        "--zero-negative",
        action="store_true",
        help="replace negative values in the extracted spectrum with zero",
    )
    parser.add_argument(
        "--force", action="store_true", help="overwrite an existing output spectrum"
    )
    return parser.parse_args()


def validate_generated(path, expected_zmax):
    if not path.is_file():
        raise ValueError("generated XSW file not found: {0}".format(path))
    if get_parameter(path, "MomentSource") != "0":
        raise ValueError("MomentSource must be Off (value=0) in {0}".format(path))
    if not parse_bool_parameter(path, "CELredshift"):
        raise ValueError("CELredshift must be true in {0}".format(path))
    actual_zmax = parse_float_parameter(path, "Zmax")
    tolerance = max(1.0, abs(expected_zmax)) * 1e-14
    if abs(actual_zmax - expected_zmax) > tolerance:
        raise ValueError(
            "unexpected Zmax in {0}: expected {1:.17g}, got {2:.17g}".format(
                path, expected_zmax, actual_zmax
            )
        )


def leading_comments(path):
    comments = []
    with open(path, "r") as spectrum:
        for line in spectrum:
            if line.startswith("#"):
                comments.append(line)
            elif line.strip():
                break
    return comments


def validate_spectrum_grids(plus_spectrum, minus_spectrum):
    plus_data = readTable(str(plus_spectrum))
    minus_data = readTable(str(minus_spectrum))
    if len(plus_data) != len(minus_data):
        raise ValueError("input spectra have different numbers of columns")
    if not plus_data or not minus_data:
        raise ValueError("input spectra contain no numeric columns")
    if plus_data[0] != minus_data[0]:
        raise ValueError("input spectra have different energy grids")


def extract(args):
    original_path, zmax = validate_original_xsw(args.xsw)
    plus_xsw, minus_xsw = generated_xsw_paths(original_path)
    zplus = parse_float_parameter(plus_xsw, "Zmax") if plus_xsw.is_file() else None
    zminus = parse_float_parameter(minus_xsw, "Zmax") if minus_xsw.is_file() else None
    if zplus is None or zminus is None:
        missing = plus_xsw if zplus is None else minus_xsw
        raise ValueError("generated XSW file not found: {0}".format(missing))
    if not zminus < zmax < zplus:
        raise ValueError(
            "generated redshifts must satisfy Zminus < Zmax < Zplus"
        )
    if abs((zplus - zmax) - (zmax - zminus)) > max(1.0, zmax) * 1e-14:
        raise ValueError("generated redshifts are not symmetric around Zmax")

    validate_generated(plus_xsw, zplus)
    validate_generated(minus_xsw, zminus)

    plus_spectrum = result_spectrum_path(original_path, plus_xsw)
    minus_spectrum = result_spectrum_path(original_path, minus_xsw)
    for spectrum in (plus_spectrum, minus_spectrum):
        if not spectrum.is_file():
            raise ValueError("spectrum not found: {0}".format(spectrum))
    validate_spectrum_grids(plus_spectrum, minus_spectrum)

    interval = time_interval(original_path, zminus, zplus)
    if interval["years"] <= 0:
        raise ValueError("calculated deltaT must be positive")

    output_path = extracted_spectrum_path(original_path)
    if output_path.exists() and not args.force:
        raise ValueError(
            "output spectrum already exists (use --force to overwrite): {0}".format(
                output_path
            )
        )

    combined = io.StringIO()
    inverse_delta_t = 1.0 / interval["years"]
    add(
        [inverse_delta_t, -inverse_delta_t],
        [str(plus_spectrum), str(minus_spectrum)],
        output=combined,
        zero_negative=args.zero_negative,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as output_file:
        for comment in leading_comments(plus_spectrum):
            output_file.write(comment)
        output_file.write(
            "# Extracted as (F_plus-F_minus)/deltaT_years; "
            "Zminus={0:.17g}; Zplus={1:.17g}; deltaT_years={2:.17g}\n".format(
                zminus, zplus, interval["years"]
            )
        )
        output_file.write(combined.getvalue())

    print("Zminus={0:.17g}, Zplus={1:.17g}".format(zminus, zplus))
    print("deltaT={0:.17g} years".format(interval["years"]))
    print("Created {0}".format(output_path))


def main():
    try:
        extract(parse_arguments())
    except (OSError, ValueError, AssertionError) as error:
        print("error: {0}".format(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
