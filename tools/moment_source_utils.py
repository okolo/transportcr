"""Shared helpers for the two-stage moment-source extraction tools."""

from __future__ import print_function

import math
from pathlib import Path

from xsw_utils import XswParameterError, get_parameter


DEFAULT_EPSILON_FRACTION = 0.01
DEFAULT_H_KM_S_MPC = 68.17
DEFAULT_LAMBDA_DENSITY = 0.6973
MPC_IN_KM = 3.0856775807e19
LIGHT_SPEED_KM_S = 2.99792458e5
JULIAN_YEAR_SECONDS = 31557600.0


def parse_float_parameter(path, name, default=None):
    try:
        raw_value = get_parameter(path, name)
    except XswParameterError:
        if default is None:
            raise
        return float(default)
    try:
        return float(raw_value)
    except ValueError:
        raise XswParameterError(
            "parameter {0} in {1} is not a number: {2}".format(name, path, raw_value)
        )


def parse_bool_parameter(path, name):
    value = get_parameter(path, name).strip().lower()
    if value == "true":
        return True
    if value == "false":
        return False
    raise XswParameterError(
        "parameter {0} in {1} must be true or false, got {2}".format(
            name, path, value
        )
    )


def validate_original_xsw(path):
    path = Path(path).resolve()
    if not path.is_file():
        raise ValueError("XSW file not found: {0}".format(path))
    if get_parameter(path, "MomentSource") != "1":
        raise ValueError("MomentSource must be Single (value=1) in {0}".format(path))
    if not parse_bool_parameter(path, "CELredshift"):
        raise ValueError("CELredshift must be true in {0}".format(path))
    zmax = parse_float_parameter(path, "Zmax")
    if zmax <= 0:
        raise ValueError("Zmax must be positive in {0}".format(path))
    return path, zmax


def generated_xsw_paths(original_path):
    original_path = Path(original_path).resolve()
    return (
        original_path.with_name(original_path.stem + "_zplus" + original_path.suffix),
        original_path.with_name(original_path.stem + "_zminus" + original_path.suffix),
    )


def runtime_directory():
    """Return the repository bin directory from which propagation must run."""
    return Path(__file__).resolve().parent.parent / "bin"


def cosmology(path):
    hubble = parse_float_parameter(
        path, "H_in_km_s_Mpc", DEFAULT_H_KM_S_MPC
    )
    lambda_density = parse_float_parameter(
        path, "Lv", DEFAULT_LAMBDA_DENSITY
    )
    if hubble <= 0:
        raise ValueError("H_in_km_s_Mpc must be positive")
    if not 0 <= lambda_density <= 1:
        raise ValueError("Lv must be between 0 and 1")
    return hubble, lambda_density


def delta_time_seconds(zminus, zplus, hubble_km_s_mpc, lambda_density):
    """Return the flat-LambdaCDM time interval between two redshifts."""
    if not 0 <= zminus < zplus:
        raise ValueError("expected 0 <= Zminus < Zplus")

    matter_density = 1.0 - lambda_density

    def integrand(redshift):
        one_plus_z = 1.0 + redshift
        expansion = math.sqrt(
            matter_density * one_plus_z * one_plus_z * one_plus_z
            + lambda_density
        )
        return 1.0 / (one_plus_z * expansion)

    # Composite Simpson integration is stable for the deliberately narrow interval.
    intervals = 256
    step = (zplus - zminus) / intervals
    total = integrand(zminus) + integrand(zplus)
    for index in range(1, intervals):
        weight = 4 if index % 2 else 2
        total += weight * integrand(zminus + index * step)
    dimensionless_interval = total * step / 3.0

    hubble_per_second = hubble_km_s_mpc / MPC_IN_KM
    return dimensionless_interval / hubble_per_second


def time_interval(path, zminus, zplus):
    hubble, lambda_density = cosmology(path)
    seconds = delta_time_seconds(zminus, zplus, hubble, lambda_density)
    return {
        "seconds": seconds,
        "years": seconds / JULIAN_YEAR_SECONDS,
        "light_mpc": seconds * LIGHT_SPEED_KM_S / MPC_IN_KM,
        "hubble": hubble,
        "lambda_density": lambda_density,
    }


def result_spectrum_path(original_path, generated_path):
    del original_path
    return runtime_directory() / "results" / generated_path.stem / "uniform" / "all"


def extracted_spectrum_path(original_path):
    original_path = Path(original_path).resolve()
    return (
        runtime_directory()
        / "results"
        / (original_path.stem + "_moment_source_fd")
        / "uniform"
        / "all"
    )
