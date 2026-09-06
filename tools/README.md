# TransportCR command-line tools

This directory contains small utilities for inspecting and editing TransportCR
XSW configuration files, combining spectra, and extracting a point-source
spectrum from two continuous-source calculations.

The tools require Python 3.8 or newer. They have no third-party runtime
dependencies.

## Installation

The commands use paths relative to this source checkout. Install them in
editable mode from the repository root:

```sh
python3 -m pip install --editable ./tools
```

A non-editable installation is not supported because `helpPar` and the
moment-source tools need to locate the repository `bin` directory.

The build uses `setuptools`, which is declared in `pyproject.toml`. A modern
version of `pip` installs this build dependency automatically.

After installation, the following commands are available:

- `getPar`
- `replacePar`
- `add`
- `addPar`
- `helpPar`
- `prepare_moment_source`
- `extract_moment_source`

The corresponding `.py` files can also be executed directly with Python 3.

## XSW parameter tools

### `getPar`

Print the value of a parameter from an XSW file:

```sh
getPar Zmax bin/switches.xsw
```

Usage:

```text
getPar <parameter-name> <file.xsw>
```

### `replacePar`

Replace one parameter value in one or more XSW files:

```sh
replacePar Zmax 1.0 run1.xsw run2.xsw
```

Usage:

```text
replacePar <parameter-name> <new-value> <file1.xsw> [file2.xsw ...]
```

If the parameter occurs more than once, the original file is left unchanged
and the proposed result is written to a sibling `.tmp` file.

Values containing `&`, `<`, a single quote, or a double quote are rejected
because these characters cannot be inserted directly into an XML attribute.

### `addPar`

Copy an XSW parameter tag from a sample configuration into one or more target
files:

```sh
addPar NewParameter sample.xsw run1.xsw run2.xsw
```

Usage:

```text
addPar <parameter-name> <sample.xsw> <file1.xsw> [file2.xsw ...]
```

This legacy utility does not check whether the target already contains the
parameter, so repeated use can create duplicate parameter tags.

### `helpPar`

Print the description of an XSW parameter. For switch parameters, also print
the available numbered values:

```sh
helpPar MomentSource bin/switches.xsw
```

Usage:

```text
helpPar <parameter-name> [file.xsw]
```

When the file is omitted, the default is `../bin/switches.xsw` relative to the
location of `helpPar.py`.

## Spectrum combination

### `add`

Calculate a linear combination of tabulated spectra:

```sh
add 1 spectrum1 -1 spectrum2 > difference
```

Usage:

```text
add <multiplier1> <spectrum1> [<multiplier2> <spectrum2> ...]
```

The first column is treated as the energy scale. Other columns are combined
with the supplied multipliers. The standalone `add` command preserves its
historical behavior of replacing negative results with zero.

The output-grid spacing is intentionally taken from the last input spectrum.
Consequently, changing the order of spectra with different grid resolutions
can change the output grid. This behavior is retained for compatibility and
may be revised in a future version.

## Two-stage moment-source extraction

This workflow approximates a point source at `Zmax` using a centered finite
difference of two continuous-source calculations:

```text
F_point = (F_plus - F_minus) / deltaT_years
```

where `Zplus = Zmax + epsilon`, `Zminus = Zmax - epsilon`, and the default
`epsilon` is `0.01 * Zmax`.

The input XSW file must have `MomentSource=1` (`Single`) and
`CELredshift=true`.

### Stage 1: prepare configurations

Run:

```sh
prepare_moment_source path/to/source.xsw
```

The command creates these files next to the input file:

```text
source_zplus.xsw
source_zminus.xsw
```

Both generated files use `MomentSource=0` (`Off`). If the original
`microStep` is greater than half the light-travel distance corresponding to
the `Zminus`--`Zplus` time interval, it is reduced to that value in both
files. A message is printed only when this adjustment is made.

Options:

```text
--epsilon-fraction VALUE  Set epsilon/Zmax; default: 0.01
--force                   Overwrite existing generated XSW files
```

The command prints the exact propagation commands and expected result paths.
Run `propagation` from the repository `bin` directory because TransportCR
resolves its runtime tables relative to the current directory.

### Stage 2: extract the spectrum

After both propagation runs finish, pass the original XSW file to:

```sh
extract_moment_source path/to/source.xsw
```

The command finds the two spectra under `bin/results`, subtracts them, divides
the result by the cosmological time interval in Julian years, and writes:

```text
bin/results/source_moment_source_fd/uniform/all
```

Negative values are preserved by default so that numerical subtraction noise
remains visible. To replace them with zero, use:

```sh
extract_moment_source --zero-negative path/to/source.xsw
```

Use `--force` to overwrite an existing extracted spectrum.

The time interval is calculated for a flat Lambda-CDM cosmology using
`H_in_km_s_Mpc` and `Lv` from the original XSW file.
