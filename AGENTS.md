# AGENTS.md

## Scope

These instructions apply to the entire repository. They are intended for Codex and other AI coding agents working on TransportCR.

## Project Overview

TransportCR is a scientific simulation of ultra-high-energy cosmic-ray and electromagnetic-cascade propagation based on one-dimensional transport equations. The codebase is primarily legacy C++, with Fortran sources for some interaction models. Runtime configuration is supplied through XML/XSW files, and many physical models depend on tabulated data stored under `bin/tables/`.

The main executable is `propagation`. It must compile and run on Linux and macOS releases no more than five years old.

## Repository Layout

- `src/`: C++, C, and Fortran sources; CMake configuration; legacy Makefiles.
- `src/pp/`: bundled Fortran code for proton-proton interactions.
- `bin/`: the conventional build and runtime directory.
- `bin/tables/`: versioned scientific input tables managed by Git LFS.
- `bin/switches.xsw`: the main example/default runtime configuration.
- `bin/switches.dtd`: the schema used to validate XML/XSW configuration files.
- `bin/results/`: generated simulation output; do not treat it as source data.

Do not assume that everything under `bin/` is disposable: it contains both generated build artifacts and tracked runtime data.

## Language Requirements

- Write `AGENTS.md`, committed documentation, commit messages, and all source-code comments in English.
- Any text that will be committed to the repository must be in English unless the task explicitly requires language-specific user data.
- Temporary working documents, such as `todo.md`, may be written in Russian, but they must remain untracked and must not be committed.

## Build Environment

The supported toolchain requires:

- CMake and Make.
- A C and C++ compiler.
- A Fortran compiler, normally `gfortran`.
- Xerces-C development libraries.
- GSL development libraries.
- Git LFS for the scientific tables.

The documented Ubuntu packages are:

```sh
sudo apt install cmake make gcc g++ gfortran git git-lfs libxerces-c-dev libgsl-dev
```

Agents may install required dependencies. Prefer the platform's standard package manager and avoid unrelated upgrades. On macOS, use mutually compatible architectures for the compiler, executable, Xerces-C, GSL, and Fortran runtime; do not mix arm64 and x86_64 libraries.

Before building a fresh checkout, ensure the LFS objects are available:

```sh
git lfs install
git lfs pull
```

The repository's established build flow is:

```sh
cd bin
cmake ../src
make
```

`bin/install.sh` performs the same configure-and-build sequence and may also be used:

```sh
cd bin
. ./install.sh
```

Run the executable from `bin/`, because runtime paths such as `tables/` and `switches.dtd` are resolved relative to the working directory.

## Required Verification

For every source or build-system change, the current minimum verification is a successful build on the available supported platform:

```sh
cd bin
cmake ../src
make
```

Report the operating system, compiler, build command, and outcome. If the build cannot be completed because a dependency or another supported platform is unavailable, state that explicitly; do not claim that the change is verified there.

The executable contains a limited `--test BLLac` path, but it is not currently the required verification suite. Do not present it as comprehensive unit-test coverage. Regression tests based on selected `.xsw` files and canonical numerical results are planned but not yet defined.

Do not start a simulation or other computation expected to take more than five minutes without first obtaining explicit user approval. Before requesting approval, identify the command, configuration, approximate duration, and expected output location.

## Change Policy

- Make the smallest change that satisfies the task.
- Preserve the existing architecture, APIs, naming, formatting, and language level unless the user explicitly requests otherwise.
- Do not perform opportunistic refactoring, broad cleanup, modernization, or migration to modern C++.
- Do not make large structural changes without explicit approval.
- Do not update CMake files without prior approval, except for the minimal edit required to add new source files to the build.
- For any other CMake change, explain why it is needed and obtain permission before editing.
- Treat the legacy Makefiles as historical build paths unless the task specifically concerns them. Do not try to synchronize or modernize them automatically.
- Keep platform-specific behavior narrowly guarded and preserve support for both Linux and macOS.
- Do not suppress compiler diagnostics or weaken runtime assertions merely to make a build pass.

When a requested change appears to require refactoring or modernization, stop before making that part of the change, explain the concrete reason and scope, and ask for permission.

## Scientific Correctness

- Preserve numerical reproducibility. Do not intentionally change physical constants, equations, interpolation, integration, binning, normalization, solver behavior, default parameters, or operation ordering unless the task explicitly requires a scientifically meaningful change.
- Treat unexplained numerical drift as a regression, even when the build succeeds.
- Avoid changes that alter floating-point evaluation order unless they are necessary and reviewed for numerical impact.
- Preserve the formats, units, column meanings, ordering, precision, and parsing behavior of scientific tables.
- Preserve XML/XSW parameter names, defaults, units, numeric meanings, and output formats unless an explicit migration is requested.
- When a task intentionally changes numerical behavior, document the reason, affected observables, validation method, and expected tolerance in the final report.
- Do not infer scientific intent from a TODO comment alone. Ask when the expected physical behavior or acceptable numerical tolerance is unclear.

## Scientific Tables and Git LFS

- Do not modify, regenerate, normalize, reformat, replace, rename, or delete files under `bin/tables/`.
- Do not change Git LFS tracking rules or table pointer files without explicit user permission.
- Reading tables for diagnosis and verification is allowed.
- Confirm that a file is not a Git LFS pointer before treating its contents as the actual scientific dataset.
- Never use cleanup commands that could remove downloaded LFS objects or tracked runtime data.

## Configuration and Output Files

- Treat `.xsw` and `.xml` files as scientific configuration, not ordinary formatting targets. Preserve parameter values and numeric spellings unless the task requires a change.
- Do not overwrite canonical or user-provided configuration files when running experiments; make a clearly named untracked copy instead.
- Simulation runs create or replace result directories derived from the configuration filename. Inspect the selected output path before running.
- Keep generated results, local build products, IDE metadata, temporary patches, and scratch notes out of commits unless the user explicitly asks to version a specific artifact.
- Generated version files such as `_autogeneratedVersionInfo.h`, `_autogeneratedDif.c`, and `_autogeneratedGitInfo.c` are build products and should not be edited by hand.

## Working Tree and Git

- Inspect `git status --short --branch` before editing and again before reporting completion.
- Existing uncommitted or untracked files belong to the user. Preserve them and do not incorporate, overwrite, delete, or reformat them.
- Do not use destructive commands such as `git reset --hard`, `git clean`, or broad recursive deletion.
- Agents may create branches and commits when useful or requested. Use an `ai/` prefix for new branches unless the user provides another name.
- Use Conventional Commits for every commit message: `<type>[optional scope]: <description>`.
- Use an appropriate type such as `feat`, `fix`, `docs`, `build`, `test`, `refactor`, `perf`, `ci`, `chore`, or `revert`. Keep the description concise, imperative, and in English.
- Mark breaking changes with `!` before the colon and explain them in a `BREAKING CHANGE:` footer.
- Keep commits focused and exclude unrelated changes, generated files, local results, IDE files, and temporary Russian-language documents.
- Do not amend, rebase, force-push, or otherwise rewrite history unless explicitly requested.

## Completion Report

At the end of a task, report:

- What changed and why.
- Which build or verification commands were run and their results.
- Which supported platforms were actually verified.
- Any unverified platform, dependency limitation, numerical risk, or follow-up validation still needed.
- Any intentional numerical or data-format effect.

Do not claim successful cross-platform support based on a build performed on only one operating system.
