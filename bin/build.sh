#!/bin/bash
# Custom build script for Mac OS

set -euo pipefail

if [[ $# -ne 1 || ( "$1" != "Debug" && "$1" != "Release" ) ]]
then
    echo "usage: build.sh Debug | Release"
    exit 1
fi

find_tool()
{
    local tool_name
    for tool_name in "$@"
    do
        if command -v "$tool_name" >/dev/null 2>&1
        then
            command -v "$tool_name"
            return 0
        fi
    done
    return 1
}

build_config="$1"
script_dir="$(cd "$(dirname "$0")" && pwd)"
source_dir="$script_dir/../src"
build_dir="$source_dir/$build_config"

c_compiler="$(find_tool gcc-15 gcc-14 gcc)" || {
    echo "error: GCC was not found" >&2
    exit 1
}
cxx_compiler="$(find_tool g++-15 g++-14 g++)" || {
    echo "error: G++ was not found" >&2
    exit 1
}
fortran_compiler="$(find_tool gfortran-15 gfortran-14 gfortran)" || {
    echo "error: GFortran was not found" >&2
    exit 1
}

echo "Building $build_config with:"
echo "  C:       $c_compiler"
echo "  C++:     $cxx_compiler"
echo "  Fortran: $fortran_compiler"

cmake --fresh \
    -S "$source_dir" \
    -B "$build_dir" \
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
    -DCMAKE_C_COMPILER="$c_compiler" \
    -DCMAKE_CXX_COMPILER="$cxx_compiler" \
    -DCMAKE_Fortran_COMPILER="$fortran_compiler" \
    -DCMAKE_BUILD_TYPE="$build_config"

cmake --build "$build_dir" --parallel

echo "Build output: $build_dir"
