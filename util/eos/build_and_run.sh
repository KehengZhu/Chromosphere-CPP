#!/bin/bash
# Build the CRASH-EOS gamma tabulator against the prebuilt SWMF CRASH library
# and run it.  Must be able to see the SWMF symlink at the project root.
#
#   Prereq: the CRASH library must be built once:
#       cd SWMF/util/CRASH/src && make LIB
#
set -euo pipefail

MODE="${1:-production}"
if [[ "$MODE" != "production" && "$MODE" != "excitation" ]]; then
    echo "usage: bash util/eos/build_and_run.sh [production|excitation]" >&2
    exit 2
fi

HERE="$(cd "$(dirname "$0")" && pwd)"      # util/eos
ROOT="$(cd "$HERE/../.." && pwd)"          # project root
SWMF="$ROOT/SWMF"                           # project symlink to the SWMF checkout
INC="$SWMF/share/include"
LIBDIR="$SWMF/lib"

if [[ ! -d "$SWMF/util/CRASH/src" ]]; then
    echo "error: $SWMF is missing or does not point to a valid SWMF checkout" >&2
    exit 1
fi

cd "$ROOT"
mkdir -p outputs/eos_gamma data/eos

# reals/doubles MUST be 8 bytes to match the library build
FFLAGS="-cpp -w -frecord-marker=4 -fdefault-real-8 -fdefault-double-8 -O3 -I$INC"

echo "compiling tabulate_gamma.f90 ..."
gfortran $FFLAGS -c util/eos/tabulate_gamma.f90 -o util/eos/tabulate_gamma.o

echo "extracting library objects ..."
TMP="$(mktemp -d)"
( cd "$TMP" && ar -x "$LIBDIR/libCRASH.a" && ar -x "$LIBDIR/libTIMING.a" && ar -x "$LIBDIR/libSHARE.a" )

echo "linking ..."
mpif90 -fdefault-real-8 -fdefault-double-8 -frecord-marker=4 \
    -o util/eos/tabulate_gamma.exe util/eos/tabulate_gamma.o "$TMP"/*.o -lc++
rm -rf "$TMP"

echo "running ..."
if [[ "$MODE" == "excitation" ]]; then
    ./util/eos/tabulate_gamma.exe excitation
    echo "output:"
    ls -l outputs/eos_gamma/gamma_hydrogen_excitation.dat \
          outputs/eos_gamma/gamma1_hydrogen_excitation.dat
    exit 0
fi

./util/eos/tabulate_gamma.exe

echo "running 2x-refined Stage-9 grid ..."
./util/eos/tabulate_gamma.exe refined

echo "validating raw and runtime tables ..."
cmake -S "$ROOT" -B "$ROOT/build"
cmake --build "$ROOT/build" --target eos_validate_saha eos_validate_refined_table -j4
"$ROOT/build/eos_validate_saha" \
    outputs/eos_gamma/gamma_hydrogen.dat data/eos/gamma1_hydrogen_v1.dat
"$ROOT/build/eos_validate_refined_table" \
    data/eos/gamma1_hydrogen_v1.dat \
    outputs/eos_gamma/gamma1_hydrogen_refined.dat \
    outputs/eos_gamma/gamma_hydrogen_refined.dat

echo "checking tracked runtime-table checksum ..."
cmake -DROOT="$ROOT" -P "$ROOT/cmake/check_eos_checksum.cmake"

if [[ -n "${PYTHON:-}" ]]; then
    PYTHON_BIN="$PYTHON"
elif [[ -x "$ROOT/.venv/bin/python" ]]; then
    PYTHON_BIN="$ROOT/.venv/bin/python"
else
    PYTHON_BIN=python3
fi
"$PYTHON_BIN" util/plot_eos_gamma.py --validate-only

echo "output:"
ls -l outputs/eos_gamma/gamma_hydrogen.dat data/eos/gamma1_hydrogen_v1.dat \
      outputs/eos_gamma/gamma_hydrogen_refined.dat \
      outputs/eos_gamma/gamma1_hydrogen_refined.dat
