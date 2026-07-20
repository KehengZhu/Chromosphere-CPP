#!/bin/bash
# Build the CRASH-EOS gamma tabulator against the prebuilt SWMF CRASH library
# and run it.  Must be able to see the SWMF symlink at the project root.
#
#   Prereq: the CRASH library must be built once:
#       cd SWMF/util/CRASH/src && make LIB
#
set -euo pipefail

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
./util/eos/tabulate_gamma.exe

echo "output:"
ls -l outputs/eos_gamma/
