#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
CHROMO_BINARY="${CHROMO_BINARY:-${PROJECT_ROOT}/build_omp/chromo_main}"

export OMP_DYNAMIC="${OMP_DYNAMIC:-FALSE}"
export OMP_MAX_ACTIVE_LEVELS="${OMP_MAX_ACTIVE_LEVELS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-12}"

if [[ ! -x "${CHROMO_BINARY}" ]]; then
    echo "OpenMP chromo_main not found or not executable: ${CHROMO_BINARY}" >&2
    echo "Build it with:" >&2
    echo "  cmake -S \"${PROJECT_ROOT}\" -B \"${PROJECT_ROOT}/build_omp\" -DCMAKE_BUILD_TYPE=Release -DCHROMO_ENABLE_OPENMP=ON" >&2
    echo "  cmake --build \"${PROJECT_ROOT}/build_omp\" -j4" >&2
    exit 1
fi

exec "${CHROMO_BINARY}" "$@"
