#!/usr/bin/env bash
# Long-run / real-time production launcher.
#
# Same OpenMP defaults as run_chromo_omp.sh, plus the low-I/O output defaults a
# long CHROMO_T_END run needs: at ~7e-4 s timesteps the legacy time_mult-derived
# snapshot stride would emit tens of thousands of ASCII snapshots and dominate
# the wall time. Everything is an override-preserving default, so an explicit
# environment assignment always wins.
#
# Low-I/O production run (the shape used for real-time acceptance):
#   CHROMO_T_END=1000 scripts/run_chromo_realtime.sh out.txt \
#       full no-ionization model_column - 20.0 no-cooling
#
# Same run WITH snapshots, on a physical-time cadence of one frame per second:
#   CHROMO_OUTPUT=1 \
#   CHROMO_GAMMA_DIAG=1 \
#   CHROMO_FRAME_DT=1 \
#   scripts/run_chromo_realtime.sh out.txt \
#       full no-ionization model_column - 20.0 no-cooling
#
# Use scripts/run_chromo_omp.sh for the generic OpenMP runtime defaults.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
CHROMO_BINARY="${CHROMO_BINARY:-${PROJECT_ROOT}/build_omp/chromo_main}"

# OpenMP runtime (see CLAUDE.md: 12 threads = the performance cores).
export OMP_DYNAMIC="${OMP_DYNAMIC:-FALSE}"
export OMP_MAX_ACTIVE_LEVELS="${OMP_MAX_ACTIVE_LEVELS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-12}"

# Output defaults: everything expensive off, profiling on.
export CHROMO_OUTPUT="${CHROMO_OUTPUT:-0}"
export CHROMO_GAMMA_DIAG="${CHROMO_GAMMA_DIAG:-0}"
export CHROMO_FACE_FLUX_DIAG="${CHROMO_FACE_FLUX_DIAG:-0}"
export CHROMO_OUTER_COND_DIAG="${CHROMO_OUTER_COND_DIAG:-0}"
export CHROMO_EOS_COUNTS="${CHROMO_EOS_COUNTS:-0}"
export CHROMO_PROFILE="${CHROMO_PROFILE:-1}"
export CHROMO_PROGRESS_STRIDE="${CHROMO_PROGRESS_STRIDE:-10000}"

# Model/grid/physics defaults belong to the model_column scenario itself. This
# launcher adds runtime policy only. CFL=0.50 is validated for the canonical
# N=500/R4, 661-cell physical-conduction-only model_column release. Other
# scenarios retain the solver's conservative CHROMO_CFL=0.25 fallback unless
# CHROMO_CFL is explicitly supplied.
if [[ "${4:-}" == "model_column" ]]; then
    export CHROMO_CFL="${CHROMO_CFL:-0.50}"
fi

if [[ ! -x "${CHROMO_BINARY}" ]]; then
    echo "OpenMP chromo_main not found or not executable: ${CHROMO_BINARY}" >&2
    echo "Build it with:" >&2
    echo "  cmake -S \"${PROJECT_ROOT}\" -B \"${PROJECT_ROOT}/build_omp\" -DCMAKE_BUILD_TYPE=Release -DCHROMO_ENABLE_OPENMP=ON" >&2
    echo "  cmake --build \"${PROJECT_ROOT}/build_omp\" -j4" >&2
    echo "Or point CHROMO_BINARY at an existing OpenMP build." >&2
    exit 1
fi

exec "${CHROMO_BINARY}" "$@"
