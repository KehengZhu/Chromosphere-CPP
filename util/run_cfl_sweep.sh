#!/usr/bin/env bash
# CFL sweep driver for the production-shaped model_column case.
#   util/run_cfl_sweep.sh perf  <cfl> [<cfl> ...]   # 3 timed low-I/O reps each
#   util/run_cfl_sweep.sh acc   <cfl> [<cfl> ...]   # 1 run each, 1 frame/s
#   util/run_cfl_sweep.sh med   <t_end> <frame_dt> <cfl> [<cfl> ...]
# Outputs land in outputs/model_column/cfl_sweep/.
set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"
EXE="${CHROMO_BINARY:-${ROOT}/build_cfl_audit/chromo_main}"
OUTDIR="${ROOT}/outputs/model_column/cfl_sweep"
mkdir -p "${OUTDIR}"

# Production-shaped model configuration. Identical for every run in the sweep;
# only CHROMO_CFL (and the output controls) change.
model_env=(
  OMP_DYNAMIC=FALSE OMP_MAX_ACTIVE_LEVELS=1 OMP_NUM_THREADS="${OMP_NUM_THREADS:-12}"
  GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat
  ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_NS=2000
  ISO_HYDRO_T_DECOUPLE=1 ISO_REFINE_PROFILE=outer ISO_REFINE_FACTOR=4
  ISO_REFINE_S_LO_KM=500 ISO_REFINE_TRANSITION_KM=20
  ISO_NUMERICAL_DIFFUSIVITY_MULT=1
  CHROMO_FACE_FLUX_DIAG=0 CHROMO_OUTER_COND_DIAG=0 CHROMO_PROFILE=1
)

tag() { printf "%s" "$(echo "$1" | tr -d '.')"; }

run_one() {  # $1=out_stem  rest=extra env
  local stem="$1"; shift
  /usr/bin/time -p env "${model_env[@]}" "$@" \
    "${EXE}" "${stem}.txt" full no-ionization model_column - 20.0 no-cooling \
    > "${stem}.log" 2>&1
}

mode="$1"; shift
case "${mode}" in
perf)
  for cfl in "$@"; do
    for rep in 1 2 3; do
      stem="${OUTDIR}/perf_cfl$(tag "${cfl}")_rep${rep}"
      run_one "${stem}" CHROMO_CFL="${cfl}" CHROMO_T_END=20 \
        CHROMO_OUTPUT=0 CHROMO_GAMMA_DIAG=0 CHROMO_EOS_COUNTS=0
      w=$(grep -E '^real' "${stem}.log" | awk '{print $2}')
      s=$(grep -E '^final_step=' "${stem}.log" | cut -d= -f2)
      echo "perf cfl=${cfl} rep=${rep} wall=${w} steps=${s}"
    done
  done
  ;;
acc)
  for cfl in "$@"; do
    stem="${OUTDIR}/cfl$(tag "${cfl}")_20s"
    run_one "${stem}" CHROMO_CFL="${cfl}" CHROMO_T_END=20 \
      CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=1 CHROMO_EOS_COUNTS=1
    echo "acc cfl=${cfl} $(grep -E '^(real|termination|final_step)' "${stem}.log" | tr '\n' ' ')"
  done
  ;;
med)
  tend="$1"; fdt="$2"; shift 2
  for cfl in "$@"; do
    stem="${OUTDIR}/med${tend}s_cfl$(tag "${cfl}")"
    run_one "${stem}" CHROMO_CFL="${cfl}" CHROMO_T_END="${tend}" \
      CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT="${fdt}" CHROMO_EOS_COUNTS=1
    echo "med cfl=${cfl} $(grep -E '^(real|termination|final_step)' "${stem}.log" | tr '\n' ' ')"
  done
  ;;
*) echo "unknown mode: ${mode}" >&2; exit 1;;
esac
