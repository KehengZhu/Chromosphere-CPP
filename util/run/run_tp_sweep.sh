#!/usr/bin/env bash
# Boundary (T_wall, p_top) sweep of the model_column release, for comparing the
# simulated evaporation mass flux with the analytic steady-TR closure
# (docs/studies/boundaries/analytic_evaporation_J_table.md).
#
# Each case is the canonical release (N=500/R4, single fixed outer temperature,
# physical conduction only, CFL 0.50) with two things changed:
#   CHROMO_BC_T_TOP_K  -- the outer ghost / conduction temperature;
#   CHROMO_IC_P_TOP_PA -- the IC is rebuilt with the same C7 T(h) and the base pressure
#                 that puts the top cell at this pressure, so both reservoirs
#                 move consistently with it (unset = the C7 column).
# One N=1000 control at a mid-grid point bounds the resolution error.
#
# Jobs run concurrently at SWEEP_THREADS OpenMP threads each through xargs -P,
# which works with the bash 3.2 that macOS ships (no `wait -n`). Default 4 jobs x
# 2 threads = the 8 performance cores of an M2 Pro; set SWEEP_JOBS/SWEEP_THREADS
# for other hardware. Cases are started most expensive first.
#
# Usage:  util/run/run_tp_sweep.sh [out_dir] [T_END] [FRAME_DT]
#         SWEEP_CASES="T_K:p_Pa:NS ..." util/run/run_tp_sweep.sh ...   # subset
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

if [[ "${1:-}" == "--one" ]]; then
    # Worker: --one OUT_DIR T_END FRAME_DT T_K:p_Pa:NS
    OUT_DIR="$2"; T_END="$3"; FRAME_DT="$4"; case_="$5"
    T="${case_%%:*}"; rest="${case_#*:}"; P="${rest%%:*}"; NS="${rest##*:}"
    tag="tp_T$((T/1000))kK_p${P}_N${NS}"
    out="${OUT_DIR}/${tag}.txt"
    if grep -qs '^termination=end_time' "${out%.txt}.console.log"; then
        echo "${tag} already complete, skipped"; exit 0
    fi
    # Face-flux cadence ~10 s of physical time at the CFL-0.50 step of the 22 kK
    # column; hotter walls take shorter steps, so their cadence is finer.
    stride=$(( 1780 * NS / 500 ))
    set +e
    OMP_NUM_THREADS="${SWEEP_THREADS:-2}" \
    CHROMO_BC_T_TOP_K="${T}" CHROMO_IC_P_TOP_PA="${P}" CHROMO_MESH_NS="${NS}" \
    CHROMO_T_END="${T_END}" \
    CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT="${FRAME_DT}" \
    CHROMO_OUTER_COND_DIAG=1 CHROMO_OUTER_COND_STRIDE=2000 \
    CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_STRIDE="${stride}" \
    scripts/run_chromo_realtime.sh "${out}" \
        full no-ionization model_column - 20.0 no-cooling \
        > "${out%.txt}.console.log" 2>&1
    echo "${tag} exit=$?"
    exit 0
fi

OUT_DIR="${1:-outputs/model_column/tp_sweep}"
T_END="${2:-3000}"
FRAME_DT="${3:-10}"
mkdir -p "${OUT_DIR}"

# 5 wall temperatures x 6 pressures. p spans the AWSoM-R transition-region top
# pressure for T_TR ~ 3e5-1e6 K at L_TR ~ 1.4 Mm; 0.0103 Pa is the C7 column.
# T stays <= 70 kK, where CFL 0.50 was validated in the T_wall sweep.
T_LIST="70000 50000 35000 25000 18000"
P_LIST="0.08 0.04 0.02 0.0103 0.005 0.003"
CASES=("50000:0.0103:1000")                     # resolution control, longest first
for T in ${T_LIST}; do for P in ${P_LIST}; do CASES+=("${T}:${P}:500"); done; done
if [[ -n "${SWEEP_CASES:-}" ]]; then
    read -r -a CASES <<< "${SWEEP_CASES}"
fi

echo "sweep: ${#CASES[@]} cases, T_END=${T_END} s, ${SWEEP_JOBS:-4} jobs x ${SWEEP_THREADS:-2} threads -> ${OUT_DIR}"
printf '%s\n' "${CASES[@]}" | \
    xargs -P "${SWEEP_JOBS:-4}" -I{} "$0" --one "${OUT_DIR}" "${T_END}" "${FRAME_DT}" {}

echo "=== sweep summary ==="
status=0
for case_ in "${CASES[@]}"; do
    T="${case_%%:*}"; rest="${case_#*:}"; P="${rest%%:*}"; NS="${rest##*:}"
    log="${OUT_DIR}/tp_T$((T/1000))kK_p${P}_N${NS}.console.log"
    term="$(grep -m1 '^termination=' "${log}" 2>/dev/null || echo 'termination=MISSING')"
    [[ "${term}" == "termination=end_time" ]] || status=1
    printf '%-28s %s | %s | %s\n' "T=$((T/1000))kK p=${P} N=${NS}" "${term}" \
        "$(grep -m1 '^final_physical_time=' "${log}" 2>/dev/null || echo 'final_physical_time=?')" \
        "$(grep -m1 '^godunov\.fallbacks=' "${log}" 2>/dev/null | cut -d' ' -f1 || echo 'fallbacks=?')"
done
exit "${status}"
