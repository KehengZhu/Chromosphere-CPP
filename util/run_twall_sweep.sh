#!/usr/bin/env bash
# Upper-boundary conduction-wall temperature sweep for the model_column release.
#
# Everything except the conduction wall is held fixed: the same C7 initial
# condition, the same N=500/R4 mesh, the release Godunov flux and reconstruction,
# physical conduction only, and ISO_HYDRO_T_DECOUPLE=1 so the HYDRO ghost
# temperature still free-floats with the live top cell. ISO_T_TOP therefore moves
# exactly one thing: grid.outer_conduction_temperature, the fixed thermal datum
# the implicit conduction stage sees at the PHYSICAL outer face.
#
# The double-precision release build is mandatory here (build_omp, i.e.
# CHROMO_STATE_FLOAT32=OFF): the float32 state quantizes the lower-chromosphere
# mass update and would contaminate a 4000 s mass-flux/velocity comparison.
#
# CFL: 0.50 is the validated production value for this configuration and is used
# for every case that tolerates it. The implicit conduction Newton iteration does
# not converge at CFL 0.50 for a wall far above the top-cell temperature (fails at
# 100 kK, marginal at 0.35, converges at 0.30), so cases listed in SLOW_CFL_CASES
# run at the conservative solver default 0.25 instead. util/run_twall_sweep.sh
# also runs a CFL control at one wall temperature so the two groups can be joined.
#
# Usage:  util/run_twall_sweep.sh [out_dir] [T_END] [FRAME_DT]
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

OUT_DIR="${1:-outputs/model_column/twall_sweep}"
T_END="${2:-4000}"
FRAME_DT="${3:-20}"
mkdir -p "${OUT_DIR}"

# T_wall [K] : CFL. Override with e.g. SWEEP_CASES="11000:0.50 12000:0.50" to run
# a refinement pass without re-running the broad sweep.
CASES=(
    "10000:0.50"
    "15000:0.50"
    "18000:0.50"
    "22000:0.50"
    "30000:0.50"
    "45000:0.50"
    "70000:0.50"
    "100000:0.25"
    "70000:0.25"     # CFL control: joins the 0.50 group to the 0.25 group
)
if [[ -n "${SWEEP_CASES:-}" ]]; then
    read -r -a CASES <<< "${SWEEP_CASES}"
fi

# 8-9 concurrent jobs at 2 OpenMP threads each. Measured on this workstation the
# solver reaches ~1.1e3 steps/s at 2 threads versus ~1.4e3 at 12, so concurrent
# 2-thread jobs are far better throughput than sequential 12-thread jobs. Every
# case uses the SAME thread count so the comparison is apples-to-apples.
THREADS="${SWEEP_THREADS:-2}"

# Conservative face mass flux f_total (<out>.faceflux), OFF by default: it is a
# read-only diagnostic that changes no numerical result (verified byte-identical
# .gamma_diag with and without it), but it costs ~230 kB per capture, so a 10 s
# cadence over 4000 s is ~90 MB per run and ~1.2 GB for the whole sweep. Turn it
# on with SWEEP_FACE_FLUX=1 when the sweep needs the conservative flux rather than
# the cell-centred product rho*V — e.g. for the overlay movie's mass-flux panel.
# The default stride 1780 is ~10 s at the CFL 0.50 step and matches the value
# scripts/release_validation.sh uses.
FACE_FLUX="${SWEEP_FACE_FLUX:-0}"
FACE_FLUX_STRIDE="${SWEEP_FACE_FLUX_STRIDE:-1780}"

pids=()
for case in "${CASES[@]}"; do
    T="${case%%:*}"
    CFL="${case##*:}"
    tag="twall_$((T/1000))kK_cfl${CFL}"
    out="${OUT_DIR}/${tag}.txt"
    (
        OMP_NUM_THREADS="${THREADS}" \
        CHROMO_CFL="${CFL}" \
        ISO_T_TOP="${T}" \
        CHROMO_T_END="${T_END}" \
        CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT="${FRAME_DT}" \
        CHROMO_OUTER_COND_DIAG=1 CHROMO_OUTER_COND_STRIDE=2000 \
        CHROMO_FACE_FLUX_DIAG="${FACE_FLUX}" \
        CHROMO_FACE_FLUX_STRIDE="${FACE_FLUX_STRIDE}" \
        scripts/run_chromo_realtime.sh "${out}" \
            full no-ionization model_column - 20.0 no-cooling \
            > "${out%.txt}.console.log" 2>&1
        echo "${tag} exit=$?"
    ) &
    pids+=($!)
done

status=0
for pid in "${pids[@]}"; do
    wait "${pid}" || status=1
done

echo "=== sweep summary ==="
for case in "${CASES[@]}"; do
    T="${case%%:*}"; CFL="${case##*:}"
    log="${OUT_DIR}/twall_$((T/1000))kK_cfl${CFL}.console.log"
    printf '%-24s %s | %s | %s\n' "$((T/1000))kK cfl=${CFL}" \
        "$(grep -m1 '^termination=' "${log}" || echo 'termination=MISSING')" \
        "$(grep -m1 '^final_time=' "${log}" || echo 'final_time=?')" \
        "$(grep -m1 '^godunov\.fallbacks=' "${log}" | cut -d' ' -f1 || echo 'fallbacks=?')"
done
exit "${status}"
