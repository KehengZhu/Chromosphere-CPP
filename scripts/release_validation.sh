#!/usr/bin/env bash
# model_column release validation.
#
# Two checks, both using the NORMAL release invocation (no ISO_RIEMANN /
# ISO_RECONSTRUCTION overrides — the SWMF exact-Riemann Godunov flux and
# (ln rho,V,ln p) are the defaults):
#
#   1. default-path equivalence — the no-override run must reproduce an explicit
#      ISO_RIEMANN=swmf-godunov ISO_RECONSTRUCTION=lnrho-v-lnp run bitwise (apart
#      from the self-referential sidecar path written into the snapshot header);
#   2. release smoke — N500 (661 actual cells), 100 s, with the face-flux and
#      outer-conduction sidecars, analysed by util/pressure_reconstruction_diag.py
#      for eps_p, Q_V, Q_M, Q_eff, F_eff, T_top, p_top and negative-velocity count.
#
# Usage:  scripts/release_validation.sh [output_dir]
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

OUT_DIR="${1:-outputs/model_column}"
TMP_DIR="tmp/release_validation"
PY="${PY:-.venv/bin/python}"
mkdir -p "${OUT_DIR}" "${TMP_DIR}"

run_column() {   # run_column <out> [extra env already exported by caller]
    CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT="${FRAME_DT:-5}" \
    scripts/run_chromo_realtime.sh "$1" \
        full no-ionization model_column - 20.0 no-cooling
}

echo "=== 1/2 default-path equivalence (T_END=20 s) ==="
CHROMO_T_END=20 run_column "${TMP_DIR}/default.txt" > "${TMP_DIR}/default.log" 2>&1
ISO_RIEMANN=swmf-godunov ISO_RECONSTRUCTION=lnrho-v-lnp \
CHROMO_T_END=20 run_column "${TMP_DIR}/explicit.txt" > "${TMP_DIR}/explicit.log" 2>&1

"${PY}" - "${TMP_DIR}" <<'PY'
import sys
d = sys.argv[1]
clean = lambda f: [l for l in open(f"{d}/{f}") if not l.startswith("# EOS_MODE=")]
same = clean("default.txt") == clean("explicit.txt")
print("snapshots identical:", same)
print("gamma_diag identical:",
      open(f"{d}/default.txt.gamma_diag","rb").read()
      == open(f"{d}/explicit.txt.gamma_diag","rb").read())
sys.exit(0 if same else 1)
PY

echo "=== 2/2 release smoke (N500, 100 s) ==="
SMOKE="${OUT_DIR}/release_godunov_lnp_N500_100s.txt"
CHROMO_T_END=100 FRAME_DT=10 \
CHROMO_FACE_FLUX_DIAG=1 CHROMO_OUTER_COND_DIAG=1 \
CHROMO_FACE_FLUX_STRIDE=1780 CHROMO_OUTER_COND_STRIDE=89 CHROMO_EOS_COUNTS=1 \
    run_column "${SMOKE}" > "${SMOKE%.txt}.console.log" 2>&1

# The sidecar header records the numerical method actually used.
grep -m1 '^# ns=' "${SMOKE}.faceflux"
grep -m1 '^termination=' "${SMOKE%.txt}.console.log" || true

# Release hygiene: the exact Riemann solve must have supplied EVERY face. A
# fallback to Rusanov is legal but is never expected on a healthy release run,
# and it would silently change the numerical method under the reported result.
grep -m1 '^godunov\.' "${SMOKE%.txt}.console.log" || true
grep -m1 '^godunov\.fallbacks=' "${SMOKE%.txt}.console.log"
if ! grep -q '^godunov\.fallbacks=0 ' "${SMOKE%.txt}.console.log"; then
    echo "FAIL: the release Godunov flux fell back to Rusanov on this run." >&2
    exit 1
fi

"${PY}" util/pressure_reconstruction_diag.py \
    --run release="${SMOKE}" --targets 100 \
    --json "${SMOKE%.txt}_metrics.json"

echo "=== RELEASE VALIDATION COMPLETE ==="
