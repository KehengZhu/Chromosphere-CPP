#!/usr/bin/env bash
# ARCHIVED — kept for provenance of outputs/model_column/acc_*, not for reuse.
#
# Controlled-variable acceptance matrix for the (ln rho,V,ln p) reconstruction,
# holding the Riemann solver fixed at Rusanov (the release solver AT THE TIME
# this was run) so the reconstruction change was the only variable. It produced
# the evidence that lnP is a reconstruction improvement rather than a
# Roe-specific patch; see docs/pressure_reconstruction_recap.md.
#
# The release now defaults to Roe + (ln rho,V,ln p), so ISO_RIEMANN=rusanov and
# ISO_RECONSTRUCTION below are reference overrides, not production settings.
# For a release check use scripts/release_validation.sh instead.
#
# Everything except ISO_NS (resolution axis) and ISO_RECONSTRUCTION (the
# variable under test) is held fixed.
#
# Diagnostic cadence is matched in PHYSICAL TIME across resolutions:
#   snapshots/face-flux : 20 s        (frame_dt=20 ; step strides below)
#   outer conduction    : ~0.5 s      (for a well-resolved int q dt)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

run_one() {
    local ns="$1" recon="$2" tag="$3" ff_stride="$4" oc_stride="$5"
    local out="outputs/model_column/acc_${tag}.txt"
    echo "=== ${tag}: ISO_NS=${ns} ISO_RECONSTRUCTION=${recon} ==="
    ISO_NS="${ns}" \
    ISO_RIEMANN=rusanov \
    ISO_RECONSTRUCTION="${recon}" \
    CHROMO_T_END=4000 \
    CHROMO_OUTPUT=1 \
    CHROMO_GAMMA_DIAG=1 \
    CHROMO_FACE_FLUX_DIAG=1 \
    CHROMO_OUTER_COND_DIAG=1 \
    CHROMO_FRAME_DT=20 \
    CHROMO_FACE_FLUX_STRIDE="${ff_stride}" \
    CHROMO_OUTER_COND_STRIDE="${oc_stride}" \
    CHROMO_EOS_COUNTS=1 \
    /usr/bin/time -p scripts/run_chromo_realtime.sh "${out}" \
        full no-ionization model_column - 20.0 no-cooling \
        > "outputs/model_column/acc_${tag}.console.log" 2>&1
    echo "--- ${tag} done ---"
}

# dt ~ 5.62e-3 s (N500) and ~2.805e-3 s (N1000) -> matched physical cadence.
run_one 500  lnrho-v-lnt lnt_N500  3560 89
run_one 500  lnrho-v-lnp lnp_N500  3560 89
run_one 1000 lnrho-v-lnt lnt_N1000 7130 178
run_one 1000 lnrho-v-lnp lnp_N1000 7130 178
echo "=== MATRIX COMPLETE ==="
