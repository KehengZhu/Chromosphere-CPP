#!/usr/bin/env bash
# Create a project-local Python venv at .venv/ and install pipeline deps.
# Idempotent: re-running just upgrades the deps.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
VENV="${REPO_ROOT}/.venv"

if [[ ! -d "${VENV}" ]]; then
    echo "[setup_venv] creating venv at ${VENV}"
    python3 -m venv "${VENV}"
fi

# shellcheck source=/dev/null
source "${VENV}/bin/activate"
python -m pip install --upgrade pip
# Try sunkit-magex first (current name); fall back to pfsspy if that fails.
pip install -r "${SCRIPT_DIR}/requirements.txt" || {
    echo "[setup_venv] full install failed; trying pfsspy-only fallback"
    pip install numpy scipy astropy sunpy pfsspy
}

echo "[setup_venv] done. activate with:  source ${VENV}/bin/activate"
