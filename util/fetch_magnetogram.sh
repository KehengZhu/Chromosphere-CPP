#!/usr/bin/env bash
# Download a GONG zero-point-corrected synoptic magnetogram into util/data/.
# Usage:
#   bash util/fetch_magnetogram.sh                     # default file
#   bash util/fetch_magnetogram.sh <fits-basename>     # override
#
# The default file (mrzqs190801t0014c2220_229.fits) is for Carrington
# rotation 2220, time 2019-08-01 00:14 UT, central Carrington longitude 229.
# GONG hosts these at https://gong2.nso.edu/oQR/zqs/yyyymm/mrzqsyymmdd/.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
mkdir -p "${DATA_DIR}"

BASENAME="${1:-mrzqs190801t0014c2220_229.fits}"
TARGET="${DATA_DIR}/${BASENAME}"

if [[ -s "${TARGET}" ]]; then
    echo "[fetch] already present: ${TARGET}"
    exit 0
fi

# Parse the date out of the filename (mrzqs YYMMDD ...).
DATESTR="${BASENAME#mrzqs}"
YYMMDD="${DATESTR:0:6}"
YYYYMM="20${YYMMDD:0:4}"
SUBDIR="mrzqs${YYMMDD}"

URLS=(
    "https://gong2.nso.edu/oQR/zqs/${YYYYMM}/${SUBDIR}/${BASENAME}.gz"
    "https://gong2.nso.edu/oQR/zqs/${YYYYMM}/${SUBDIR}/${BASENAME}"
    "https://nispdata.nso.edu/ftp/oQR/zqs/${YYYYMM}/${SUBDIR}/${BASENAME}.gz"
)

for url in "${URLS[@]}"; do
    echo "[fetch] trying ${url}"
    if curl -sSfL -o "${TARGET}.tmp" "${url}"; then
        if [[ "${url}" == *.gz ]]; then
            gunzip -c "${TARGET}.tmp" > "${TARGET}"
            rm -f "${TARGET}.tmp"
        else
            mv "${TARGET}.tmp" "${TARGET}"
        fi
        echo "[fetch] saved ${TARGET}"
        exit 0
    fi
    rm -f "${TARGET}.tmp"
done

echo "[fetch] ERROR: all candidate URLs failed for ${BASENAME}" >&2
echo "  GONG may have reorganized; check https://gong2.nso.edu/oQR/zqs/" >&2
exit 1
