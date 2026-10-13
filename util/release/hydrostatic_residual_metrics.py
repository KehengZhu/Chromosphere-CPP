#!/usr/bin/env python3
"""Hydrostatic-residual metrics of conduction-off model_column runs.

For each run prefix (the snapshot path given to the driver) reads the
`.faceflux` and `.gamma_diag` sidecars and prints:

  * t = 0 max|F^rho| [kg m^-2 s^-1] over all faces, its face height, and the
    base-face split f_central / f_diff / f_total;
  * the last-frame max|V| [m/s], its cell (by upper-face height, as the
    .gamma_diag header lists it), and the maximum excluding the top
    cell.

Run with CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_TOP=0 CHROMO_GAMMA_DIAG=1 so
both sidecars cover every cell. Used by
docs/studies/numerics/hydrostatic_reconstruction_and_mc3_fix.md.

Usage:
  util/release/hydrostatic_residual_metrics.py RUN.txt [RUN.txt ...]
"""
from __future__ import annotations

import sys

import numpy as np

# .faceflux column indices (see docs/doxygen/pages/io_formats.md)
FACE_KM, F_CENTRAL, F_DIFF, F_TOTAL = 1, 14, 15, 16


def frames(path: str) -> list[tuple[float, np.ndarray]]:
    """All `# t = ...` frames of a sidecar as (time [s], rows) pairs."""
    out, rows, t = [], None, None
    with open(path) as handle:
        for line in handle:
            if line.startswith("# t = "):
                if rows is not None:
                    out.append((t, np.array(rows, float)))
                t, rows = float(line.split()[3]), []
            elif rows is not None and not line.startswith("#") and line.strip():
                rows.append(line.split())
    if rows:
        out.append((t, np.array(rows, float)))
    return out


def main(runs: list[str]) -> None:
    for run in runs:
        _, face = frames(run + ".faceflux")[0]
        flux = np.abs(face[:, F_TOTAL])
        k = int(np.argmax(flux))
        base = face[0]
        with open(run + ".gamma_diag") as handle:
            handle.readline()
            heights = np.array(handle.readline().split(), float)
        t_last, cells = frames(run + ".gamma_diag")[-1]
        speed = np.abs(cells[:, 1])
        kv = int(np.argmax(speed))
        kvi = int(np.argmax(speed[:-1]))
        print(f"{run}\n"
              f"  t=0   max|F^rho| = {flux[k]:.4e} at face {face[k, FACE_KM]:.1f} km; "
              f"base face f_central {base[F_CENTRAL]:+.3e} f_diff {base[F_DIFF]:+.3e} "
              f"f_total {base[F_TOTAL]:+.3e}\n"
              f"  t={t_last:<5.1f} max|V| = {speed[kv]:.4e} m/s at {heights[kv]:.1f} km "
              f"(excluding top cell {speed[kvi]:.4e} at {heights[kvi]:.1f} km)")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
