#!/usr/bin/env python3
"""Space-time / propagation figures for the model_column low-chromosphere mode.

Panel set 1 (`--out ..._spacetime.png`): V(h, t) heatmaps for two runs plus the
high-cadence early transient taken from a `.faceflux` sidecar, which resolves the
conduction-launched downgoing acoustic front and its reflection off the inner
V = 0 wall.

Panel set 2 (`--out ..._modes.png`): V(t) at one height for every supplied run,
the mode amplitude profile vs height, and the conservative evaporation-window
mass flux.

Usage:
    python util/plot_column_sloshing.py \
        --run LABEL=<run>.txt.gamma_diag [--run ...] \
        [--faceflux <run>.txt.faceflux] [--height 1800] \
        --out visualization/model_column/roe_n1000_sloshing.png
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from column_spacetime import load, col

FACEFLUX_COLS = ("cell_km face_km rho_cell v_cell T_cell rho_L rho_R v_L v_R T_L T_R "
                 "cs_L cs_R a_face f_central f_diff f_total r_rho "
                 "phi_plus_rho r_ip1_rho phi_minus_rho r_v phi_plus_v r_T phi_plus_T "
                 "p_cell p_L p_R").split()


def read_faceflux_v(path, t_max):
    """-> (cell_km, t, V) for records with t <= t_max (streamed)."""
    iv, ik = FACEFLUX_COLS.index("v_cell"), FACEFLUX_COLS.index("cell_km")
    ts, vs, km, cur, ct = [], [], None, None, None
    with open(path) as fh:
        for line in fh:
            if line[0] == "#":
                if line.startswith("# t "):
                    if cur:
                        a = np.asarray(cur, dtype=float)
                        km = a[:, ik] if km is None else km
                        ts.append(ct)
                        vs.append(a[:, iv])
                    ct = float(line.split()[3])
                    if ct > t_max:
                        break
                    cur = []
                continue
            if cur is not None:
                cur.append([float(x) for x in line.split()])
    return km, np.asarray(ts), np.asarray(vs)


def detrend(t, x, win=1400.0):
    dt = t[1] - t[0]
    w = int(round(win / dt)) | 1
    return x - np.convolve(np.pad(x, (w // 2, w // 2), mode="edge"),
                           np.ones(w) / w, mode="valid")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, metavar="LABEL=PATH")
    ap.add_argument("--faceflux", default=None)
    ap.add_argument("--height", type=float, default=1800.0)
    ap.add_argument("--evap", nargs=2, type=float, default=(2130.0, 2150.0))
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    runs = []
    for spec in args.run:
        label, _, path = spec.partition("=")
        runs.append((label, *load(path)))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    ncol = 3 if args.faceflux else 2
    fig, ax = plt.subplots(2, ncol, figsize=(5.2 * ncol, 8.4))

    # --- row 0: space-time heatmaps of the first two runs ------------------
    vmax = 0.0
    maps = []
    for label, h, t, d in runs[:2]:
        v = col(d, "v")
        m = (h >= 1600) & (h <= 2130)     # exclude the TR boundary layer
        maps.append((label, h[m], t, v[:, m]))
        vmax = max(vmax, np.abs(v[:, m]).max())
    for k, (label, hh, tt, vv) in enumerate(maps):
        a = ax[0, k]
        im = a.pcolormesh(tt, hh, vv.T, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                          shading="auto")
        a.set_title(f"{label}:  V(h, t)")
        a.set_xlabel("t [s]")
        a.set_ylabel("height [km]")
        fig.colorbar(im, ax=a, label="V [m/s]")

    # --- row 0 col 2: high-cadence early transient ------------------------
    if args.faceflux:
        km, te, ve = read_faceflux_v(args.faceflux, 220.0)
        m = km <= 2130
        a = ax[0, 2]
        lim = np.abs(ve[:, m]).max()
        im = a.pcolormesh(te, km[m], ve[:, m].T, cmap="RdBu_r",
                          vmin=-lim, vmax=lim, shading="auto")
        a.set_title("conduction turn-on: downgoing front + reflection\n"
                    "(high-cadence .faceflux capture)")
        a.set_xlabel("t [s]")
        a.set_ylabel("height [km]")
        fig.colorbar(im, ax=a, label="V [m/s]")

    # --- row 1 col 0: V(t) at the probe height ----------------------------
    a = ax[1, 0]
    for label, h, t, d in runs:
        i = int(np.argmin(np.abs(h - args.height)))
        a.plot(t, col(d, "v")[:, i], lw=1.2, label=label)
    a.axhline(0, c="k", lw=0.6)
    a.set_xlabel("t [s]")
    a.set_ylabel("V [m/s]")
    a.set_title(f"V at {args.height:.0f} km")
    a.legend(fontsize=7)
    a.grid(alpha=0.3)

    # --- row 1 col 1: mode amplitude vs height ----------------------------
    a = ax[1, 1]
    for label, h, t, d in runs:
        m = t >= 300
        amp = np.array([detrend(t, col(d, "v")[:, i])[m].std()
                        for i in range(h.size)])
        sel = h <= 2152
        a.plot(h[sel], amp[sel], lw=1.2, label=label)
    a.set_xlabel("height [km]")
    a.set_ylabel("rms of high-pass V [m/s]")
    a.set_title("mode amplitude vs height (t $\\geq$ 300 s)")
    a.legend(fontsize=7)
    a.grid(alpha=0.3)

    # --- row 1 col 2: evaporation-window mass flux ------------------------
    if ncol == 3:
        a = ax[1, 2]
        for label, h, t, d in runs:
            m = (h >= args.evap[0]) & (h <= args.evap[1])
            a.plot(t, (col(d, "rho") * col(d, "v"))[:, m].mean(1), lw=1.2, label=label)
        a.axhline(0, c="k", lw=0.6)
        a.set_xlabel("t [s]")
        a.set_ylabel(r"$\langle \rho V \rangle$ [kg m$^{-2}$ s$^{-1}$]")
        a.set_title(f"evaporation window {args.evap[0]:.0f}-{args.evap[1]:.0f} km")
        a.legend(fontsize=7)
        a.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(args.out, dpi=140)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
