#!/usr/bin/env python3
"""Relaxation / steady-state convergence monitor for a gentle-evaporation run.

Reads a chromo_main multi-snapshot output .txt (the same format plot_output.py
reads) and tracks, frame by frame, the diagnostics that decide whether the
atmosphere has reached the V≈0 steady preflare state (docs/gentle_evaporation_plan.md
Phase 2): max|V|(t), the density-weighted mean |V|(t), the per-cell temperature
drift max|∂T/∂t|/T, and the apex / coronal-mean temperature. The fixed point is
declared when max|V| sits far below the expected evaporation speed (~tens–hundreds
of km/s) and the profiles stop drifting.

This is the time-series companion to check_loop_ic.py (which checks a single static
.dat IC for hydrostatic balance before the run).

Usage:
  python util/check_relaxation.py outputs/model_column/gentle_relax.txt [visualization/model_column/gentle_relax_convergence.png]
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M_I = 1.6726219e-27
M_N = M_I
K_B = 1.380649e-23
CNI, CNN, CNV, CNU, CEI, CEN = 0, 1, 2, 3, 4, 5

# A relaxed preflare loop should drift far slower than the gentle-evaporation
# upflow it will later drive (tens of km/s). Declare convergence when the bulk
# flow is below this and the temperature is no longer drifting appreciably.
V_STEADY_MS = 2.0e3          # density-weighted mean |V| target [m/s]
DTDT_STEADY = 3.0e-3         # max fractional temperature drift target [1/s]


def read_frames(path):
    """Return (arc_km, frames) with frames = list of (t, step, xn[ns,neq])."""
    with open(path) as f:
        raw = [ln for ln in f.read().splitlines() if ln.strip()]
    hdr = raw[0].split()
    ns, neq = int(hdr[0]), int(hdr[1])
    arc_km = np.fromstring(raw[1], sep=" ")
    frames, i = [], 2
    while i < len(raw):
        line = raw[i]
        if not line.startswith("#"):
            raise ValueError(f"expected '# t = ...' at line {i+1}, got: {line!r}")
        toks = line.split()
        t, stp = float(toks[3]), int(toks[6])
        i += 1
        if i + ns > len(raw):
            break  # truncated final frame
        data = np.zeros((ns, neq))
        for j in range(ns):
            data[j, :] = np.fromstring(raw[i + j], sep=" ")
        i += ns
        frames.append((t, stp, data))
    return arc_km, frames


def primitives(xn):
    ni = xn[:, CNI] / M_I
    nn = xn[:, CNN] / M_N
    v  = xn[:, CNV] / (M_I * ni)
    u  = xn[:, CNU] / (M_N * nn)
    p_i = 2.0 / 3.0 * xn[:, CEI] - 1.0 / 3.0 * M_I * ni * v * v
    Ti = p_i / (2.0 * ni * K_B)
    rho = xn[:, CNI] + xn[:, CNN]
    return ni, nn, v, u, Ti, rho


def main():
    in_path = sys.argv[1] if len(sys.argv) > 1 else "outputs/model_column/gentle_relax.txt"
    out = sys.argv[2] if len(sys.argv) > 2 else "visualization/model_column/gentle_relax_convergence.png"
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    arc_km, frames = read_frames(in_path)
    if len(frames) < 2:
        sys.exit(f"need ≥2 frames to measure drift; found {len(frames)}")

    t      = np.array([fr[0] for fr in frames])
    vmax   = np.zeros(len(frames))
    vmean  = np.zeros(len(frames))     # density-weighted mean |V|
    Tapex  = np.zeros(len(frames))     # max T_i (loop apex / hottest cell)
    Tcor   = np.zeros(len(frames))     # mean T_i over coronal cells (T>0.5 MK)
    dTdt   = np.full(len(frames), np.nan)   # max fractional |dT/dt| vs previous frame

    Ti_prev, t_prev = None, None
    for k, (tk, _, xn) in enumerate(frames):
        ni, nn, v, u, Ti, rho = primitives(xn)
        av = np.abs(v)
        vmax[k]  = av.max()
        vmean[k] = np.sum(rho * av) / np.sum(rho)
        Tapex[k] = Ti.max()
        cor = Ti > 5.0e5
        Tcor[k] = Ti[cor].mean() if cor.any() else np.nan
        if Ti_prev is not None and tk > t_prev:
            dTdt[k] = np.max(np.abs(Ti - Ti_prev) / np.maximum(Ti, 1.0)) / (tk - t_prev)
        Ti_prev, t_prev = Ti, tk

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(f"Gentle-evaporation relaxation — {os.path.basename(in_path)} "
                 f"({len(frames)} frames, t→{t[-1]:.0f} s)", fontsize=13)

    ax[0, 0].semilogy(t, np.maximum(vmax, 1e-3), "C0", label="max |V|")
    ax[0, 0].semilogy(t, np.maximum(vmean, 1e-3), "C1", label="ρ-weighted mean |V|")
    ax[0, 0].axhline(V_STEADY_MS, color="k", ls="--", lw=0.8,
                     label=f"steady target {V_STEADY_MS/1e3:.0f} km/s")
    ax[0, 0].set_ylabel("|V| [m/s]"); ax[0, 0].set_xlabel("t [s]")
    ax[0, 0].legend(fontsize=8); ax[0, 0].grid(alpha=0.3)
    ax[0, 0].set_title("bulk-flow decay (steady ⇒ both fall and flatten)")

    ax[0, 1].semilogy(t, np.maximum(dTdt, 1e-6), "C3")
    ax[0, 1].axhline(DTDT_STEADY, color="k", ls="--", lw=0.8,
                     label=f"steady target {DTDT_STEADY:.0e} /s")
    ax[0, 1].set_ylabel(r"max$_s$ |∂lnT/∂t| [1/s]"); ax[0, 1].set_xlabel("t [s]")
    ax[0, 1].legend(fontsize=8); ax[0, 1].grid(alpha=0.3)
    ax[0, 1].set_title("temperature drift (→0 at fixed point)")

    ax[1, 0].plot(t, Tapex / 1e6, "C2", label="apex T (max)")
    ax[1, 0].plot(t, Tcor / 1e6, "C4", label="coronal mean T (>0.5 MK)")
    ax[1, 0].set_ylabel("T [MK]"); ax[1, 0].set_xlabel("t [s]")
    ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=0.3)
    ax[1, 0].set_title("coronal temperature (settles at the relaxed T_max)")

    # Final-frame field-aligned |V| profile vs arc length.
    _, _, vL, _, TiL, _ = primitives(frames[-1][2])
    ax[1, 1].plot(arc_km / 1e3, vL / 1e3, "C0")
    ax[1, 1].axhline(0, color="k", lw=0.6)
    ax[1, 1].set_ylabel("V [km/s]"); ax[1, 1].set_xlabel("arc length [Mm]")
    ax[1, 1].grid(alpha=0.3)
    ax[1, 1].set_title(f"final |V| profile (t = {t[-1]:.0f} s)")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")
    print(f"  frames: {len(frames)}  t: {t[0]:.1f} → {t[-1]:.1f} s")
    print(f"  final max|V|       = {vmax[-1]:.3e} m/s")
    print(f"  final mean|V| (ρw) = {vmean[-1]:.3e} m/s   (target < {V_STEADY_MS:.0e})")
    print(f"  final max|dlnT/dt| = {dTdt[-1]:.3e} /s     (target < {DTDT_STEADY:.0e})")
    print(f"  final apex T       = {Tapex[-1]/1e6:.3f} MK")
    print(f"  final coronal⟨T⟩   = {Tcor[-1]/1e6:.3f} MK")
    converged = (vmean[-1] < V_STEADY_MS) and (dTdt[-1] < DTDT_STEADY)
    print(f"  CONVERGED: {converged}")


if __name__ == "__main__":
    main()
