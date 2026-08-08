#!/usr/bin/env python3
"""Three-panel (rho, V, p) evolution movie for a Gamma/Saha model_column run.

Plots the physical fields straight out of the <run>.gamma_diag sidecar
(columns: rho_total v_cm T x_eq n_e n_HI p_total Gamma1 kappa_phys kappa_solver),
so density, velocity and total gas pressure are the EOS-consistent x_eq-derived
fields rather than a two-fluid reconstruction.

The model_column domain (~1600-2153 km) is thin enough that a linear height axis
is honest -- no split/compressed axis is needed (that convention applies to runs
spanning the chromosphere through a multi-Mm corona). The 2130-2150 km analysis
window used by the release diagnostics is shaded for reference.

Usage:
  python util/animate_column_rvp.py <run.txt> <out.mp4> [fps]
Env:
  ANIM_MAX_FRAMES  frame budget (default: all frames)
  ANIM_STRIDE      take every Nth frame (default 1)
"""
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _anim_parallel import save_frames_parallel

WINDOW = (2130.0, 2150.0)          # release velocity/roughness analysis window
ZOOM = (2100.0, 2153.0)            # refined top / evaporation region
C_RHO = "#2f3640"                  # neutral
C_UP = "#1f6fb2"                   # signal blue  (upflow)
C_DN = "#c0392b"                   # accent red   (downflow)
C_P = "#2f3640"


def _apply_style():
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "figure.dpi": 130,
    })


def read_frames(path):
    """Stream the gamma_diag sidecar into (heights, [(t, array), ...])."""
    with open(path) as f:
        raw = [ln for ln in f.read().splitlines() if ln.strip()]
    ns, _ncol = (int(x) for x in raw[0].split())
    h_km = np.fromstring(raw[1], sep=" ")
    frames, i = [], 2
    while i < len(raw):
        if not raw[i].startswith("# t ="):
            i += 1
            continue
        t = float(raw[i].split()[3])
        i += 1
        if i + ns > len(raw):
            break
        d = np.array([np.fromstring(raw[i + j], sep=" ") for j in range(ns)])
        i += ns
        frames.append((t, d))
    return h_km, frames


def _render(k, tmpdir, h_km, rho, v, p, t, lims, zlims, label):
    import matplotlib.pyplot as plt

    _apply_style()
    fig, axes = plt.subplots(3, 2, figsize=(9.6, 7.2), sharex="col",
                             gridspec_kw={"width_ratios": [1.45, 1.0]})

    for col, (xlo, xhi, ylims) in enumerate(
            ((h_km[0], h_km[-1], lims), (ZOOM[0], ZOOM[1], zlims))):
        ax_r, ax_v, ax_p = axes[0, col], axes[1, col], axes[2, col]
        (rho_lo, rho_hi), (v_lo, v_hi), (p_lo, p_hi) = ylims

        for ax in (ax_r, ax_v, ax_p):
            ax.axvspan(*WINDOW, color="#f0c419", alpha=0.13, lw=0, zorder=0)
            ax.set_xlim(xlo, xhi)

        ax_r.semilogy(h_km, rho, color=C_RHO, lw=1.1)
        ax_r.set_ylim(rho_lo, rho_hi)

        ax_v.axhline(0.0, color="#95a5a6", lw=0.6, zorder=1)
        ax_v.fill_between(h_km, 0.0, v, where=(v >= 0), color=C_UP, alpha=0.22, lw=0)
        ax_v.fill_between(h_km, 0.0, v, where=(v < 0), color=C_DN, alpha=0.22, lw=0)
        ax_v.plot(h_km, v, color=C_UP, lw=1.1)
        ax_v.set_ylim(v_lo, v_hi)

        ax_p.semilogy(h_km, p, color=C_P, lw=1.1)
        ax_p.set_ylim(p_lo, p_hi)
        ax_p.set_xlabel("height  [km]")

    axes[0, 0].set_ylabel(r"$\rho$  [kg m$^{-3}$]")
    axes[1, 0].set_ylabel(r"$V$  [m s$^{-1}$]")
    axes[2, 0].set_ylabel(r"$p_{\rm total}$  [Pa]")

    axes[0, 0].set_title(f"{label}   full column", fontsize=9, loc="left",
                         color="#57606f")
    axes[0, 1].set_title(f"zoom {ZOOM[0]:.0f}-{ZOOM[1]:.0f} km   "
                         f"(shaded = {WINDOW[0]:.0f}-{WINDOW[1]:.0f} km window)",
                         fontsize=8, loc="left", color="#57606f")
    axes[0, 0].text(0.985, 0.92, f"t = {t:7.1f} s", transform=axes[0, 0].transAxes,
                    ha="right", va="top", fontsize=12, color="#2f3640",
                    family="monospace")

    fig.tight_layout(h_pad=0.7, w_pad=1.4)
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"))
    plt.close(fig)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    run, out = sys.argv[1], sys.argv[2]
    fps = int(sys.argv[3]) if len(sys.argv) > 3 else 20

    diag = run + ".gamma_diag"
    if not os.path.exists(diag):
        sys.exit(f"missing gamma_diag sidecar: {diag}")
    h_km, frames = read_frames(diag)

    stride = int(os.environ.get("ANIM_STRIDE", 1))
    frames = frames[::stride]
    budget = os.environ.get("ANIM_MAX_FRAMES")
    if budget and len(frames) > int(budget):
        sel = np.linspace(0, len(frames) - 1, int(budget)).astype(int)
        frames = [frames[i] for i in sel]
    print(f"{len(frames)} frames, t = {frames[0][0]:.1f} -> {frames[-1][0]:.1f} s, "
          f"{len(h_km)} cells")

    rho = [d[:, 0] for _, d in frames]
    vel = [d[:, 1] for _, d in frames]
    prs = [d[:, 6] for _, d in frames]

    def lim(arrs, log, pad=0.05):
        lo = min(a.min() for a in arrs)
        hi = max(a.max() for a in arrs)
        if log:
            return lo * 0.7, hi * 1.4
        span = hi - lo
        return lo - pad * span, hi + pad * span

    lims = (lim(rho, True), lim(vel, False), lim(prs, True))
    # The zoom column gets its own limits, taken over the zoom region only, so the
    # late-time evaporation structure is not squashed by the early top transient.
    zm = (h_km >= ZOOM[0]) & (h_km <= ZOOM[1])
    zlims = (lim([a[zm] for a in rho], True),
             lim([a[zm] for a in vel], False),
             lim([a[zm] for a in prs], True))
    label = os.path.basename(run)

    save_frames_parallel(
        _render,
        [(k, h_km, rho[k], vel[k], prs[k], frames[k][0], lims, zlims, label)
         for k in range(len(frames))],
        out, fps,
    )


if __name__ == "__main__":
    main()
