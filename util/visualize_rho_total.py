#!/usr/bin/env python3
"""Visualize evolution of total mass density ρ_tot = ρ_i + ρ_n along the
field for a chromo_main run.

Produces (for the same run file):
  - <stem>_rho_total.png : space-time pcolormesh + snapshot lines + mass budget
  - <stem>_rho_total.mp4 : animation of ρ_tot(s) sweeping in time

Usage:
    python util/visualize_rho_total.py [output.txt] [--png ...] [--mp4 ...] [--fps 24]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colormaps

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plot_output import read_frames
from _anim_parallel import save_frames_parallel

CNI, CNN = 0, 1   # conserved-variable indices for ρ_i, ρ_n (see plot_output.py:20)


def load_rho_total(path):
    xx, frames = read_frames(path)
    t   = np.array([fr[0] for fr in frames])
    rho = np.stack([fr[2][:, CNI] + fr[2][:, CNN] for fr in frames])  # (nT, ns) kg/m^3
    return xx, t, rho


def render_png(out_png, xx, t, rho, ds_m_estimate):
    """Three-panel diagnostic: space-time map, snapshot lines, mass over time."""
    fig, axes = plt.subplots(
        3, 1, figsize=(11, 11),
        gridspec_kw=dict(height_ratios=[1.6, 1.0, 0.9], hspace=0.32),
    )

    # 1) Space-time pcolormesh (log) -----------------------------------------
    ax = axes[0]
    pos = rho > 0
    vmin = float(rho[pos].min())
    vmax = float(rho.max())
    pcm = ax.pcolormesh(
        xx, t, rho, shading="auto",
        norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
        cmap="viridis",
    )
    cb = fig.colorbar(pcm, ax=ax, pad=0.01)
    cb.set_label(r"$\rho_i + \rho_n$  (kg/m$^3$)")
    ax.set_xlabel("height s (km)")
    ax.set_ylabel("time (s)")
    ax.set_title(r"Space-time evolution of $\rho_\mathrm{tot} = \rho_i + \rho_n$")

    # 2) Snapshot lines ------------------------------------------------------
    ax = axes[1]
    n_show = 10
    idxs = np.linspace(0, len(t) - 1, n_show).astype(int)
    cmap_t = colormaps["viridis"]
    for k, i in enumerate(idxs):
        c = cmap_t(k / max(1, n_show - 1))
        ax.plot(xx, rho[i], color=c, lw=1.0)
    ax.plot(xx, rho[0],  "k:",  lw=1.6, label=f"t = {t[0]:.0f} s")
    ax.plot(xx, rho[-1], "C3-", lw=1.8, label=f"t = {t[-1]:.0f} s")
    ax.set_yscale("log")
    ax.set_xlabel("height s (km)")
    ax.set_ylabel(r"$\rho_\mathrm{tot}$  (kg/m$^3$)")
    ax.set_title(rf"$\rho_\mathrm{{tot}}(s)$ at {n_show} evenly-spaced times")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    sm = matplotlib.cm.ScalarMappable(
        cmap=cmap_t, norm=matplotlib.colors.Normalize(vmin=t[0], vmax=t[-1])
    )
    sm.set_array([])
    fig.colorbar(sm, ax=ax, pad=0.01, label="time (s)")

    # 3) Domain-integrated mass over time ------------------------------------
    # ds estimate: use mean cell spacing (ns cells, span = xx[-1] - xx[0] km).
    ds_km = (xx[-1] - xx[0]) / max(1, xx.size - 1)
    ds_m  = ds_km * 1.0e3
    M_per_area = rho.sum(axis=1) * ds_m  # kg/m^2 (column mass)
    ax = axes[2]
    ax.plot(t, M_per_area, lw=1.6, color="C0")
    ax.set_xlabel("time (s)")
    ax.set_ylabel(r"$\int \rho_\mathrm{tot}\,ds$  (kg/m$^2$)")
    ax.set_title("Domain-integrated mass column")
    ax.grid(True, alpha=0.3)
    drift = (M_per_area[-1] - M_per_area[0]) / max(M_per_area[0], 1e-30) * 100.0
    ax.text(
        0.98, 0.04,
        f"Δ(column mass) = {drift:+.3f}% over the run",
        transform=ax.transAxes, ha="right", va="bottom",
        fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
    )

    fig.suptitle(
        f"Total mass density evolution — {os.path.basename(out_png)[:-4]}",
        fontsize=12, y=0.995,
    )
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[viz] wrote {out_png}  (column-mass drift {drift:+.3f}%)")


def _render_rho_frame(k, tmpdir, xx, rho_k, rho_0, t_k, nT, vmin, vmax):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.set_xlim(xx.min(), xx.max())
    ax.set_ylim(vmin, vmax)
    ax.set_yscale("log")
    ax.set_xlabel("height s (km)")
    ax.set_ylabel(r"$\rho_i + \rho_n$  (kg/m$^3$)")
    ax.grid(True, alpha=0.3)
    ax.plot(xx, rho_0, lw=1.2, color="k", ls=":", alpha=0.7, label="ρ_tot(s, t=0)")
    ax.plot(xx, rho_k, lw=2.0, color="C3", label="ρ_tot(s, t)")
    ax.legend(loc="upper right", fontsize=9)
    ax.set_title(rf"$\rho_\mathrm{{tot}}(s)$  |  t = {t_k:7.2f} s   ({k+1}/{nT})")
    fig.tight_layout()
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=120)
    plt.close(fig)


def render_mp4(out_mp4, xx, t, rho, fps=24):
    pos = rho > 0
    vmin = float(rho[pos].min()) * 0.8
    vmax = float(rho.max()) * 1.2
    nT = len(t)
    save_frames_parallel(
        _render_rho_frame,
        [(k, xx, rho[k], rho[0], t[k], nT, vmin, vmax) for k in range(nT)],
        out_mp4, fps,
    )
    print(f"[viz] wrote {out_mp4}  ({nT} frames @ {fps} fps = {nT/fps:.1f} s)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", nargs="?", default="out_pfss_10x.txt",
                    help="path to chromo_main output.txt (default: out_pfss_10x.txt)")
    ap.add_argument("--png")
    ap.add_argument("--mp4")
    ap.add_argument("--fps", type=int, default=24)
    args = ap.parse_args()

    stem = os.path.splitext(os.path.basename(args.input))[0]
    viz_dir = os.path.join(HERE, "visualization")
    os.makedirs(viz_dir, exist_ok=True)
    out_png = args.png or os.path.join(viz_dir, f"{stem}_rho_total.png")
    out_mp4 = args.mp4 or os.path.join(viz_dir, f"{stem}_rho_total.mp4")

    print(f"[viz] reading {args.input}")
    xx, t, rho = load_rho_total(args.input)
    print(f"  ns={xx.size}, frames={t.size}, t in [{t[0]:.1f}, {t[-1]:.1f}] s")
    print(f"  rho_tot global range: [{rho[rho>0].min():.3e}, {rho.max():.3e}] kg/m^3")

    render_png(out_png, xx, t, rho, ds_m_estimate=None)
    render_mp4(out_mp4, xx, t, rho, fps=args.fps)


if __name__ == "__main__":
    main()
