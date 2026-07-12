#!/usr/bin/env python3
"""Compare ρ_tot = ρ_i + ρ_n evolution between two chromo_main runs
(e.g. ionization on vs off, or different scenarios).

Usage:
    python util/visualize_rho_total_compare.py <runA.txt> <runB.txt> \\
        [--labelA "ionization ON"] [--labelB "ionization OFF"] \\
        [--out-png util/<...>.png] [--out-mp4 util/<...>.mp4] [--fps 24]

Defaults compare the 10x PFSS runs that we just produced.
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

CNI, CNN = 0, 1


def load(path):
    xx, frames = read_frames(path)
    t   = np.array([fr[0] for fr in frames])
    rho = np.stack([fr[2][:, CNI] + fr[2][:, CNN] for fr in frames])
    return xx, t, rho


def render_png(out_png, xx, tA, rA, tB, rB, labelA, labelB):
    """4-panel layout: two pcolormeshes on top, difference + column-mass below."""
    assert xx.size == rA.shape[1] == rB.shape[1], "grid mismatch"

    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.0, 1.05], hspace=0.32, wspace=0.30)

    # Common log scale spanning both runs.
    pos_min = float(min(rA[rA > 0].min(), rB[rB > 0].min()))
    pos_max = float(max(rA.max(), rB.max()))
    norm = matplotlib.colors.LogNorm(vmin=pos_min, vmax=pos_max)

    # --- Pcolormesh A ---
    axA = fig.add_subplot(gs[0, 0])
    pcmA = axA.pcolormesh(xx, tA, rA, shading="auto", norm=norm, cmap="viridis")
    axA.set_xlabel("height s (km)"); axA.set_ylabel("time (s)")
    axA.set_title(rf"$\rho_\mathrm{{tot}}$ — {labelA}")

    # --- Pcolormesh B ---
    axB = fig.add_subplot(gs[0, 1], sharey=axA)
    pcmB = axB.pcolormesh(xx, tB, rB, shading="auto", norm=norm, cmap="viridis")
    axB.set_xlabel("height s (km)")
    axB.set_title(rf"$\rho_\mathrm{{tot}}$ — {labelB}")
    cbar = fig.colorbar(pcmB, ax=[axA, axB], pad=0.02)
    cbar.set_label(r"$\rho_i + \rho_n$  (kg/m$^3$)")

    # --- Snapshot lines at start, mid, end for both runs ---
    axS = fig.add_subplot(gs[0, 2])
    show_idxsA = [0, len(tA) // 2, -1]
    show_idxsB = [0, len(tB) // 2, -1]
    styles = [("k:", 1.3), ("C0-", 1.4), ("C3-", 1.7)]
    labelsA = [f"{labelA}, t=0", f"{labelA}, t={tA[len(tA)//2]:.0f}", f"{labelA}, t={tA[-1]:.0f}"]
    labelsB = [f"{labelB}, t=0", f"{labelB}, t={tB[len(tB)//2]:.0f}", f"{labelB}, t={tB[-1]:.0f}"]
    for k, (i, (sty, lw), lab) in enumerate(zip(show_idxsA, styles, labelsA)):
        axS.plot(xx, rA[i], sty, lw=lw, label=lab)
    # plot B as same colours but dashed
    for k, (i, (sty, lw), lab) in enumerate(zip(show_idxsB, styles, labelsB)):
        color = "k" if k == 0 else ("C0" if k == 1 else "C3")
        axS.plot(xx, rB[i], color=color, ls="--", lw=lw, alpha=0.85, label=lab)
    axS.set_yscale("log")
    axS.set_xlabel("height s (km)"); axS.set_ylabel(r"$\rho_\mathrm{tot}$  (kg/m$^3$)")
    axS.set_title("snapshots (solid = A, dashed = B)")
    axS.grid(True, alpha=0.3)
    axS.legend(loc="best", fontsize=7)

    # --- Relative difference at run B's frame timing (resampled to A) ---
    axD = fig.add_subplot(gs[1, 0:2], sharex=axA)
    # Resample B onto A's time axis (both ran with similar dt -> safe nearest).
    idxBforA = np.array([int(np.argmin(np.abs(tB - t))) for t in tA])
    rB_at_A  = rB[idxBforA]
    rel = (rA - rB_at_A) / np.maximum(rB_at_A, 1e-30) * 100.0  # %
    vmax = float(np.nanmax(np.abs(rel)))
    pcmD = axD.pcolormesh(
        xx, tA, rel, shading="auto",
        norm=matplotlib.colors.Normalize(vmin=-vmax, vmax=vmax),
        cmap="RdBu_r",
    )
    cbD = fig.colorbar(pcmD, ax=axD, pad=0.02)
    cbD.set_label(rf"$(\rho_{{A}} - \rho_{{B}}) / \rho_{{B}}$  (%) — A={labelA}, B={labelB}")
    axD.set_xlabel("height s (km)"); axD.set_ylabel("time (s)")
    axD.set_title("Relative difference  A − B  (red = A heavier)")

    # --- Domain-integrated column mass over time ---
    axM = fig.add_subplot(gs[1, 2])
    ds_km = (xx[-1] - xx[0]) / max(1, xx.size - 1)
    ds_m  = ds_km * 1.0e3
    MA = rA.sum(axis=1) * ds_m  # kg/m^2
    MB = rB.sum(axis=1) * ds_m
    axM.plot(tA, MA, "C0-", lw=1.6, label=labelA)
    axM.plot(tB, MB, "C3--", lw=1.6, label=labelB)
    axM.set_xlabel("time (s)")
    axM.set_ylabel(r"$\int \rho_\mathrm{tot}\,ds$  (kg/m$^2$)")
    axM.set_title("column mass")
    axM.grid(True, alpha=0.3)
    axM.legend(loc="best", fontsize=9)
    drA = (MA[-1] - MA[0]) / max(MA[0], 1e-30) * 100.0
    drB = (MB[-1] - MB[0]) / max(MB[0], 1e-30) * 100.0
    txt = (f"Δcolumn ({labelA}):  {drA:+.2f}%\n"
           f"Δcolumn ({labelB}):  {drB:+.2f}%")
    axM.text(0.98, 0.04, txt, transform=axM.transAxes, ha="right", va="bottom",
             fontsize=9, bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85))

    fig.suptitle(rf"$\rho_\mathrm{{tot}} = \rho_i + \rho_n$ evolution: {labelA}  vs  {labelB}",
                 fontsize=13, y=0.995)
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[viz] wrote {out_png}")
    print(f"  column-mass drift: A {drA:+.2f}%,  B {drB:+.2f}%")
    return rB_at_A


def _render_rho_compare_frame(k, tmpdir, xx, rA_k, rB_k, rA_0, t_A, labelA, labelB, nT, vmin, vmax):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.set_xlim(xx.min(), xx.max())
    ax.set_ylim(vmin, vmax)
    ax.set_yscale("log")
    ax.set_xlabel("height s (km)")
    ax.set_ylabel(r"$\rho_i + \rho_n$  (kg/m$^3$)")
    ax.grid(True, alpha=0.3)
    ax.plot(xx, rA_0, lw=1.0, color="k", ls=":", alpha=0.6, label="t = 0 reference (both)")
    ax.plot(xx, rA_k, lw=2.0, color="C0", label=labelA)
    ax.plot(xx, rB_k, lw=2.0, color="C3", ls="--", label=labelB)
    ax.legend(loc="upper right", fontsize=9)
    ax.set_title(
        rf"$\rho_\mathrm{{tot}}(s)$  |  t = {t_A:7.2f} s   "
        f"({k+1}/{nT})    [{labelA} solid, {labelB} dashed]"
    )
    fig.tight_layout()
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=120)
    plt.close(fig)


def render_mp4(out_mp4, xx, tA, rA, tB, rB, labelA, labelB, fps=24):
    idxBforA = np.array([int(np.argmin(np.abs(tB - t))) for t in tA])
    vmin = float(min(rA[rA > 0].min(), rB[rB > 0].min())) * 0.8
    vmax = float(max(rA.max(), rB.max())) * 1.2
    nT = len(tA)
    save_frames_parallel(
        _render_rho_compare_frame,
        [(k, xx, rA[k], rB[idxBforA[k]], rA[0], tA[k], labelA, labelB, nT, vmin, vmax)
         for k in range(nT)],
        out_mp4, fps,
    )
    print(f"[viz] wrote {out_mp4}  ({nT} frames @ {fps} fps = {nT/fps:.1f} s)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runA", nargs="?", default="out_pfss_10x.txt")
    ap.add_argument("runB", nargs="?", default="out_pfss_10x_noioniz.txt")
    ap.add_argument("--labelA", default="ionization ON")
    ap.add_argument("--labelB", default="ionization OFF")
    ap.add_argument("--out-png", default=None)
    ap.add_argument("--out-mp4", default=None)
    ap.add_argument("--fps", type=int, default=24)
    args = ap.parse_args()

    viz_dir = os.path.join(HERE, "..", "visualization")
    os.makedirs(viz_dir, exist_ok=True)
    out_png = args.out_png or os.path.join(viz_dir, "pfss_10x_rho_total_compare.png")
    out_mp4 = args.out_mp4 or os.path.join(viz_dir, "pfss_10x_rho_total_compare.mp4")

    print(f"[viz] reading A: {args.runA}")
    xxA, tA, rA = load(args.runA)
    print(f"  ns={xxA.size}, frames={tA.size}, t in [{tA[0]:.1f}, {tA[-1]:.1f}] s")
    print(f"[viz] reading B: {args.runB}")
    xxB, tB, rB = load(args.runB)
    print(f"  ns={xxB.size}, frames={tB.size}, t in [{tB[0]:.1f}, {tB[-1]:.1f}] s")
    assert xxA.size == xxB.size and np.allclose(xxA, xxB), "grid mismatch between runs"

    render_png(out_png, xxA, tA, rA, tB, rB, args.labelA, args.labelB)
    render_mp4(out_mp4, xxA, tA, rA, tB, rB, args.labelA, args.labelB, fps=args.fps)


if __name__ == "__main__":
    main()
