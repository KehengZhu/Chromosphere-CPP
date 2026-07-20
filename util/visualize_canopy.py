#!/usr/bin/env python3
"""Visualize the analytic_canopy scenario: 3D flux-tube geometry + time
evolution of physical quantities along the field.

Canopy recipe (analytic, no magnetogram needed):
    B(z) = B_inf + (B_0 - B_inf) * exp(-(z - z_base) / H_B)

with defaults from scenarios/analytic_canopy.cpp:
    B_0   = 100 G  = 1.0e-2 T   (footpoint network field)
    B_inf =  15 G  = 1.5e-3 T   (canopy-merged value)
    H_B   = 300 km                (magnetic scale height)

By flux conservation B(z) * A(z) = const, so the tube cross-section grows
as A(z) = A_0 * B_0 / B(z). The 3D "trumpet" we draw is this axisymmetric
flux tube; field lines decorate its outer surface.

Outputs:
    util/canopy_geometry_evolution.png  (static, BC + 3D + evolution panels)
    util/canopy_evolution.mp4           (animation of the evolution panels)

Usage (from repo root with .venv active):
    python util/visualize_canopy.py
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import colormaps

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plot_output import read_frames, primitives
from _anim_parallel import save_frames_parallel


# Match scenarios/analytic_canopy.cpp:
B0_T   = 1.0e-2     # 100 G footpoint
BINF_T = 1.5e-3     # 15 G canopy
HB_M   = 3.0e5      # 300 km magnetic scale height
G_SI   = 274.0
H_BASE_KM = 1003.0  # C7 photospheric base height


def B_canopy(z_m):
    """B(z) along the field. z is height above the C7 base (m)."""
    return BINF_T + (B0_T - BINF_T) * np.exp(-z_m / HB_M)


def draw_flux_tube(ax, z_m, r0_Mm=0.5, n_phi=72, alpha=0.30, cmap_name="plasma"):
    """Draw an axisymmetric flux tube whose cross-sectional radius follows
    flux conservation A(z) ∝ 1/B(z) -> r(z) = r0 * sqrt(B_0 / B(z)).

    z_m is in metres (height above base), r0_Mm sets the photospheric tube
    radius in Mm.

    Tube is rendered as a surface coloured by B(z), and field lines run along
    constant azimuth on the surface.
    """
    z_Mm = z_m / 1.0e6
    Bz = B_canopy(z_m)
    rz_Mm = r0_Mm * np.sqrt(B0_T / Bz)

    phi = np.linspace(0.0, 2.0 * np.pi, n_phi)
    PHI, ZZ = np.meshgrid(phi, z_Mm)
    RR     = np.broadcast_to(rz_Mm[:, None], PHI.shape)
    Xs = RR * np.cos(PHI)
    Ys = RR * np.sin(PHI)
    Zs = ZZ

    norm = matplotlib.colors.LogNorm(
        vmin=float(min(BINF_T, B0_T) * 0.9),
        vmax=float(max(BINF_T, B0_T) * 1.1),
    )
    cmap = colormaps[cmap_name]
    facecolors = cmap(norm(np.broadcast_to(Bz[:, None] * np.ones_like(PHI[0]), PHI.shape)))
    surf = ax.plot_surface(
        Xs, Ys, Zs,
        facecolors=facecolors,
        linewidth=0, antialiased=False, shade=False, alpha=alpha,
        rcount=len(z_Mm), ccount=n_phi, zorder=2,
    )

    # Field lines on the tube surface, every ~36 deg
    for phi0 in np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False):
        ax.plot(rz_Mm * np.cos(phi0), rz_Mm * np.sin(phi0), z_Mm,
                color="C3", lw=1.4, alpha=0.9, zorder=5)

    # A baseline disk at z=0 with the footpoint radius
    th = np.linspace(0.0, 2.0 * np.pi, 200)
    ax.plot(r0_Mm * np.cos(th), r0_Mm * np.sin(th), np.zeros_like(th),
            color="0.3", lw=1.0, zorder=6)

    return norm, cmap


def render_combined(args, xx_sim, frames_sim):
    """Static figure: 3D tube (top-left), B(z) profile (top-mid), BC table
    (top-right), and 4-panel evolution (bottom)."""
    z_m = np.linspace(0.0, args.top_height_km * 1000.0, 200)
    z_Mm = z_m / 1.0e6
    Bz_G = B_canopy(z_m) * 1.0e4

    prims = [primitives(fr[2]) for fr in frames_sim]
    t_arr = np.array([fr[0] for fr in frames_sim])
    ni = np.stack([p[0] for p in prims])
    nn = np.stack([p[1] for p in prims])
    v  = np.stack([p[2] for p in prims])
    Ti = np.stack([p[6] for p in prims])
    f  = ni / (ni + nn)

    fig = plt.figure(figsize=(19, 11))
    gs = gridspec.GridSpec(
        3, 4,
        height_ratios=[1.5, 1.0, 0.4],
        width_ratios=[1.4, 1.0, 1.0, 1.0],
        hspace=0.42, wspace=0.32,
        left=0.04, right=0.98, top=0.94, bottom=0.06,
    )

    # ---- 3D flux tube --------------------------------------------------------
    ax3d = fig.add_subplot(gs[0:2, 0], projection="3d")
    ax3d.set_box_aspect((1.0, 1.0, 1.4))
    norm, cmap = draw_flux_tube(ax3d, z_m, r0_Mm=0.5)
    # Axes
    ax3d.set_xlabel("x (Mm)"); ax3d.set_ylabel("y (Mm)")
    ax3d.set_zlabel("z above base (Mm)")
    ax3d.set_title(
        "Analytic-canopy flux tube  |  "
        rf"B(z) = $B_\infty + (B_0 - B_\infty)\,e^{{-z/H_B}}$"
        "\n"
        f"B_0 = {B0_T*1e4:.0f} G,  B_∞ = {BINF_T*1e4:.0f} G,  H_B = {HB_M/1e3:.0f} km",
        fontsize=11)
    ax3d.view_init(elev=18, azim=-55)
    # Colour bar for B(z)
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax = fig.add_axes([0.01, 0.30, 0.012, 0.40])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("B(z) (T) on tube surface")

    # ---- B(z) profile --------------------------------------------------------
    axB = fig.add_subplot(gs[0, 1])
    axB.plot(Bz_G, z_Mm, color="C3", lw=2.0)
    axB.set_xlabel("|B| (G)"); axB.set_ylabel("z above base (Mm)")
    axB.set_title("B(z) — exponential canopy", fontsize=11)
    axB.axhline(HB_M / 1e6, color="0.6", ls="--", lw=0.8, label=f"z = H_B = {HB_M/1e3:.0f} km")
    axB.legend(loc="best", fontsize=9)
    axB.grid(True, alpha=0.3)

    # ---- A(z) cross-section --------------------------------------------------
    axA = fig.add_subplot(gs[0, 2])
    A_rel = B0_T / B_canopy(z_m)
    axA.plot(A_rel, z_Mm, color="C2", lw=2.0)
    axA.set_xlabel("A(z) / A_0  (cross-section, normalized)")
    axA.set_ylabel("z above base (Mm)")
    axA.set_title("Tube cross-section\n(flux conservation: A · B = const)", fontsize=11)
    axA.grid(True, alpha=0.3)

    # ---- BC summary table ----------------------------------------------------
    axt = fig.add_subplot(gs[0, 3])
    axt.axis("off")
    # the canopy BC mirrors model_c7's pinned-Dirichlet + doubled T_e
    base_T = 6225.0    # C7 base
    top_T  = 6674.0    # C7 top
    rows = [
        ["",                              "B (G)",         "n_e (m⁻³)",      "n_n (m⁻³)",      "T_n (K)",     "T_e (K)"],
        ["photospheric base (z = 0)",     f"{B0_T*1e4:.0f}",   "1.90e+17",  "2.69e+19",     f"{base_T:.0f}", f"{base_T:.0f}"],
        ["chromospheric top (z = top)",   f"{BINF_T*1e4:.0f} (asymp.)",   "4.83e+16",  "3.45e+16",     f"{top_T:.0f}",  f"{2*top_T:.0f}  (×2 mimic-corona)"],
    ]
    tbl = axt.table(cellText=rows, loc="center", cellLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(9)
    tbl.scale(1, 1.6)
    for j in range(len(rows[0])):
        tbl[(0, j)].set_facecolor("#dddddd")
        tbl[(0, j)].set_text_props(weight="bold")
    axt.set_title("BC summary  (Dirichlet, halved-velocity damping; same as model_c7)", fontsize=10)

    # ---- Evolution panels (row 1, four columns) ------------------------------
    cmap_t = colormaps["viridis"]
    n_show = 12
    idxs = np.linspace(0, len(t_arr) - 1, n_show).astype(int)
    panels = [("f", f, "ionization fraction f", False),
              ("Ti", Ti, r"$T_i$ (K)", False),
              ("V",  v,  r"$V$ (m/s, ion velocity)", False),
              ("nn", nn, r"$n_n$ (m$^{-3}$)", True)]
    for j, (key, arr, label, logy) in enumerate(panels):
        ax = fig.add_subplot(gs[1, j])
        for k, idx in enumerate(idxs):
            color = cmap_t(k / max(1, n_show - 1))
            ax.plot(xx_sim, arr[idx], color=color, lw=0.9, alpha=0.95)
        ax.plot(xx_sim, arr[0],  color="black", lw=1.4, ls=":",  alpha=0.95, label="t=0 s")
        ax.plot(xx_sim, arr[-1], color="C3",    lw=1.6, alpha=0.95, label=f"t={t_arr[-1]:.0f} s")
        ax.grid(True, alpha=0.3)
        ax.set_xlabel("height s (km)")
        ax.set_ylabel(label)
        if logy:
            ax.set_yscale("log")
        if key == "V":
            ax.axhline(0, color="0.7", lw=0.5)
        if j == 0:
            ax.legend(loc="best", fontsize=9)
    # Time colour-bar on the far right
    sm_t = matplotlib.cm.ScalarMappable(cmap=cmap_t, norm=matplotlib.colors.Normalize(vmin=0, vmax=t_arr[-1]))
    sm_t.set_array([])
    cax_t = fig.add_axes([0.985, 0.30, 0.008, 0.18])
    cbar_t = fig.colorbar(sm_t, cax=cax_t)
    cbar_t.set_label("simulation time (s)", fontsize=9)

    fig.suptitle(
        "Analytic-canopy scenario — flux-tube geometry, BC, and chromospheric relaxation",
        fontsize=13, y=0.985,
    )
    fig.savefig(args.out_png, dpi=140, bbox_inches="tight")
    print(f"[viz] wrote {args.out_png}")
    plt.close(fig)


def _render_canopy_frame(k, tmpdir, xx_sim, f_k, Ti_k, v_k, nn_k,
                         f_0, Ti_0, Bz_G, z_km, A_rel, t_k, nT, ylims_bot):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 4, figsize=(15, 7))
    fig.suptitle(
        f"Analytic-canopy scenario  |  t = {t_k:6.2f} s   step {k}/{nT-1}",
        fontsize=12,
    )

    # Top row — static geometry + t=0 reference
    axes[0, 0].plot(Bz_G, z_km, "C3-", lw=1.8)
    axes[0, 0].set_xlabel("|B| (G)"); axes[0, 0].set_ylabel("z above base (km)")
    axes[0, 0].set_title("B(z)"); axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].plot(A_rel, z_km, "C2-", lw=1.8)
    axes[0, 1].set_xlabel("A(z)/A_0"); axes[0, 1].set_title("Cross-section")
    axes[0, 1].grid(True, alpha=0.3)

    axes[0, 2].plot(xx_sim, f_0, "k:", lw=1.4, label="t=0")
    axes[0, 2].plot(xx_sim, f_k, color="C3", lw=1.4)
    axes[0, 2].set_xlabel("height s (km)"); axes[0, 2].set_ylabel("f")
    axes[0, 2].set_title("ionization fraction f")
    axes[0, 2].grid(True, alpha=0.3); axes[0, 2].legend(loc="best", fontsize=9)

    axes[0, 3].plot(xx_sim, Ti_0, "k:", lw=1.4, label="t=0")
    axes[0, 3].plot(xx_sim, Ti_k, color="C3", lw=1.4)
    axes[0, 3].set_xlabel("height s (km)"); axes[0, 3].set_ylabel(r"$T_i$ (K)")
    axes[0, 3].set_title("ion temperature")
    axes[0, 3].grid(True, alpha=0.3); axes[0, 3].legend(loc="best", fontsize=9)

    # Bottom row — animated quantities
    bot_panels = [
        (f_k,  "ionization fraction f", False),
        (Ti_k, r"$T_i$ (K)",            False),
        (v_k,  r"$V$ (m/s)",            False),
        (nn_k, r"$n_n$ (m$^{-3}$)",     True),
    ]
    for j, (arr, label, logy) in enumerate(bot_panels):
        ax = axes[1, j]
        ax.set_xlim(xx_sim.min(), xx_sim.max())
        ax.set_ylim(*ylims_bot[j])
        if logy:
            ax.set_yscale("log")
        if label == r"$V$ (m/s)":
            ax.axhline(0, color="0.7", lw=0.5)
        ax.set_xlabel("height s (km)"); ax.set_ylabel(label); ax.set_title(label)
        ax.grid(True, alpha=0.3)
        ax.plot(xx_sim, arr, color="C3", lw=1.4)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=120)
    plt.close(fig)


def render_movie(args, xx_sim, frames_sim, fps=15):
    prims = [primitives(fr[2]) for fr in frames_sim]
    t_arr = np.array([fr[0] for fr in frames_sim])
    ni = np.stack([p[0] for p in prims])
    nn = np.stack([p[1] for p in prims])
    v  = np.stack([p[2] for p in prims])
    Ti = np.stack([p[6] for p in prims])
    f  = ni / (ni + nn)

    z_m   = np.linspace(0.0, args.top_height_km * 1000.0, 200)
    z_km  = z_m / 1.0e3
    Bz_G  = B_canopy(z_m) * 1.0e4
    A_rel = B0_T / B_canopy(z_m)

    def ylim(arr, log=False):
        if log:
            a = arr[arr > 0]
            return float(a.min()) / 1.5, float(a.max()) * 1.5
        lo, hi = float(arr.min()), float(arr.max())
        pad = 0.06 * (hi - lo if hi > lo else max(abs(hi), 1.0))
        return lo - pad, hi + pad

    ylims_bot = [ylim(f), ylim(Ti), ylim(v), ylim(nn, log=True)]

    nT = len(t_arr)
    save_frames_parallel(
        _render_canopy_frame,
        [(k, xx_sim, f[k], Ti[k], v[k], nn[k],
          f[0], Ti[0], Bz_G, z_km, A_rel, t_arr[k], nT, ylims_bot)
         for k in range(nT)],
        args.out_mp4, fps,
    )
    print(f"[viz] wrote {args.out_mp4}  ({nT} frames @ {fps} fps)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim-output", default="out_canopy.txt")
    ap.add_argument("--out-png", default="visualization/analytic_canopy/canopy_geometry_evolution.png")
    ap.add_argument("--out-mp4", default="visualization/analytic_canopy/canopy_evolution.mp4")
    ap.add_argument("--top-height-km", type=float, default=986.0,
                    help="vertical extent of the flux-tube render (Mm)")
    ap.add_argument("--fps", type=int, default=15)
    args = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out_png)), exist_ok=True)

    print(f"[viz] reading {args.sim_output}")
    xx, frames = read_frames(args.sim_output)
    print(f"  ns={xx.size}, frames={len(frames)}, t in [{frames[0][0]:.1f}, {frames[-1][0]:.1f}] s")

    print(f"[viz] rendering combined PNG")
    render_combined(args, xx, frames)
    print(f"[viz] rendering movie")
    render_movie(args, xx, frames, fps=args.fps)
    print(f"[viz] done.")


if __name__ == "__main__":
    main()
