#!/usr/bin/env python3
"""Visualize 2-3 PFSS field lines side by side: where they are on the Sun, what
boundary state they're pinned to, and how the chromospheric simulation evolves
along each.

Layout of the output PNG:
  Row 1 (full width):
      Left  - 3D PFSS sphere with the chosen field lines highlighted
      Right - boundary-condition summary table (footpoint, |B|, topology, n_e/n_n/T at outer ghost)
  Rows 2-4 (one per line):
      4 panels per row showing time evolution along height of:
          ionization fraction f, T_i, V (ion velocity), n_n  (or P_n, n_i...)

And a companion MP4 (per-line strip animation) is also written.

Usage (from repo root, .venv active):
    python util/visualize_pfss_lines_evolution.py
"""

from __future__ import annotations

import argparse
import os
import re
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
from extract_field_line import _import_pfsspy, load_magnetogram
from visualize_pfss_3d import make_tracer, trace_seeds, R_SUN_M

import astropy.units as u
from astropy.coordinates import SkyCoord
import sunpy.coordinates  # noqa: F401


HEADER_RX = re.compile(
    r"# Footpoint:\s*lon=([\-0-9.]+)\s*deg,\s*lat=([\-0-9.]+)\s*deg,\s*\|B_LOS\|=([\-0-9.]+)\s*G"
)


def parse_dat_header(path):
    fp_lon = fp_lat = fp_B = None
    topology = "unknown"
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                # also pull ghost rows.
                continue
            m = HEADER_RX.search(line)
            if m:
                fp_lon, fp_lat, fp_B = float(m.group(1)), float(m.group(2)), float(m.group(3))
    # extract ghost row from [GHOSTS] section.
    outer_ne = outer_nn = outer_T = outer_Tfac = None
    inner_ne = inner_nn = inner_T = None
    in_ghosts = False
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if "[GHOSTS]" in line:
                in_ghosts = True
                continue
            if not in_ghosts or not line or line.startswith("#"):
                continue
            toks = line.split()
            if len(toks) < 7:
                continue
            tag = toks[0]
            B   = float(toks[1])
            phi = float(toks[2])
            ne  = float(toks[3])
            nn  = float(toks[4])
            T   = float(toks[5])
            tfac= float(toks[6])
            if tag == "outer_0":
                outer_ne, outer_nn, outer_T, outer_Tfac = ne, nn, T, tfac
            elif tag == "inner_0":
                inner_ne, inner_nn, inner_T = ne, nn, T
    return dict(
        fp_lon=fp_lon, fp_lat=fp_lat, fp_B_LOS=fp_B,
        outer_ne=outer_ne, outer_nn=outer_nn, outer_T=outer_T, outer_Tfac=outer_Tfac,
        inner_ne=inner_ne, inner_nn=inner_nn, inner_T=inner_T,
    )


def trace_full_lines(output, footpoints, tracing):
    """Trace each footpoint upward through the entire PFSS volume and return
    (xyz_Rsun_list, is_open_list, n_pts_list)."""
    lons = np.array([fp["fp_lon"] for fp in footpoints])
    lats = np.array([fp["fp_lat"] for fp in footpoints])
    r_seed_m = R_SUN_M * 1.005
    lines, is_open = trace_seeds(output, lons, lats, r_seed_m, tracing)
    return lines, is_open


def render_combined(out_png, dats, frames_per_line, fp_meta, full_traces, openness, fig_w=20):
    n_lines = len(dats)
    colors = ["C0", "C2", "C3"]

    # Time grid common to all three (driven by line 0).
    xx, frames0 = frames_per_line[0]
    t_arr = np.array([fr[0] for fr in frames0])
    step_arr = np.array([fr[1] for fr in frames0])

    # Layout: left column holds 3D plot (rows 0..n_lines-1) + BC table (row n_lines).
    # Right four columns hold the per-line evolution panels (one row per line).
    fig = plt.figure(figsize=(fig_w, 3.4 * n_lines + 2.0))
    gs = gridspec.GridSpec(
        n_lines + 1, 5,
        height_ratios=[1.0] * n_lines + [0.55],
        width_ratios=[2.2, 1.0, 1.0, 1.0, 1.0],
        hspace=0.45, wspace=0.40,
        left=0.03, right=0.965, top=0.94, bottom=0.06,
    )

    # ---- left column, top n_lines rows: 3D sphere with the traced lines ------
    ax3d = fig.add_subplot(gs[:n_lines, 0], projection="3d")
    ax3d.set_box_aspect((1, 1, 1))
    # tiny sphere stub for context
    u_ang = np.linspace(0, 2 * np.pi, 80)
    v_ang = np.linspace(-np.pi / 2, np.pi / 2, 40)
    U, V = np.meshgrid(u_ang, v_ang)
    Xs = np.cos(V) * np.cos(U)
    Ys = np.cos(V) * np.sin(U)
    Zs = np.sin(V)
    ax3d.plot_surface(Xs, Ys, Zs, color="0.85", alpha=0.35, linewidth=0,
                      shade=False, rcount=40, ccount=80, zorder=1)

    # source-surface reference rings
    rss = 2.5
    th = np.linspace(0, 2 * np.pi, 200)
    ax3d.plot(rss * np.cos(th), rss * np.sin(th), np.zeros_like(th),
              color="0.6", alpha=0.3, lw=0.5)
    ax3d.plot(rss * np.cos(th), np.zeros_like(th), rss * np.sin(th),
              color="0.6", alpha=0.3, lw=0.5)

    # footpoints + full traces
    for i, (fp, xyz, op) in enumerate(zip(fp_meta, full_traces, openness)):
        c = colors[i]
        # footpoint dot
        lon = np.deg2rad(fp["fp_lon"])
        lat = np.deg2rad(fp["fp_lat"])
        xp = np.cos(lat) * np.cos(lon)
        yp = np.cos(lat) * np.sin(lon)
        zp = np.sin(lat)
        ax3d.scatter([xp], [yp], [zp], color=c, s=60, edgecolor="black",
                     linewidth=0.6, zorder=10)
        # trace through volume
        ax3d.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], color=c, lw=1.6, zorder=5,
                  label=f"line {chr(ord('A')+i)}  ({'open' if op else 'closed'})")
    # Zoom around the loops (closed loops don't reach the source surface, so
    # crop the view to ~1.4 R_sun for visibility). Centre on the cluster of
    # footpoints we picked.
    fp_xyz = np.array([
        [np.cos(np.deg2rad(fp["fp_lat"])) * np.cos(np.deg2rad(fp["fp_lon"])),
         np.cos(np.deg2rad(fp["fp_lat"])) * np.sin(np.deg2rad(fp["fp_lon"])),
         np.sin(np.deg2rad(fp["fp_lat"]))]
        for fp in fp_meta
    ])
    cx, cy, cz = fp_xyz.mean(axis=0)
    R = 1.0
    ax3d.set_xlim(cx - R, cx + R)
    ax3d.set_ylim(cy - R, cy + R)
    ax3d.set_zlim(cz - R, cz + R)
    ax3d.set_axis_off()
    ax3d.legend(loc="upper left", fontsize=9, framealpha=0.85)
    ax3d.set_title("PFSS field lines traced from footpoints\n(grey sphere = Sun, faint rings = source surface r = 2.5 R⊙)", fontsize=11)
    ax3d.view_init(elev=18, azim=-50)

    # ---- bottom-left: BC summary table (spans first 2 columns) ---------------
    ax_table = fig.add_subplot(gs[n_lines, :])
    ax_table.axis("off")
    rows = []
    rows.append(["line", "lon(°)", "lat(°)", "|B_LOS|(G)", "topo",
                 "outer n_e (m⁻³)", "outer n_n (m⁻³)", "outer T_n (K)",
                 "outer T_e (K)", "f_outer"])
    for i, (fp, op) in enumerate(zip(fp_meta, openness)):
        f_outer = fp["outer_ne"] / (fp["outer_ne"] + fp["outer_nn"])
        T_e = fp["outer_T"] * fp["outer_Tfac"]
        T_n = fp["outer_T"]
        rows.append([
            chr(ord("A") + i),
            f"{fp['fp_lon']:.2f}",
            f"{fp['fp_lat']:.2f}",
            f"{fp['fp_B_LOS']:+.2f}",
            ("open" if op else "closed"),
            f"{fp['outer_ne']:.2e}",
            f"{fp['outer_nn']:.2e}",
            f"{T_n:.0f}",
            f"{T_e:.0f}",
            f"{f_outer:.3f}",
        ])
    tbl = ax_table.table(cellText=rows, loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10.0)
    tbl.scale(1, 1.9)
    # colour header row + line-label cells
    for j in range(len(rows[0])):
        tbl[(0, j)].set_facecolor("#dddddd")
        tbl[(0, j)].set_text_props(weight="bold")
    for i in range(n_lines):
        tbl[(i + 1, 0)].set_facecolor(matplotlib.colors.to_hex(colors[i]))
        tbl[(i + 1, 0)].set_text_props(color="white", weight="bold")
    ax_table.set_title(
        "Outer-boundary state (Dirichlet ghost cell, T_e doubled per file)\n"
        "n_e, n_n, T_n from C7 spline extrapolation at the line top; T_e = T_n × T_e_factor",
        fontsize=10)

    # ---- per-line evolution rows ---------------------------------------------
    field_specs = [
        ("f",   "ionization fraction f",      False,  None),
        ("Ti",  r"$T_i$ (K)",                False,  None),
        ("V",   r"$V$ (m/s, ion velocity)",  False,  None),
        ("nn",  r"$n_n$ (m$^{-3}$)",          True,   None),
    ]
    field_labels = {k: lbl for k, lbl, _, _ in field_specs}
    # height axis is shared but each line has its own xx (same length, similar range).
    # build cell-centered height axis: dat 1 km grid centered on midpoints. We
    # already have xx from frames file (cumulative km, offset).
    for i, ((xx, frames), fp) in enumerate(zip(frames_per_line, fp_meta)):
        # compute all primitives once per frame, then pick chosen fields
        prims = [primitives(fr[2]) for fr in frames]
        ni_all = np.stack([p[0] for p in prims], axis=0)
        nn_all = np.stack([p[1] for p in prims], axis=0)
        v_all  = np.stack([p[2] for p in prims], axis=0)
        Ti_all = np.stack([p[6] for p in prims], axis=0)
        f_all  = ni_all / (ni_all + nn_all)
        t_arr  = np.array([fr[0] for fr in frames])

        # use a sequential colormap to display time evolution as lines
        # (use ~12 traces from early to late)
        n_show = 12
        idxs = np.linspace(0, len(prims) - 1, n_show).astype(int)
        cmap = colormaps["viridis"]

        panels = [("f", f_all), ("Ti", Ti_all), ("V", v_all), ("nn", nn_all)]
        for j, (key, arr) in enumerate(panels):
            # right four columns of row i
            ax = fig.add_subplot(gs[i, j + 1])
            for k, idx in enumerate(idxs):
                t = t_arr[idx]
                color = cmap(k / max(1, n_show - 1))
                ax.plot(xx, arr[idx], color=color, lw=0.9, alpha=0.95)
            # overlay first + last bold in line's signature color
            ax.plot(xx, arr[0],  color="black", lw=1.4, ls=":",  alpha=0.9, label=f"t=0 s")
            ax.plot(xx, arr[-1], color=colors[i], lw=1.6, alpha=0.95, label=f"t={t_arr[-1]:.0f} s")
            ax.grid(True, alpha=0.3)
            label = field_labels[key]
            if j == 0:
                ax.set_ylabel(f"line {chr(ord('A')+i)}", fontsize=11,
                              color=colors[i], weight="bold")
            if i == n_lines - 1:
                ax.set_xlabel("height along s (km)")
            # column header on top row
            if i == 0:
                ax.set_title(label, fontsize=11)
            if key == "nn":
                ax.set_yscale("log")
            if key == "V":
                ax.axhline(0, color="0.7", lw=0.6)
            if i == 0 and j == 0:
                ax.legend(loc="best", fontsize=8)
        # add a colourbar for time on the right of the last panel
        if i == 0:
            sm = matplotlib.cm.ScalarMappable(
                cmap=cmap, norm=matplotlib.colors.Normalize(vmin=0, vmax=t_arr[-1]))
            cbar_ax = fig.add_axes([0.978, 0.30, 0.008, 0.50])
            cbar = fig.colorbar(sm, cax=cbar_ax)
            cbar.set_label("simulation time (s)", fontsize=10)

    fig.suptitle(
        "PFSS field-line scenarios — geometry, BC, and chromospheric relaxation",
        fontsize=13, y=0.985,
    )
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    print(f"[viz] wrote {out_png}")
    plt.close(fig)


def _render_pfss_lines_frame(k, tmpdir, t_k, nF, xx_all, sliced_arrays,
                              ylims, fp_meta, colors, field_specs):
    """Render one frame; sliced_arrays[p] is a 1-D array for panel p."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_lines = len(xx_all)
    fig, axes = plt.subplots(4, n_lines, figsize=(4.6 * n_lines, 10.5))
    if n_lines == 1:
        axes = [[axes[r]] for r in range(4)]
    fig.suptitle(
        f"PFSS field-line chromospheric evolution  |  t = {t_k:6.2f} s   ({k+1}/{nF})",
        fontsize=13,
    )

    for j, (key, label, logy) in enumerate(field_specs):
        for i in range(n_lines):
            panel_idx = j * n_lines + i
            ax = axes[j][i]
            ax.set_xlim(xx_all[i].min(), xx_all[i].max())
            ax.set_ylim(*ylims[panel_idx])
            if logy:
                ax.set_yscale("log")
            ax.grid(True, alpha=0.3)
            if j == 0:
                ax.set_title(
                    f"line {chr(ord('A')+i)}  |  "
                    f"lon={fp_meta[i]['fp_lon']:.1f}°, lat={fp_meta[i]['fp_lat']:.1f}°  "
                    f"|B_LOS|={fp_meta[i]['fp_B_LOS']:+.1f} G",
                    fontsize=10, color=colors[i],
                )
            if i == 0:
                ax.set_ylabel(label)
            if j == 3:
                ax.set_xlabel("height s (km)")
            ax.plot(xx_all[i], sliced_arrays[panel_idx], color=colors[i], lw=1.5)
            if key == "V":
                ax.axhline(0, color="0.6", lw=0.5)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=120)
    plt.close(fig)


def render_movie(out_mp4, dats, frames_per_line, fp_meta, fps=15):
    """Animation: n_lines columns × 4 rows of fields, one column per line."""
    n_lines = len(dats)
    colors  = ["C0", "C2", "C3"]
    field_specs = [
        ("f",  "ionization fraction f",  False),
        ("Ti", r"$T_i$ (K)",             False),
        ("V",  r"$V$ (m/s)",             False),
        ("nn", r"$n_n$ (m$^{-3}$)",      True),
    ]

    xx0, frames0 = frames_per_line[0]
    t0 = np.array([fr[0] for fr in frames0])
    nF = len(t0)

    # Pre-compute stacked field arrays and time axes per line.
    stacked = []   # stacked[i] = {key: (nT, ns) array} for line i
    t_all  = []
    xx_all = []
    for xx, frames in frames_per_line:
        prims = [primitives(fr[2]) for fr in frames]
        ni = np.stack([p[0] for p in prims])
        nn = np.stack([p[1] for p in prims])
        v  = np.stack([p[2] for p in prims])
        Ti = np.stack([p[6] for p in prims])
        stacked.append({"f": ni / (ni + nn), "Ti": Ti, "V": v, "nn": nn})
        t_all.append(np.array([fr[0] for fr in frames]))
        xx_all.append(xx)

    # Build flat ylims list (row-major: row 0 all lines, row 1 all lines, ...)
    ylims = []
    for key, _, logy in field_specs:
        for i in range(n_lines):
            arr = stacked[i][key]
            if logy:
                a = arr[arr > 0]
                ylims.append((float(a.min()) / 1.5, float(a.max()) * 1.5))
            else:
                lo, hi = float(arr.min()), float(arr.max())
                pad = 0.06 * (hi - lo if hi > lo else max(abs(hi), 1.0))
                ylims.append((lo - pad, hi + pad))

    # Build per-frame args: pre-slice to 1-D arrays so IPC sends minimal data.
    frame_args = []
    for k in range(nF):
        tn = t0[k]
        idx = [int(np.argmin(np.abs(t_all[i] - tn))) for i in range(n_lines)]
        sliced = [stacked[p % n_lines][fld][idx[p % n_lines]]
                  for p, (fld, _, _) in enumerate(
                      f for f in field_specs for _ in range(n_lines))]
        frame_args.append((k, t0[k], nF, xx_all, sliced, ylims, fp_meta, colors, field_specs))

    save_frames_parallel(_render_pfss_lines_frame, frame_args, out_mp4, fps)
    print(f"[viz] wrote {out_mp4}  ({nF} frames @ {fps} fps)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--magnetogram", default="util/data/mrzqs190801t0014c2220_229.fits")
    ap.add_argument("--dats", nargs="+", default=[
        "scenarios/data/pfss_lineA.dat",
        "scenarios/data/pfss_lineB.dat",
        "scenarios/data/pfss_lineC.dat",
    ])
    ap.add_argument("--outs", nargs="+", default=[
        "out_lineA.txt", "out_lineB.txt", "out_lineC.txt",
    ])
    ap.add_argument("--out-png", default="visualization/pfss/pfss_lines_evolution.png")
    ap.add_argument("--out-mp4", default="visualization/pfss/pfss_lines_evolution.mp4")
    ap.add_argument("--fps", type=int, default=15)
    ap.add_argument("--nrho", type=int, default=60)
    ap.add_argument("--rss",  type=float, default=2.5)
    args = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out_png)), exist_ok=True)

    assert len(args.dats) == len(args.outs), "dats/outs length mismatch"
    print(f"[viz] loading {len(args.dats)} PFSS data files")
    fp_meta = [parse_dat_header(p) for p in args.dats]
    for i, fp in enumerate(fp_meta):
        print(f"  line {chr(ord('A')+i)}: lon={fp['fp_lon']:.2f}, lat={fp['fp_lat']:.2f}, "
              f"|B_LOS|={fp['fp_B_LOS']:+.2f} G; outer (n_e, n_n, T) = "
              f"({fp['outer_ne']:.2e}, {fp['outer_nn']:.2e}, {fp['outer_T']:.0f}); "
              f"T_e_factor={fp['outer_Tfac']}")

    print(f"[viz] solving PFSS for the geometry view")
    pfss, tracing, pkg = _import_pfsspy()
    smap = load_magnetogram(args.magnetogram)
    pfss_in = pfss.Input(smap, args.nrho, args.rss)
    output  = pfss.pfss(pfss_in)
    full_traces, openness = trace_full_lines(output, fp_meta, tracing)
    print(f"[viz] traced {len(full_traces)} full PFSS lines from the footpoints "
          f"({sum(openness)} open, {len(full_traces) - sum(openness)} closed)")

    print(f"[viz] reading simulation outputs")
    frames_per_line = []
    for p in args.outs:
        xx, frames = read_frames(p)
        frames_per_line.append((xx, frames))
        print(f"  {p}: ns={xx.size}, frames={len(frames)}, t in [{frames[0][0]:.1f}, {frames[-1][0]:.1f}] s")

    print(f"[viz] rendering combined PNG")
    render_combined(args.out_png, args.dats, frames_per_line, fp_meta, full_traces, openness)
    print(f"[viz] rendering evolution movie")
    render_movie(args.out_mp4, args.dats, frames_per_line, fp_meta, fps=args.fps)
    print("[viz] done.")


if __name__ == "__main__":
    main()
