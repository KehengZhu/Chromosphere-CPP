#!/usr/bin/env python3
"""Render the SIMULATED event ensemble on a 3-D solar sphere.

Each field line is drawn along its real PFSS 3-D path, coloured by the temperature
from its own chromo_main run — the 1-D solution T(s,t) is mapped onto the 3-D path
by arc length (closed loops: the full loop; open lines: the simulated lower column,
with the un-simulated corona above the cap drawn faded). This is the campaign's
1-D→3-D step (Phase P5, single-frame + time animation), analogous to Noraz et al.
(2026) Fig. 1 but with our flare temperatures instead of a static field.

Inputs: the ensemble index from extract_event_ensemble.py (lines + 3-D paths + the
.dat/.txt of each run) and the ADAPT magnetogram (for the photosphere colour).

Outputs (util/visualization/):
  event_ensemble_3d_peak.png   — snapshot at --snap-time (flare-filled loops)
  event_ensemble_3d.mp4        — time evolution: the flare lighting up the lines

Usage:
  python util/visualize_event_3d.py \
    --lines-json outputs/event_ensemble_lines.json \
    --magnetogram util/data/adapt40311_044012_202408010600_i00053600n1.fts.gz
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from extract_event_loop import load_magnetogram, R_SUN_M           # noqa: E402
from animate_loop_flare import parse_dat, read_frames, primitives  # noqa: E402

R_SUN_MM = R_SUN_M / 1.0e6
T_LO, T_HI = 8.0e3, 1.2e7        # K, LogNorm range (chromosphere → flare corona)
RSS = 2.5                         # source-surface radius [R_sun]
FADE = (0.62, 0.62, 0.66, 0.20)   # colour for un-simulated open-line corona
# Field lines are tens of Mm tall = a few % of R_sun, invisible at true scale on the
# sphere. Exaggerate height nonlinearly for VISUALISATION ONLY: r → 1 + G·(r−1)^P.
# The sqrt (P<1) lifts the short closed loops a lot while keeping the long open lines
# bounded. Arc-length (hence the T mapping) uses the TRUE geometry; only the plotted
# radius is stretched. Stated in the caption.
RAD_GAIN, RAD_P = 1.3, 0.5


def lonlat_to_xyz(lon_deg, lat_deg, r=1.0):
    lo, la = np.deg2rad(lon_deg), np.deg2rad(lat_deg)
    cl = np.cos(la)
    return np.array([r * cl * np.cos(lo), r * cl * np.sin(lo), r * np.sin(la)])


def exaggerate(xyz):
    """Stretch radius for visibility (visualisation only)."""
    r = np.linalg.norm(xyz, axis=-1, keepdims=True)
    rn = 1.0 + RAD_GAIN * np.clip(r - 1.0, 0.0, None) ** RAD_P
    return xyz / np.clip(r, 1e-9, None) * rn


def resample_smooth(xyz, n=120):
    """Smooth, evenly-resample a 3-D polyline (coarse PFSS traces → clean curves)."""
    if xyz.shape[0] < 3:
        return xyz
    try:
        from scipy.interpolate import splprep, splev
        k = min(3, xyz.shape[0] - 1)
        (tck, _u), _ = splprep(xyz.T, s=0, k=k), None
        uu = np.linspace(0, 1, n)
        return np.array(splev(uu, tck)).T
    except Exception:
        # Linear fallback: resample by cumulative chord length.
        seg = np.linalg.norm(np.diff(xyz, axis=0), axis=1)
        u = np.concatenate([[0], np.cumsum(seg)]); u /= u[-1]
        uu = np.linspace(0, 1, n)
        return np.stack([np.interp(uu, u, xyz[:, j]) for j in range(3)], axis=1)


def build_surface(magnetogram):
    """Sphere (X,Y,Z) + Br facecolors from the ADAPT magnetogram (downsampled)."""
    import matplotlib
    smap = load_magnetogram(magnetogram)
    data = np.asarray(smap.data, float)
    ny, nx = data.shape
    sx, sy = max(1, nx // 360), max(1, ny // 180)
    Bg = data[::sy, ::sx]
    ny2, nx2 = Bg.shape
    lon = np.linspace(0.0, 2 * np.pi, nx2, endpoint=False) + np.pi / nx2
    sinlat = np.linspace(-1.0, 1.0, ny2, endpoint=False) + 1.0 / ny2
    lat = np.arcsin(sinlat)
    LON, LAT = np.meshgrid(lon, lat)
    cl = np.cos(LAT)
    Xs, Ys, Zs = cl * np.cos(LON), cl * np.sin(LON), np.sin(LAT)
    clip = max(20.0, float(np.nanpercentile(np.abs(Bg), 99.5)))
    norm = matplotlib.colors.Normalize(-clip, clip)
    fc = matplotlib.colormaps["RdBu_r"](norm(np.clip(Bg, -clip, clip)))
    return Xs, Ys, Zs, fc, (ny2, nx2)


def load_line_temperatures(lines, t_grid):
    """For each line read its run, map T(s) onto the 3-D path at every t in t_grid.

    Returns a list of dicts: xyz (N,3, R_sun), seg_arc_Mm (N-1,) midpoint arc,
    T_grid (n_t, N-1) temperature per segment per time, sim_cap_Mm."""
    out = []
    for L in lines:
        ds, s_arc_m, phi_c, h = parse_dat(L["dat"])
        s_sim_Mm = s_arc_m / 1.0e6
        frames = read_frames(L["out"])
        ft = np.array([f[0] for f in frames])
        xyz0 = np.array(L["path_xyz_rsun"], float)             # (N,3) in R_sun
        if xyz0.shape[0] > 2000:                               # finely-traced/odd lines
            sel = np.linspace(0, xyz0.shape[0] - 1, 800).round().astype(int)
            xyz0 = xyz0[sel]
        xyz_real = resample_smooth(xyz0, n=140)                # smooth, even sampling
        seg = np.linalg.norm(np.diff(xyz_real, axis=0), axis=1)
        arc = np.concatenate([[0.0], np.cumsum(seg)]) * R_SUN_MM   # TRUE path arc [Mm]
        mid_arc = 0.5 * (arc[:-1] + arc[1:])                       # per-segment
        cap = float(L["sim_length_Mm"])
        Tg = np.full((len(t_grid), len(mid_arc)), np.nan)
        for it, tg in enumerate(t_grid):
            k = int(np.argmin(np.abs(ft - tg)))                    # nearest frame
            _ni, _nn, _v, Ti = primitives(frames[k][1], phi_c)
            T_at = np.interp(mid_arc, s_sim_Mm, Ti,
                             left=Ti[0], right=np.nan)             # beyond grid → nan
            T_at[mid_arc > cap] = np.nan                           # beyond sim cap
            Tg[it] = T_at
        out.append({"name": L["name"], "kind": L["kind"],
                    "xyz": exaggerate(xyz_real),               # stretched for plotting
                    "T_grid": Tg, "cap": cap})
        print(f"  [{L['name']:12s}] {L['kind']:6s} peak T {np.nanmax(Tg)/1e6:.1f} MK")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lines-json", default="outputs/event_ensemble_lines.json")
    ap.add_argument("--magnetogram", required=True)
    ap.add_argument("--out-dir", default="util/visualization")
    ap.add_argument("--snap-time", type=float, default=40.0, help="static PNG time [s]")
    ap.add_argument("--t-max", type=float, default=160.0)
    ap.add_argument("--n-frames", type=int, default=90)
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--zoom-rsun", type=float, default=0.55,
                    help="half-width [R_sun] of the AR zoom box (loops are tens of Mm "
                         "= small fractions of R_sun, invisible on the full sphere)")
    ap.add_argument("--no-anim", action="store_true")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    idx = json.load(open(args.lines_json))
    lines = idx["lines"]
    ar_lon, ar_lat = idx["ar_carrington_lon_deg"], idx["ar_carrington_lat_deg"]
    os.makedirs(args.out_dir, exist_ok=True)

    t_grid = np.linspace(0.0, args.t_max, args.n_frames)
    print(f"[viz3d] reading {len(lines)} runs, mapping T(s,t) onto 3-D paths …")
    L = load_line_temperatures(lines, t_grid)

    print("[viz3d] building sphere + scene …")
    Xs, Ys, Zs, fc, (ny2, nx2) = build_surface(args.magnetogram)
    cmap = matplotlib.colormaps["inferno"].copy()
    norm = matplotlib.colors.LogNorm(T_LO, T_HI)

    fig = plt.figure(figsize=(10, 9))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_box_aspect((1, 1, 1))
    ax.plot_surface(Xs, Ys, Zs, rcount=ny2 - 1, ccount=nx2 - 1,
                    facecolors=fc, shade=False, linewidth=0, antialiased=False, zorder=1)
    rss_plot = 1.0 + RAD_GAIN * (RSS - 1.0) ** RAD_P          # exaggerated source surf
    th = np.linspace(0, 2 * np.pi, 200)
    ca, sa = np.cos(th), np.sin(th)
    ax.plot(rss_plot * ca, rss_plot * sa, 0 * ca, color="0.6", alpha=0.3, lw=0.6, zorder=2)
    ax.plot(rss_plot * ca, 0 * ca, rss_plot * sa, color="0.6", alpha=0.3, lw=0.6, zorder=2)
    # AR marker
    p = lonlat_to_xyz(ar_lon, ar_lat, 1.02)
    ax.scatter([p[0]], [p[1]], [p[2]], color="cyan", s=60, marker="*",
               edgecolor="k", linewidth=0.5, zorder=10, label="AR 13768")

    # One Line3DCollection per line; recoloured each frame.
    collections = []
    for ln in L:
        xyz = ln["xyz"]
        segs = np.stack([xyz[:-1], xyz[1:]], axis=1)              # (N-1,2,3)
        lc = Line3DCollection(segs, zorder=6)
        ax.add_collection3d(lc)
        collections.append(lc)

    def seg_colors(ln, it):
        T = ln["T_grid"][it]
        cols = cmap(norm(T))
        cols[np.isnan(T)] = FADE
        return cols

    def seg_widths(ln, it):
        T = ln["T_grid"][it]
        return np.where(np.isnan(T), 0.8, 2.9)

    def draw(it):
        for lc, ln in zip(collections, L):
            lc.set_color(seg_colors(ln, it))
            lc.set_linewidth(seg_widths(ln, it))
        ax.set_title(f"2024-08-01 M8.2 / AR 13768 — simulated ensemble on PFSS field\n"
                     f"t = {t_grid[it]:6.1f} s   (colour = electron temperature)",
                     fontsize=12)

    ax.set_axis_off()

    def set_view():
        R = 1.75                                           # tight frame: AR loops fill view
        ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_zlim(-R, R)
        ax.view_init(elev=14, azim=ar_lon + 45)            # 3/4 view, loops profile off the limb

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("electron temperature T [K]")
    fig.text(0.5, 0.05, f"{len(L)} field lines: "
             f"{sum(1 for x in L if x['kind']=='closed')} closed (full loop) + "
             f"{sum(1 for x in L if x['kind']=='open')} open (lower column simulated; "
             f"corona above faded). Photosphere colour = B$_r$. "
             f"Loop heights exaggerated ×(r−1)$^{{0.5}}$ for visibility.",
             ha="center", fontsize=9)

    it_snap = int(np.argmin(np.abs(t_grid - args.snap_time)))
    draw(it_snap)
    set_view()
    png = os.path.join(args.out_dir, "event_ensemble_3d_peak.png")
    fig.savefig(png, dpi=145, bbox_inches="tight")
    print(f"[viz3d] wrote {png}  (t={t_grid[it_snap]:.0f}s)")

    # Time-evolution movie: the flare lighting up the loops.
    if not args.no_anim:
        set_view()
        mp4 = os.path.join(args.out_dir, "event_ensemble_3d.mp4")
        writer = animation.FFMpegWriter(fps=args.fps, codec="libx264",
                                        extra_args=["-pix_fmt", "yuv420p"])
        anim = animation.FuncAnimation(fig, draw, frames=len(t_grid),
                                       interval=1000.0 / args.fps, blit=False)
        anim.save(mp4, writer=writer, dpi=120)
        print(f"[viz3d] wrote {mp4}  ({len(t_grid)} frames @ {args.fps} fps = "
              f"{len(t_grid)/args.fps:.1f} s)")
    plt.close(fig)


if __name__ == "__main__":
    main()
