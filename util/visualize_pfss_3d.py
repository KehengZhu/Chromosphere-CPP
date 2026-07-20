#!/usr/bin/env python3
"""3D visualisation of a PFSS extrapolation of a GONG synoptic magnetogram,
analogous to Fig. 1 of Noraz et al. 2026 (A&A 705, A86).

Differences from Noraz Fig. 1:
  - Noraz colours field lines by temperature from a Bifrost rMHD run; PFSS gives
    only a potential field, so we colour field lines by altitude (h = r - R_sun)
    instead. Pink/red ≈ low heights (chromosphere), green ≈ source surface.
  - Noraz overlays vertical v_z on a side slice; PFSS has no velocity field, so
    we draw the source-surface (r=2.5 R_sun) shell as a transparent reference.
  - Noraz boxes 12 × 12 × 10.5 Mm of one QS region; the "local" rendering here
    boxes a similar patch around the same footpoint used by extract_field_line.py.

Outputs (under util/):
  - pfss_3d_global.png  : full-sphere view of the photosphere + PFSS lines
  - pfss_3d_local.png   : Noraz-style box around the chosen footpoint
  - pfss_3d_rotate.mp4  : rotating-camera animation of the global view

Usage (from repo root, with .venv active):
    python util/visualize_pfss_3d.py \\
        --magnetogram util/data/mrzqs190801t0014c2220_229.fits \\
        --nseeds 96 --nrho 60 --rss 2.5

Tunable knobs (--help for the full list).
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings

import numpy as np

R_SUN_M = 6.957e8

# Re-use helpers from extract_field_line.py.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from extract_field_line import _import_pfsspy, load_magnetogram, pick_network_footpoint


def lonlat_grid_seeds(n, lat_window=(-60.0, 60.0)):
    """Return n approximately-uniform (lon, lat) seed pairs inside the window.

    Uses a Fibonacci lattice on the lat-restricted band so seeds aren't clumped
    at the poles. Returns degrees.
    """
    lo, hi = lat_window
    # Map [0,1] -> [lo, hi] in sin(lat) so seeds are uniform per solid angle.
    sin_lo, sin_hi = np.sin(np.deg2rad(lo)), np.sin(np.deg2rad(hi))
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    ks = np.arange(n)
    sin_lat = sin_lo + (sin_hi - sin_lo) * (ks + 0.5) / n
    lats = np.rad2deg(np.arcsin(sin_lat))
    lons = (ks * 360.0 / phi) % 360.0
    return lons, lats


def make_tracer(tracing):
    if hasattr(tracing, "PerformanceTracer"):
        return tracing.PerformanceTracer()
    if hasattr(tracing, "FortranTracer"):
        return tracing.FortranTracer()
    return tracing.PythonTracer()


def trace_seeds(output, lons_deg, lats_deg, r_seed_m, tracing):
    """Trace PFSS field lines from (lon, lat) seeds at radius r_seed_m.

    Returns a list of (N_k, 3) Cartesian arrays in units of R_sun, each line as
    (x, y, z) along the trace.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    import sunpy.coordinates  # noqa: F401

    obstime = output.input_map.date if hasattr(output, "input_map") else None
    seeds = SkyCoord(
        lons_deg * u.deg,
        lats_deg * u.deg,
        np.full_like(lons_deg, r_seed_m) * u.m,
        frame="heliographic_carrington",
        obstime=obstime,
        observer="earth",
    )
    tracer = make_tracer(tracing)
    fls = tracer.trace(seeds, output)
    lines = []
    is_open = []
    for fl in fls:
        try:
            xyz = np.stack(
                [
                    fl.coords.cartesian.x.to_value(u.m) / R_SUN_M,
                    fl.coords.cartesian.y.to_value(u.m) / R_SUN_M,
                    fl.coords.cartesian.z.to_value(u.m) / R_SUN_M,
                ],
                axis=-1,
            )
            if xyz.shape[0] < 2:
                continue
        except Exception:
            continue
        lines.append(xyz)
        is_open.append(bool(getattr(fl, "is_open", False)))
    return lines, is_open


def lonlat_to_xyz(lon_deg, lat_deg, r_over_rsun=1.0):
    lon = np.deg2rad(lon_deg)
    lat = np.deg2rad(lat_deg)
    cl = np.cos(lat)
    return np.stack(
        [r_over_rsun * cl * np.cos(lon),
         r_over_rsun * cl * np.sin(lon),
         r_over_rsun * np.sin(lat)],
        axis=-1,
    )


def render_global(smap, lines, is_open, out_path, dpi=140, figsize=(9, 9),
                   rotate_mp4=None, fps=24, frames=72):
    """Render a 3D sphere of the photosphere with PFSS field lines on top."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    import matplotlib.animation as animation

    # Build a sphere surface coloured by B_LOS (G). Downsample the magnetogram
    # to a manageable resolution.
    data = np.asarray(smap.data, dtype=float)
    ny, nx = data.shape
    stride_x = max(1, nx // 512)
    stride_y = max(1, ny // 256)
    Bg = data[::stride_y, ::stride_x]
    ny2, nx2 = Bg.shape
    lon = np.linspace(0.0, 2.0 * np.pi, nx2, endpoint=False) + (np.pi / nx2)
    sin_lat = np.linspace(-1.0, 1.0, ny2, endpoint=False) + (1.0 / ny2)
    lat = np.arcsin(sin_lat)
    LON, LAT = np.meshgrid(lon, lat)
    cl = np.cos(LAT)
    Xs, Ys, Zs = cl * np.cos(LON), cl * np.sin(LON), np.sin(LAT)

    # Quiet-Sun B_LOS is dominated by ≤30 G structure; clip tight for contrast.
    Bclip = np.clip(Bg, -30.0, 30.0)
    norm = matplotlib.colors.Normalize(vmin=-30.0, vmax=30.0)
    cmap = matplotlib.colormaps["RdBu_r"]
    facecolors = cmap(norm(Bclip))

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_box_aspect((1, 1, 1))
    # plot_surface ignores facecolors past rcount × ccount, so force full
    # resolution by passing rcount = ny2-1, ccount = nx2-1 (cells, not vertices).
    ax.plot_surface(
        Xs, Ys, Zs,
        rcount=ny2 - 1, ccount=nx2 - 1,
        facecolors=facecolors, shade=False, linewidth=0, antialiased=False,
        zorder=1,
    )

    # Source-surface as a faint equatorial+meridional reference ring, not a full
    # wireframe (a full wireframe obscures the photospheric map).
    rss = 2.5
    th = np.linspace(0.0, 2.0 * np.pi, 200)
    ax.plot(rss * np.cos(th), rss * np.sin(th), np.zeros_like(th),
            color="0.55", alpha=0.35, lw=0.6, zorder=4)
    ax.plot(rss * np.cos(th), np.zeros_like(th), rss * np.sin(th),
            color="0.55", alpha=0.35, lw=0.6, zorder=4)

    # Field lines coloured by mean altitude (h = r - R_sun) along the trace, in Mm.
    cmap_h = matplotlib.colormaps["plasma"]
    h_lo, h_hi = 0.0, (rss - 1.0) * R_SUN_M / 1.0e6  # Mm
    norm_h = matplotlib.colors.Normalize(vmin=h_lo, vmax=h_hi)

    n_open = 0
    for xyz, openness in zip(lines, is_open):
        rs = np.linalg.norm(xyz, axis=-1)
        h_Mm = (rs - 1.0) * R_SUN_M / 1.0e6
        h_mean = float(np.mean(h_Mm))
        color = cmap_h(norm_h(h_mean))
        ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2],
                color=color, lw=0.9, alpha=0.95, zorder=5)
        if openness:
            n_open += 1

    R = rss * 1.05
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    ax.set_zlim(-R, R)
    ax.set_axis_off()
    ax.view_init(elev=22, azim=210)  # rotate to bring the active band into view

    title = (f"PFSS extrapolation of GONG synoptic — {len(lines)} field lines\n"
             f"({n_open} open, {len(lines) - n_open} closed; source surface at "
             f"r={rss} R⊙; photosphere colour = B$_{{LOS}}$, ±30 G)")
    fig.suptitle(title, fontsize=11)

    # Colour bars: photospheric B and field-line altitude.
    cax1 = fig.add_axes([0.08, 0.07, 0.36, 0.018])
    matplotlib.colorbar.ColorbarBase(
        cax1, cmap=cmap, norm=norm, orientation="horizontal",
        label="photospheric B$_{LOS}$ (G, clipped at ±30)",
    )
    cax2 = fig.add_axes([0.56, 0.07, 0.36, 0.018])
    matplotlib.colorbar.ColorbarBase(
        cax2, cmap=cmap_h, norm=norm_h, orientation="horizontal",
        label="mean altitude along field line (Mm)",
    )

    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    print(f"[viz] wrote {out_path}")

    if rotate_mp4 is not None:
        ax.view_init(elev=20, azim=-60)

        def update(k):
            ax.view_init(elev=20, azim=-60 + 360.0 * k / frames)
            return []

        writer = animation.FFMpegWriter(fps=fps, codec="libx264",
                                        extra_args=["-pix_fmt", "yuv420p"])
        anim = animation.FuncAnimation(fig, update, frames=frames, interval=1000.0 / fps, blit=False)
        anim.save(rotate_mp4, writer=writer, dpi=dpi)
        print(f"[viz] wrote {rotate_mp4}  ({frames} frames @ {fps} fps = {frames/fps:.1f} s)")

    plt.close(fig)


def render_local(smap, output, fp_lon, fp_lat, tracing, out_path,
                 box_half_deg=4.0, top_r=1.005, n_seeds_side=12, dpi=140,
                 figsize=(9, 9)):
    """Noraz-Fig1-style box: photospheric B map at the base + traced field lines
    over a small lat/lon patch around the chosen footpoint.

    Coordinates: x = (lon - fp_lon) * R_sun_km, y = (lat - fp_lat) * R_sun_km in
    Mm; z = height above photosphere in Mm.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    import sunpy.coordinates  # noqa: F401

    deg2km = 2.0 * np.pi * R_SUN_M / 360.0 / 1000.0   # km per deg of arc on solar surface
    half_km = box_half_deg * deg2km
    half_Mm = half_km / 1000.0

    # Patch of the photospheric map inside ±box_half_deg around the footpoint.
    data = np.asarray(smap.data, dtype=float)
    ny, nx = data.shape
    lon_axis = np.linspace(0.0, 360.0, nx, endpoint=False) + (180.0 / nx)
    sin_lat = np.linspace(-1.0, 1.0, ny, endpoint=False) + (1.0 / ny)
    lat_axis = np.degrees(np.arcsin(sin_lat))

    lon_lo, lon_hi = fp_lon - box_half_deg, fp_lon + box_half_deg
    lat_lo, lat_hi = fp_lat - box_half_deg, fp_lat + box_half_deg
    ix_lo = int(np.searchsorted(lon_axis, lon_lo))
    ix_hi = int(np.searchsorted(lon_axis, lon_hi))
    iy_lo = int(np.searchsorted(lat_axis, lat_lo))
    iy_hi = int(np.searchsorted(lat_axis, lat_hi))
    ix_lo, ix_hi = max(0, ix_lo), min(nx, max(ix_lo + 2, ix_hi))
    iy_lo, iy_hi = max(0, iy_lo), min(ny, max(iy_lo + 2, iy_hi))
    Bpatch = data[iy_lo:iy_hi, ix_lo:ix_hi]
    lon_grid, lat_grid = np.meshgrid(lon_axis[ix_lo:ix_hi], lat_axis[iy_lo:iy_hi])
    Xg = (lon_grid - fp_lon) * deg2km / 1000.0  # Mm
    Yg = (lat_grid - fp_lat) * deg2km / 1000.0  # Mm
    Zg = np.zeros_like(Xg)

    # Seed grid for tracing inside the patch.
    lon_seeds = np.linspace(lon_lo, lon_hi, n_seeds_side)
    lat_seeds = np.linspace(lat_lo, lat_hi, n_seeds_side)
    LS, AS = np.meshgrid(lon_seeds, lat_seeds)
    lons_flat = LS.ravel()
    lats_flat = AS.ravel()
    r_seed = R_SUN_M * top_r  # ~ half a percent above photosphere; PFSS-resolvable

    obstime = output.input_map.date if hasattr(output, "input_map") else None
    seeds = SkyCoord(
        lons_flat * u.deg, lats_flat * u.deg, np.full_like(lons_flat, r_seed) * u.m,
        frame="heliographic_carrington", obstime=obstime, observer="earth",
    )
    tracer = make_tracer(tracing)
    fls = tracer.trace(seeds, output)

    # Convert traces to box-local (x_Mm, y_Mm, z_Mm) so the rendering matches a
    # Bifrost-style Cartesian box. We approximate the local tangent plane:
    # x ≈ (lon - fp_lon)*cos(fp_lat)*deg2km , y ≈ (lat - fp_lat)*deg2km , z = r - R_sun.
    cosfp = np.cos(np.deg2rad(fp_lat))
    lines_local = []
    line_colors = []
    cmap_T_like = cm.get_cmap("magma")
    z_lo, z_hi = 0.0, 8.0  # Mm — Noraz Fig 1's vertical extent for visual parity

    for fl in fls:
        try:
            # heliographic Carrington spherical -> get lon, lat, r at each step
            sph = fl.coords.represent_as("spherical")
            lon_deg = sph.lon.to_value(u.deg)
            lat_deg = sph.lat.to_value(u.deg)
            r_m = sph.distance.to_value(u.m)
        except Exception:
            continue
        if lon_deg.size < 2:
            continue
        # Unwrap longitude across 0/360 boundary.
        dl = np.diff(lon_deg)
        if np.any(np.abs(dl) > 180.0):
            lon_unw = np.unwrap(np.deg2rad(lon_deg))
            lon_deg = np.rad2deg(lon_unw)
        x_Mm = (lon_deg - fp_lon) * cosfp * deg2km / 1000.0
        y_Mm = (lat_deg - fp_lat) * deg2km / 1000.0
        z_Mm = (r_m - R_SUN_M) / 1.0e6
        xyz = np.stack([x_Mm, y_Mm, z_Mm], axis=-1)
        # Clip to the box.
        in_box = ((np.abs(xyz[:, 0]) <= half_Mm) &
                  (np.abs(xyz[:, 1]) <= half_Mm) &
                  (xyz[:, 2] >= z_lo) & (xyz[:, 2] <= z_hi))
        if np.count_nonzero(in_box) < 2:
            continue
        xyz = xyz[in_box]
        # Colour by mean z (analog to Noraz's "temperature increases with height").
        zmean = float(np.mean(xyz[:, 2]))
        color = cmap_T_like((zmean - z_lo) / max(z_hi - z_lo, 1e-9))
        lines_local.append(xyz)
        line_colors.append(color)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_box_aspect((1, 1, 0.6))

    # Base plane: photospheric B_LOS patch (clipped to ±200 G for contrast).
    norm = matplotlib.colors.Normalize(vmin=-200.0, vmax=200.0)
    cmap_B = cm.get_cmap("RdBu_r")
    facecolors = cmap_B(norm(np.clip(Bpatch, -200.0, 200.0)))
    ax.plot_surface(Xg, Yg, Zg, facecolors=facecolors, shade=False,
                    linewidth=0, antialiased=False, rstride=1, cstride=1, zorder=1)

    for xyz, color in zip(lines_local, line_colors):
        ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], color=color, lw=0.9, alpha=0.9, zorder=5)

    ax.set_xlim(-half_Mm, half_Mm)
    ax.set_ylim(-half_Mm, half_Mm)
    ax.set_zlim(z_lo, z_hi)
    ax.set_xlabel("x (Mm)")
    ax.set_ylabel("y (Mm)")
    ax.set_zlabel("z (Mm)")
    ax.view_init(elev=22, azim=-55)
    ax.set_title(
        f"Local PFSS volume around the chosen footpoint\n"
        f"(lon={fp_lon:.2f}°, lat={fp_lat:.2f}°, ±{box_half_deg}° ≈ ±{half_Mm:.1f} Mm; "
        f"{len(lines_local)} traces)",
        fontsize=10,
    )

    cax1 = fig.add_axes([0.08, 0.06, 0.36, 0.018])
    matplotlib.colorbar.ColorbarBase(
        cax1, cmap=cmap_B, norm=norm, orientation="horizontal",
        label="photospheric B$_{LOS}$ (G)",
    )
    cax2 = fig.add_axes([0.56, 0.06, 0.36, 0.018])
    matplotlib.colorbar.ColorbarBase(
        cax2, cmap=cmap_T_like,
        norm=matplotlib.colors.Normalize(vmin=z_lo, vmax=z_hi),
        orientation="horizontal",
        label="mean altitude along field line (Mm) — proxy for T in Noraz Fig 1",
    )

    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    print(f"[viz] wrote {out_path}")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--magnetogram", default="util/data/mrzqs190801t0014c2220_229.fits")
    ap.add_argument("--nrho", type=int, default=60)
    ap.add_argument("--rss",  type=float, default=2.5)
    ap.add_argument("--nseeds", type=int, default=96, help="seeds for global trace")
    ap.add_argument("--lat-window", type=float, nargs=2, default=[-60.0, 60.0],
                    help="lat range for seeds in global trace (degrees)")
    ap.add_argument("--b-min", type=float, default=20.0)
    ap.add_argument("--b-max", type=float, default=75.0)
    ap.add_argument("--fp-lat-window", type=float, nargs=2, default=[-15.0, 15.0])
    ap.add_argument("--fp-lon-window", type=float, nargs=2, default=[0.0, 360.0])
    ap.add_argument("--out-dir", default=os.path.join(HERE, "..", "visualization", "pfss"))
    ap.add_argument("--frames", type=int, default=72, help="rotation movie frame count")
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--no-movie", action="store_true")
    ap.add_argument("--box-half-deg", type=float, default=4.0,
                    help="half-width of the local box in degrees of arc")
    ap.add_argument("--n-local-seeds", type=int, default=12,
                    help="seeds per side in the local box (so n^2 traces)")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    warnings.filterwarnings("ignore", message=".*WCS.*")

    pfss, tracing, pkg = _import_pfsspy()
    print(f"[viz] using {pkg}")
    print(f"[viz] loading magnetogram {args.magnetogram}")
    smap = load_magnetogram(args.magnetogram)

    print(f"[viz] solving PFSS (nrho={args.nrho}, rss={args.rss})")
    pfss_in = pfss.Input(smap, args.nrho, args.rss)
    output = pfss.pfss(pfss_in)

    print(f"[viz] picking footpoint (same band as extract_field_line.py)")
    rng = np.random.default_rng(0)
    fp_lon, fp_lat, fp_B = pick_network_footpoint(
        smap, args.fp_lat_window, args.fp_lon_window, args.b_min, args.b_max, rng)
    print(f"[viz] footpoint: lon={fp_lon:.3f}°, lat={fp_lat:.3f}°, |B_LOS|={fp_B:.2f} G")

    print(f"[viz] tracing {args.nseeds} field lines for the global view")
    lons, lats = lonlat_grid_seeds(args.nseeds, lat_window=tuple(args.lat_window))
    # Seed slightly above the photosphere so the tracer can step.
    r_seed_m = R_SUN_M * 1.005
    lines, is_open = trace_seeds(output, lons, lats, r_seed_m, tracing)
    print(f"[viz] traced {len(lines)} usable lines ({sum(is_open)} open)")

    os.makedirs(args.out_dir, exist_ok=True)
    global_png = os.path.join(args.out_dir, "pfss_3d_global.png")
    local_png  = os.path.join(args.out_dir, "pfss_3d_local.png")
    rotate_mp4 = None if args.no_movie else os.path.join(args.out_dir, "pfss_3d_rotate.mp4")

    render_global(smap, lines, is_open, global_png,
                  rotate_mp4=rotate_mp4, fps=args.fps, frames=args.frames)
    render_local(smap, output, fp_lon, fp_lat, tracing, local_png,
                 box_half_deg=args.box_half_deg, n_seeds_side=args.n_local_seeds)

    print("[viz] done.")


if __name__ == "__main__":
    main()
