#!/usr/bin/env python3
"""
Extract a single field line from a GONG synoptic magnetogram via PFSS, sample
B(s) and (n_e, n_n, T) along it (from Model C7), and emit the tabulated data
file that scenarios/pfss_field_line.{hpp,cpp} reads.

Usage (from the repo root, with the venv active):
    python util/extract_field_line.py \
        --magnetogram util/data/mrzqs190801t0014c2220_229.fits \
        --output      scenarios/data/pfss_qs_20190801.dat \
        --ns          200

Footpoint strategy: scan the magnetogram in the chosen lat/lon window, find
pixels with |B_LOS| in the "network" band (50-150 G by default), and pick the
one whose |B_LOS| is closest to the band's median. This gives a representative
quiet-sun network flux-tube footpoint.

Caveat — also in the writeup: PFSS resolves r in [1, 2.5] R_sun. Our entire
1.5-Mm chromospheric domain occupies <0.3% of that range, so B(s) within the
chromosphere will be nearly constant. PFSS pins the footpoint |B| and the
field-line topology; chromospheric canopy expansion belongs to the
analytic_canopy scenario.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

# Constants (SI).
R_SUN_M = 6.957e8
G_SI    = 274.0
H_C7_BASE_KM = 1003.0  # base height of Model C7 table relative to tau=1

# Model C7 table, mirrored from scenarios/model_c7.cpp:46.
# Columns: height [km], n_e [m^-3], n_n [m^-3], T [K]. Keep in sync.
MODEL_C7 = np.array([
    [1.003e+03, 1.903e+17, 2.693e+19, 6.225e+03],
    [1.032e+03, 2.021e+17, 2.179e+19, 6.315e+03],
    [1.065e+03, 2.091e+17, 1.722e+19, 6.400e+03],
    [1.101e+03, 2.104e+17, 1.340e+19, 6.474e+03],
    [1.143e+03, 1.995e+17, 1.009e+19, 6.531e+03],
    [1.214e+03, 1.760e+17, 6.322e+18, 6.576e+03],
    [1.299e+03, 1.489e+17, 3.641e+18, 6.598e+03],
    [1.398e+03, 1.438e+17, 1.880e+18, 6.610e+03],
    [1.520e+03, 1.423e+17, 7.956e+17, 6.623e+03],
    [1.617e+03, 1.267e+17, 4.012e+17, 6.633e+03],
    [1.722e+03, 1.027e+17, 1.916e+17, 6.643e+03],
    [1.820e+03, 7.969e+16, 9.843e+16, 6.652e+03],
    [1.894e+03, 6.450e+16, 6.005e+16, 6.660e+03],
    [1.946e+03, 5.526e+16, 4.371e+16, 6.667e+03],
    [1.989e+03, 4.826e+16, 3.447e+16, 6.674e+03],
])


def _import_pfsspy():
    """Defensive import — package was renamed from pfsspy to sunkit-magex."""
    try:
        import sunkit_magex.pfss as pfss
        from sunkit_magex.pfss import tracing
        return pfss, tracing, "sunkit_magex"
    except ImportError:
        pass
    try:
        import pfsspy as pfss
        from pfsspy import tracing
        return pfss, tracing, "pfsspy"
    except ImportError:
        sys.stderr.write(
            "ERROR: neither sunkit-magex nor pfsspy is installed.\n"
            "  bash util/setup_venv.sh && source .venv/bin/activate\n"
        )
        sys.exit(1)


def load_magnetogram(path: str):
    """Read a GONG synoptic FITS map into a sunpy Map with the right WCS."""
    import sunpy.map
    smap = sunpy.map.Map(path)
    # GONG synoptic maps don't always carry the right WCS metadata for pfsspy;
    # patch the key fields if missing.
    md = smap.meta
    md.setdefault("ctype1", "CRLN-CAR")
    md.setdefault("ctype2", "CRLT-CAR")
    md.setdefault("cunit1", "deg")
    md.setdefault("cunit2", "deg")
    return sunpy.map.Map(smap.data, md)


def pick_network_footpoint(smap, lat_window, lon_window, b_min, b_max, rng):
    """Find a (carrington_lon_deg, lat_deg) inside the lat/lon window with
    |B_LOS| in [b_min, b_max] G, closest to the median of values in that band.
    Returns (lon_deg, lat_deg, b_value_G)."""
    import astropy.units as u
    data = np.asarray(smap.data, dtype=float)
    ny, nx = data.shape
    # GONG synoptic maps: x=Carrington longitude, y=sin(lat).
    lon = np.linspace(0.0, 360.0, nx, endpoint=False) + (180.0 / nx)
    sin_lat = np.linspace(-1.0, 1.0, ny, endpoint=False) + (1.0 / ny)
    lat = np.degrees(np.arcsin(sin_lat))
    LON, LAT = np.meshgrid(lon, lat)

    in_box = (LAT >= lat_window[0]) & (LAT <= lat_window[1]) \
           & (LON >= lon_window[0]) & (LON <= lon_window[1])
    absB = np.abs(data)
    in_band = (absB >= b_min) & (absB <= b_max) & in_box & np.isfinite(absB)
    if not np.any(in_band):
        raise RuntimeError(
            f"no pixels with |B_LOS| in [{b_min},{b_max}] G within "
            f"lat={lat_window} lon={lon_window} — relax the window or change the band."
        )
    band_vals = absB[in_band]
    target = float(np.median(band_vals))
    flat_idx = np.argmin(np.abs(absB - target) + np.where(in_band, 0.0, 1e30))
    iy, ix = np.unravel_index(flat_idx, absB.shape)
    return float(LON[iy, ix]), float(LAT[iy, ix]), float(data[iy, ix])


def probe_field_topology(output, lon_deg, lat_deg, tracing):
    """Run a PFSS trace from a seed ~0.01 R_sun above the footpoint (one
    radial PFSS cell up). Returns (is_open: bool, B_footpoint_T: float).
    The chromospheric domain itself (<2 Mm) sits entirely inside the first
    PFSS radial cell, so we cannot resolve B(s) in the chromosphere from
    PFSS alone — we just take the footpoint magnitude and report the
    topology.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    import sunpy.coordinates  # noqa: F401

    obstime = output.input_map.date if hasattr(output, "input_map") else None
    # Seed must be at least ~7 Mm above r=R_sun for the tracer to step.
    seed = SkyCoord(
        lon_deg * u.deg, lat_deg * u.deg, (R_SUN_M + 0.01 * R_SUN_M) * u.m,
        frame="heliographic_carrington", obstime=obstime, observer="earth",
    )
    if hasattr(tracing, "PerformanceTracer"):
        tracer = tracing.PerformanceTracer()
    elif hasattr(tracing, "FortranTracer"):
        tracer = tracing.FortranTracer()
    else:
        tracer = tracing.PythonTracer()
    fls = tracer.trace(seed, output)
    fl = fls[0]
    is_open = bool(getattr(fl, "is_open", True))

    # |B| at the lowest trace point (closest to the photosphere PFSS resolves).
    b_vec = getattr(fl, "b_along_fline", None)
    if b_vec is None or len(b_vec) == 0:
        # Fallback: take |B_LOS| from the input map at the footpoint.
        ny, nx = output.input_map.data.shape
        ix = int((lon_deg / 360.0) * nx) % nx
        iy = int((np.sin(np.deg2rad(lat_deg)) + 1.0) * 0.5 * ny)
        iy = max(0, min(ny - 1, iy))
        return is_open, abs(float(output.input_map.data[iy, ix])) * 1.0e-4
    b_vec_T = b_vec.to(u.T) if hasattr(b_vec, "unit") else u.Quantity(b_vec, u.G).to(u.T)
    b_mags = np.linalg.norm(b_vec_T.value, axis=-1).astype(float)
    # Find the trace point nearest r=R_sun.
    try:
        rs = np.linalg.norm(np.stack([
            fl.coords.cartesian.x.to_value(u.m),
            fl.coords.cartesian.y.to_value(u.m),
            fl.coords.cartesian.z.to_value(u.m),
        ], axis=-1), axis=-1)
        idx = int(np.argmin(np.abs(rs - R_SUN_M)))
    except Exception:
        idx = 0
    return is_open, float(b_mags[idx])


def build_chromospheric_grid(B_footpoint_T, ns, top_height_m):
    """Build a uniform-in-s grid of ns cells over [0, top_height_m] with B
    held constant at B_footpoint_T. This is the honest representation of what
    PFSS gives us in the chromosphere (one radial cell, no gradient
    information). The analytic_canopy scenario is the place to overlay
    canopy expansion if desired.

    Returns face arrays (length ns+1) and cell ds (length ns)."""
    s_faces = np.linspace(0.0, top_height_m, ns + 1)
    r_faces = R_SUN_M + s_faces  # vertical-field assumption
    b_faces = np.full_like(s_faces, B_footpoint_T)
    ds_cells = np.diff(s_faces)
    return s_faces, r_faces, b_faces, ds_cells


def interp_c7(h_m: np.ndarray):
    """Interpolate (n_e, n_n, T) onto height array h [m]. Below the C7 base
    we clamp to base values; above the top we exponentially extend density."""
    h_km = h_m / 1000.0
    h_tab = MODEL_C7[:, 0]
    ne_t  = MODEL_C7[:, 1]
    nn_t  = MODEL_C7[:, 2]
    T_t   = MODEL_C7[:, 3]
    # Clamp below; spline-like via numpy.interp (linear) is fine for the
    # 15-point table.
    ne = np.interp(h_km, h_tab, ne_t, left=ne_t[0], right=ne_t[-1])
    nn = np.interp(h_km, h_tab, nn_t, left=nn_t[0], right=nn_t[-1])
    T  = np.interp(h_km, h_tab, T_t,  left=T_t[0],  right=T_t[-1])
    # Above the top, exponentially extend n with H=200 km (rough chromospheric
    # scale height) so the upper ghost has sensible values.
    above = h_km > h_tab[-1]
    if np.any(above):
        z = (h_km[above] - h_tab[-1]) * 1000.0
        H = 200.0e3
        ne[above] = ne_t[-1] * np.exp(-z / H)
        nn[above] = nn_t[-1] * np.exp(-z / H)
    return ne, nn, T


def write_data_file(path, args, fp_lon, fp_lat, fp_B_G,
                    ds_cells, B_imh, B_iph,
                    phi_g_imh, phi_g_iph,
                    ne, nn, T,
                    ghost_outer, ghost_inner):
    ns = ds_cells.size
    with open(path, "w") as f:
        f.write("# pfss_field_line.dat — generated by util/extract_field_line.py\n")
        f.write(f"# Magnetogram: {os.path.basename(args.magnetogram)}\n")
        f.write(f"# Footpoint: lon={fp_lon:.3f} deg, lat={fp_lat:.3f} deg, "
                f"|B_LOS|={fp_B_G:.2f} G  (strategy=network, band=[{args.b_min},{args.b_max}] G)\n")
        f.write(f"# PFSS: rss={args.rss} R_sun, nrho={args.nrho}\n")
        f.write(f"# Resampled to ns={ns} cells in arc length s\n")
        f.write("# Units: SI throughout (m, T, K, m^-3, J/kg)\n")
        f.write("\n[META]\n")
        f.write(f"ns={ns}\n")
        f.write(f"g_si={G_SI}\n")
        f.write(f"B_outer_T={float(B_iph[-1]):.6e}\n")
        f.write(f"phi_g_offset_Jpkg=0.0\n")
        f.write("\n[CELLS]\n")
        f.write("# i  ds_m         B_imh_T      B_iph_T      "
                "phi_g_imh    phi_g_iph    ne_im3       nn_im3       T_K\n")
        for i in range(ns):
            f.write(f"{i}  {ds_cells[i]:.6e}  {B_imh[i]:.6e}  {B_iph[i]:.6e}  "
                    f"{phi_g_imh[i]:.6e}  {phi_g_iph[i]:.6e}  "
                    f"{ne[i]:.6e}  {nn[i]:.6e}  {T[i]:.6e}\n")
        f.write("\n[GHOSTS]\n")
        f.write("# tag      B_T          phi_g_Jpkg   ne_im3       nn_im3       T_K          T_e_factor\n")
        for tag, g in [("outer_0", ghost_outer[0]),
                       ("outer_1", ghost_outer[1]),
                       ("inner_0", ghost_inner[0]),
                       ("inner_1", ghost_inner[1])]:
            f.write(f"{tag}  {g['B']:.6e}  {g['phi_g']:.6e}  "
                    f"{g['ne']:.6e}  {g['nn']:.6e}  {g['T']:.6e}  {g['T_e_factor']:.3f}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--magnetogram", required=True, help="Path to GONG synoptic FITS")
    ap.add_argument("--output",      required=True, help="Path to write the .dat file")
    ap.add_argument("--ns",   type=int,   default=100)
    ap.add_argument("--rss",  type=float, default=2.5)
    ap.add_argument("--nrho", type=int,   default=60)
    # Default matches the Model C7 table span (h=1003 → 1989 km = 986 km of
    # domain), so the (n_e, n_n, T) interpolation stays inside C7's table.
    # Extend later when adding a corona-side region; for the chromospheric
    # comparison runs this is the natural choice.
    ap.add_argument("--top-height-km", type=float, default=986.0)
    ap.add_argument("--lat-window", type=float, nargs=2, default=[-15.0, 15.0])
    ap.add_argument("--lon-window", type=float, nargs=2, default=[0.0, 360.0])
    ap.add_argument("--b-min",      type=float, default=50.0)
    ap.add_argument("--b-max",      type=float, default=150.0)
    ap.add_argument("--seed",       type=int,   default=0)
    args = ap.parse_args()

    pfss, tracing, pkg = _import_pfsspy()
    rng = np.random.default_rng(args.seed)

    print(f"[extract] using {pkg}")
    print(f"[extract] reading magnetogram {args.magnetogram}")
    smap = load_magnetogram(args.magnetogram)

    print(f"[extract] picking network footpoint, |B_LOS| in [{args.b_min},{args.b_max}] G")
    fp_lon, fp_lat, fp_B_G = pick_network_footpoint(
        smap, args.lat_window, args.lon_window, args.b_min, args.b_max, rng)
    print(f"[extract] footpoint: lon={fp_lon:.3f}°, lat={fp_lat:.3f}°, |B_LOS|={fp_B_G:.2f} G")

    print(f"[extract] solving PFSS (rss={args.rss}, nrho={args.nrho})")
    pfss_in = pfss.Input(smap, args.nrho, args.rss)
    output = pfss.pfss(pfss_in)

    print(f"[extract] probing topology + footpoint |B| via PFSS")
    is_open, B_footpoint_T = probe_field_topology(output, fp_lon, fp_lat, tracing)
    print(f"[extract] topology: {'open' if is_open else 'CLOSED (loop)'} field line")
    print(f"[extract] footpoint |B| (PFSS-smoothed) = {B_footpoint_T:.3e} T "
          f"= {B_footpoint_T*1e4:.2f} G")
    print(f"[extract] (Caveat: PFSS does not resolve B(s) inside the chromosphere; "
          f"B is held uniform at the footpoint value across the {args.top_height_km}-km domain.)")

    # Build a uniform-in-s chromospheric grid with B = footpoint value.
    top_height_m = args.top_height_km * 1000.0
    s_faces, r_faces, b_faces, ds_cells = build_chromospheric_grid(
        B_footpoint_T, args.ns, top_height_m)
    B_imh = b_faces[:-1].astype(np.float32)
    B_iph = b_faces[1:].astype(np.float32)

    # Heights above the photosphere along the trace (so h = r - R_sun).
    h_faces_m = r_faces - R_SUN_M
    # φ_g(s) = g · h(s), gauge-zeroed at the footpoint.
    phi_g_faces = G_SI * (h_faces_m - h_faces_m[0])
    phi_g_imh = phi_g_faces[:-1].astype(np.float32)
    phi_g_iph = phi_g_faces[1:].astype(np.float32)

    # Cell-centered heights (for C7 interpolation): map photospheric base to
    # h=h_tab[0]=1003 km so the densities and temperatures land in the
    # tabulated range. I.e. add the C7 base offset to the trace heights.
    h_cell_m = 0.5 * (h_faces_m[:-1] + h_faces_m[1:]) + H_C7_BASE_KM * 1000.0
    ne, nn, T = interp_c7(h_cell_m)

    # Ghost cells.
    # Outer ghosts sit just past s_max; replicate the upper-end physics with
    # T_e_factor=2 to mimic the corona, matching scenarios/model_c7.cpp.
    h_outer_m = h_faces_m[-1] + ds_cells[-1] / 2.0 + H_C7_BASE_KM * 1000.0
    ne_o, nn_o, T_o = interp_c7(np.array([h_outer_m]))
    ghost_outer = [
        dict(B=float(B_iph[-1]), phi_g=float(phi_g_iph[-1]),
             ne=float(ne_o[0]), nn=float(nn_o[0]), T=float(T_o[0]), T_e_factor=2.0),
        dict(B=float(B_iph[-1]), phi_g=float(phi_g_iph[-1]),
             ne=float(ne_o[0]), nn=float(nn_o[0]), T=float(T_o[0]), T_e_factor=2.0),
    ]
    # Inner ghosts replicate cell 0.
    ghost_inner = [
        dict(B=float(B_imh[0]), phi_g=float(phi_g_imh[0]),
             ne=float(ne[0]), nn=float(nn[0]), T=float(T[0]), T_e_factor=1.0),
        dict(B=float(B_imh[0]), phi_g=float(phi_g_imh[0]),
             ne=float(ne[0]), nn=float(nn[0]), T=float(T[0]), T_e_factor=1.0),
    ]

    print(f"[extract] writing {args.output}")
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    write_data_file(args.output, args, fp_lon, fp_lat, fp_B_G,
                    ds_cells, B_imh, B_iph, phi_g_imh, phi_g_iph,
                    ne, nn, T, ghost_outer, ghost_inner)
    print(f"[extract] done. {args.ns} cells, footpoint |B|={fp_B_G:.1f} G, "
          f"top |B|={B_iph[-1]*1e4:.2f} G")


if __name__ == "__main__":
    main()
