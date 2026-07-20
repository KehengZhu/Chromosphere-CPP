#!/usr/bin/env python3
"""Extract ONE flare-loop geometry from a real-event magnetogram via PFSS.

This is the campaign "starter" (Phase P1, single line): take a global magnetogram
(ADAPT-GONG ensemble map, or a GONG synoptic map), solve PFSS, seed a cluster of
field lines around a target active region, pick a representative CLOSED loop, and
report exactly the three numbers make_loop_dat.py needs to build a runnable
chromosphere→TR→corona loop:

    * loop half-length  L  [Mm]   (footpoint→apex arc length = total_arc / 2)
    * topology                    (closed loop vs open)
    * footpoint |B|     [G]       (the flaring-footpoint field strength)

PFSS is current-free, so it gives connectivity + geometry + footpoint |B| ONLY —
not the CME, not the non-potential flaring core (see the writeup). We use it to
pin the loop the chromospheric model runs ON; the pre-flare stratification is the
standard C7 atmosphere that make_loop_dat.py lays down.

The target AR is specified in Stonyhurst heliographic coords (as reported by HEK)
plus the observation time; we convert to Carrington internally with astropy so
there is no hand-rolled L0 arithmetic.

Usage (2024-08-01 07:09 UT M8.2, AR 13768 at Stonyhurst lon=+74, lat=-16):
    python util/extract_event_loop.py \
        --magnetogram util/data/adapt40311_044012_202408010600_i00005600n1.fts.gz \
        --ar-lon 74 --ar-lat -16 --obstime 2024-08-01T07:09:00 \
        --out-json outputs/events/event_20240801_AR13768_geometry.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

R_SUN_M = 6.957e8


def _import_pfss():
    try:
        import sunkit_magex.pfss as pfss
        from sunkit_magex.pfss import tracing
        return pfss, tracing, "sunkit_magex"
    except ImportError:
        pass
    import pfsspy as pfss
    from pfsspy import tracing
    return pfss, tracing, "pfsspy"


def _date_from_filename(path):
    """Pull a YYYYMMDDHHMM stamp out of an ADAPT filename → ISO date string."""
    import re
    m = re.search(r"_(\d{12})_", os.path.basename(path))
    if not m:
        return None
    s = m.group(1)
    return f"{s[0:4]}-{s[4:6]}-{s[6:8]}T{s[8:10]}:{s[10:12]}:00"


def load_magnetogram(path):
    """Return a flux-balanced sunpy Map in Carrington **CEA** coords (pfsspy needs
    sine-latitude / CEA, not the plate-carrée grid ADAPT ships).

    Handles ADAPT ensemble cubes (12 realizations → take realization 0, the
    standard single-map choice; CAR grid → reprojected to CEA) and ordinary 2-D
    GONG synoptic maps (already CEA).
    """
    import sunpy.map
    from astropy.io import fits
    from sunkit_magex.pfss import utils as pfss_utils

    with fits.open(path) as hdul:
        hdul.verify("silentfix")
        data = hdul[0].data
        header = hdul[0].header
    if data is None:
        with fits.open(path) as hdul:
            data = hdul[1].data
            header = hdul[1].header

    is_adapt = (data.ndim == 3) or str(header.get("MODEL", "")).upper().startswith("ADAPT")
    if data.ndim == 3:                       # ADAPT ensemble cube (12,nlat,nlon)
        slab = np.asarray(data[0, :, :], dtype=float)
        print(f"[extract] ADAPT ensemble cube {data.shape}; using realization 0")
    else:
        slab = np.asarray(data, dtype=float)

    if is_adapt:
        # ADAPT is plate-carrée (uniform lon × uniform lat). Build a clean
        # Carrington-CAR map, then reproject to CEA for PFSS.
        ny, nx = slab.shape
        date = _date_from_filename(path) or "2024-01-01T00:00:00"
        meta = {
            "wcsaxes": 2,
            "ctype1": "CRLN-CAR", "ctype2": "CRLT-CAR",
            "cunit1": "deg", "cunit2": "deg",
            "crval1": 180.0, "crval2": 0.0,
            "cdelt1": 360.0 / nx, "cdelt2": 180.0 / ny,
            "crpix1": (nx + 1) / 2.0, "crpix2": (ny + 1) / 2.0,
            "date-obs": date,
            "naxis": 2, "naxis1": nx, "naxis2": ny,
        }
        car_map = sunpy.map.Map(slab, meta)
        smap = pfss_utils.car_to_cea(car_map)
    else:
        # GONG synoptic: 2-D, already CEA. Patch any missing unit keys.
        hdr = header.copy()
        for key in list(hdr.keys()):
            if key.endswith("3"):
                del hdr[key]
        hdr["NAXIS"] = 2
        hdr["NAXIS1"] = slab.shape[1]
        hdr["NAXIS2"] = slab.shape[0]
        hdr.setdefault("cunit1", "deg")
        hdr.setdefault("cunit2", "deg")
        smap = sunpy.map.Map(slab, hdr)

    # pfsspy wants flux balance; ADAPT maps are close but not exact. Remove the
    # residual mean so the monopole term is zero (standard pre-processing).
    arr = np.asarray(smap.data, dtype=float)
    arr = arr - np.nanmean(arr)
    return sunpy.map.Map(arr, smap.meta)


def ar_to_carrington(ar_lon_deg, ar_lat_deg, obstime):
    """Convert Stonyhurst heliographic (lon,lat) at obstime → Carrington (lon,lat)."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    import sunpy.coordinates  # noqa: F401

    c = SkyCoord(ar_lon_deg * u.deg, ar_lat_deg * u.deg,
                 frame="heliographic_stonyhurst", obstime=obstime,
                 observer="earth")
    cc = c.transform_to("heliographic_carrington")
    return float(cc.lon.to_value(u.deg)) % 360.0, float(cc.lat.to_value(u.deg))


def make_tracer(tracing):
    for name in ("PerformanceTracer", "FortranTracer", "PythonTracer"):
        if hasattr(tracing, name):
            cls = getattr(tracing, name)
            try:
                return cls(max_steps=20000)   # AR loops are long; avoid truncation
            except TypeError:
                return cls()
    raise RuntimeError("no usable pfsspy tracer found")


def field_line_length_m(fl):
    """Arc length [m] of a traced field line from its 3-D Carrington coords."""
    import astropy.units as u
    xyz = np.stack([
        fl.coords.cartesian.x.to_value(u.m),
        fl.coords.cartesian.y.to_value(u.m),
        fl.coords.cartesian.z.to_value(u.m),
    ], axis=-1)
    seg = np.linalg.norm(np.diff(xyz, axis=0), axis=-1)
    return float(seg.sum()), xyz


def footpoint_B_T(fl):
    """Max |B| [T] along the line (the strong photospheric footpoint).

    pfsspy returns b_along_fline in the input map's units (Gauss for ADAPT/GONG),
    as a dimensionless array, so we interpret the magnitude as Gauss → Tesla."""
    b = getattr(fl, "b_along_fline", None)
    if b is None or len(b) == 0:
        return float("nan")
    arr = np.asarray(getattr(b, "value", b), dtype=float)   # strip any (dimensionless) unit
    bmag = np.linalg.norm(arr, axis=-1)
    return float(np.nanmax(bmag)) * 1.0e-4                   # Gauss -> Tesla


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--magnetogram", required=True)
    ap.add_argument("--ar-lon", type=float, required=True,
                    help="AR Stonyhurst longitude [deg] (+W) at --obstime")
    ap.add_argument("--ar-lat", type=float, required=True,
                    help="AR Stonyhurst latitude [deg] (+N) at --obstime")
    ap.add_argument("--obstime", required=True, help="ISO time, e.g. 2024-08-01T07:09:00")
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--rss", type=float, default=2.5)
    ap.add_argument("--nrho", type=int, default=60)
    ap.add_argument("--patch-deg", type=float, default=6.0,
                    help="half-width of the seed cluster around the AR [deg]")
    ap.add_argument("--n-seed", type=int, default=7, help="seeds per axis (n_seed^2 total)")
    ap.add_argument("--seed-height-rsun", type=float, default=1.003,
                    help="r/R_sun to seed the tracer at (just above the photosphere)")
    args = ap.parse_args()

    pfss, tracing, pkg = _import_pfss()
    print(f"[extract] using {pkg}")
    print(f"[extract] loading magnetogram {os.path.basename(args.magnetogram)}")
    smap = load_magnetogram(args.magnetogram)

    car_lon, car_lat = ar_to_carrington(args.ar_lon, args.ar_lat, args.obstime)
    print(f"[extract] AR Stonyhurst (lon={args.ar_lon:+.1f}, lat={args.ar_lat:+.1f}) "
          f"@ {args.obstime}  ->  Carrington (lon={car_lon:.2f}, lat={car_lat:.2f})")

    print(f"[extract] solving PFSS (rss={args.rss} R_sun, nrho={args.nrho})")
    pfss_in = pfss.Input(smap, args.nrho, args.rss)
    output = pfss.pfss(pfss_in)

    # Seed a cluster of points around the AR footpoint, just above the photosphere.
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    import sunpy.coordinates  # noqa: F401

    obstime = smap.date
    half = args.patch_deg
    grid = np.linspace(-half, half, args.n_seed)
    # Shrink the longitude span by cos(lat) so the patch is ~square on the sphere.
    coslat = max(0.2, np.cos(np.deg2rad(car_lat)))
    lons = (car_lon + grid / coslat)
    lats = np.clip(car_lat + grid, -89.0, 89.0)
    LON, LAT = np.meshgrid(lons, lats)
    seeds = SkyCoord(LON.ravel() * u.deg, LAT.ravel() * u.deg,
                     args.seed_height_rsun * R_SUN_M * u.m,
                     frame="heliographic_carrington", obstime=obstime, observer="earth")

    tracer = make_tracer(tracing)
    print(f"[extract] tracing {seeds.shape[0]} seeds over a "
          f"±{half:.0f}° patch around the AR")
    flines = tracer.trace(seeds, output)

    closed = []
    for fl in flines:
        is_open = bool(getattr(fl, "is_open", True))
        if is_open:
            continue
        L_total, xyz = field_line_length_m(fl)
        if L_total <= 0 or not np.isfinite(L_total):
            continue
        closed.append((L_total, footpoint_B_T(fl), xyz))

    n_open = len(flines) - len(closed)
    print(f"[extract] traced: {len(closed)} closed, {n_open} open/failed")
    if not closed:
        sys.stderr.write(
            "ERROR: no closed loops found near the AR. Widen --patch-deg, lower the "
            "seed height, or check the AR coordinates / source-surface radius.\n")
        sys.exit(2)

    lengths = np.array([c[0] for c in closed])
    # Representative loop: the median-length closed line (robust to short/grazing
    # field lines and to the one giant overarching connection).
    order = np.argsort(lengths)
    pick = closed[order[len(order) // 2]]
    L_total, B_fp_T, xyz = pick
    L_half_Mm = (L_total / 2.0) / 1.0e6
    B_fp_G = B_fp_T * 1.0e4

    print(f"[extract] closed-loop length distribution [Mm]: "
          f"min={lengths.min()/1e6:.1f}, median={np.median(lengths)/1e6:.1f}, "
          f"max={lengths.max()/1e6:.1f}  (n={len(lengths)})")
    print(f"[extract] SELECTED representative loop:")
    print(f"            total arc length = {L_total/1e6:.2f} Mm  ->  half-length L = {L_half_Mm:.2f} Mm")
    print(f"            topology         = closed")
    print(f"            footpoint |B|    = {B_fp_G:.1f} G")
    apex_r = np.linalg.norm(xyz, axis=-1).max() / R_SUN_M
    print(f"            apex height      = {(apex_r-1.0)*R_SUN_M/1e6:.1f} Mm above photosphere")

    geom = {
        "event": "2024-08-01 ~07:09 UT M8.2",
        "AR": "NOAA 13768",
        "magnetogram": os.path.basename(args.magnetogram),
        "ar_stonyhurst_lon_deg": args.ar_lon,
        "ar_stonyhurst_lat_deg": args.ar_lat,
        "ar_carrington_lon_deg": round(car_lon, 3),
        "ar_carrington_lat_deg": round(car_lat, 3),
        "obstime": args.obstime,
        "pfss_rss_rsun": args.rss,
        "pfss_nrho": args.nrho,
        "n_closed": int(len(lengths)),
        "n_open_or_failed": int(n_open),
        "loop_total_arc_Mm": round(L_total / 1e6, 3),
        "loop_half_length_Mm": round(L_half_Mm, 3),
        "topology": "closed",
        "footpoint_B_G": round(B_fp_G, 2),
        "apex_height_Mm": round((apex_r - 1.0) * R_SUN_M / 1e6, 3),
        "path_xyz_rsun": (xyz / R_SUN_M).round(5).tolist(),
    }
    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w") as f:
        json.dump(geom, f, indent=2)
    print(f"[extract] wrote geometry -> {args.out_json}")
    print(f"[extract] next: python util/make_loop_dat.py --topology full "
          f"--half-length-mm {L_half_Mm:.1f} --fp-B-G {B_fp_G:.0f} "
          f"--apex-T-MK 2.0 --output scenarios/data/event_20240801_AR13768.dat")


if __name__ == "__main__":
    main()
