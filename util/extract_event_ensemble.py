#!/usr/bin/env python3
"""Extract an ENSEMBLE of real PFSS field lines (closed AND open) around an active
region, build a runnable loop .dat for each, and emit (a) a run_ensemble.py job
manifest and (b) a lines index (with 3-D paths) for the 3-D sphere visualization.

This generalizes the single-line starter util/extract_event_loop.py (campaign
Phase P1 → ensemble). Each traced line becomes one independent chromo_main job:

  * CLOSED loop  → make_loop_dat.py --topology full  (footpoint A→apex→footpoint B,
                   reflecting walls). Sim half-length = loop_arc / 2.
  * OPEN line    → make_loop_dat.py --topology open   (vertical column, coronal
                   outflow BC at the top). The PFSS open line runs to the source
                   surface (~hundreds of Mm) which is far beyond where chromospheric
                   evaporation matters and where this model is valid, so we cap the
                   SIMULATED column at --open-cap-mm; the full PFSS path is still
                   recorded for the 3-D plot.

Selection: closed loops are filtered to flare-relevant lengths
[--len-min, --len-max] Mm (drops the giant overarching connections) and sampled
evenly across that length range; open lines are taken by descending footpoint |B|
(the strongest-rooted open fields near the AR).

Usage (2024-08-01 M8.2, AR 13768):
  python util/extract_event_ensemble.py \
    --magnetogram util/data/adapt40311_044012_202408010600_i00053600n1.fts.gz \
    --ar-lon 74 --ar-lat -16 --obstime 2024-08-01T07:09:00 \
    --n-closed 10 --n-open 4
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from extract_event_loop import (                       # noqa: E402
    R_SUN_M, _import_pfss, load_magnetogram, ar_to_carrington,
    make_tracer, field_line_length_m, footpoint_B_T,
)


def trace_cluster(output, car_lon, car_lat, patch_deg, n_seed, seed_h, tracing):
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    import sunpy.coordinates  # noqa: F401

    obstime = None
    try:
        obstime = output.input_map.date
    except Exception:
        pass
    grid = np.linspace(-patch_deg, patch_deg, n_seed)
    coslat = max(0.2, np.cos(np.deg2rad(car_lat)))
    lons = car_lon + grid / coslat
    lats = np.clip(car_lat + grid, -89.0, 89.0)
    LON, LAT = np.meshgrid(lons, lats)
    seeds = SkyCoord(LON.ravel() * u.deg, LAT.ravel() * u.deg,
                     seed_h * R_SUN_M * u.m,
                     frame="heliographic_carrington", obstime=obstime, observer="earth")
    tracer = make_tracer(tracing)
    return tracer.trace(seeds, output)


def pick_even(values, n):
    """Indices of n entries spread evenly across the sorted value range."""
    order = np.argsort(values)
    if len(order) <= n:
        return list(order)
    targets = np.linspace(0, len(order) - 1, n).round().astype(int)
    return [int(order[t]) for t in sorted(set(targets))]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--magnetogram", required=True)
    ap.add_argument("--ar-lon", type=float, required=True)
    ap.add_argument("--ar-lat", type=float, required=True)
    ap.add_argument("--obstime", required=True)
    ap.add_argument("--n-closed", type=int, default=10)
    ap.add_argument("--n-open", type=int, default=4)
    ap.add_argument("--len-min", type=float, default=15.0, help="min closed loop length [Mm]")
    ap.add_argument("--len-max", type=float, default=200.0, help="max closed loop length [Mm]")
    ap.add_argument("--open-cap-mm", type=float, default=80.0,
                    help="cap the SIMULATED open-column length [Mm]")
    ap.add_argument("--apex-T-MK", type=float, default=2.0)
    ap.add_argument("--rss", type=float, default=2.5)
    ap.add_argument("--nrho", type=int, default=60)
    ap.add_argument("--patch-deg", type=float, default=12.0)
    ap.add_argument("--n-seed", type=int, default=15)
    ap.add_argument("--seed-height-rsun", type=float, default=1.003)
    ap.add_argument("--dat-dir", default="scenarios/data/event_ensemble")
    ap.add_argument("--out-dir", default="outputs/event_ensemble")
    ap.add_argument("--manifest", default="outputs/event_ensemble_jobs.json")
    ap.add_argument("--lines-json", default="outputs/event_ensemble_lines.json")
    # Shared beam params for the ensemble (uniform for now; ribbon-mask per-line
    # weighting is Phase 2). M8.2 STIX/literature-anchored values.
    ap.add_argument("--flare-flux", default="3e7")
    ap.add_argument("--flare-delta", default="4")
    ap.add_argument("--flare-ecut", default="20")
    args = ap.parse_args()

    pfss, tracing, pkg = _import_pfss()
    print(f"[ensemble] using {pkg}; loading {os.path.basename(args.magnetogram)}")
    smap = load_magnetogram(args.magnetogram)
    car_lon, car_lat = ar_to_carrington(args.ar_lon, args.ar_lat, args.obstime)
    print(f"[ensemble] AR Carrington (lon={car_lon:.2f}, lat={car_lat:.2f})")
    print(f"[ensemble] PFSS rss={args.rss} nrho={args.nrho}")
    output = pfss.pfss(pfss.Input(smap, args.nrho, args.rss))

    flines = trace_cluster(output, car_lon, car_lat, args.patch_deg,
                           args.n_seed, args.seed_height_rsun, tracing)
    closed, openl = [], []
    for fl in flines:
        L_total, xyz = field_line_length_m(fl)
        if L_total <= 0 or not np.isfinite(L_total):
            continue
        rec = {"len_Mm": L_total / 1e6, "B_G": footpoint_B_T(fl) * 1e4, "xyz": xyz}
        (openl if bool(getattr(fl, "is_open", True)) else closed).append(rec)
    print(f"[ensemble] traced {len(flines)} seeds -> {len(closed)} closed, {len(openl)} open")

    # --- select ---------------------------------------------------------------
    cl_ok = [c for c in closed if args.len_min <= c["len_Mm"] <= args.len_max]
    cl_sel = [cl_ok[i] for i in pick_even([c["len_Mm"] for c in cl_ok], args.n_closed)]
    op_sel = sorted(openl, key=lambda o: -o["B_G"])[:args.n_open]
    print(f"[ensemble] selected {len(cl_sel)} closed (len {args.len_min}-{args.len_max} Mm) "
          f"+ {len(op_sel)} open")

    os.makedirs(args.dat_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.manifest) or ".", exist_ok=True)

    jobs, lines = [], []

    def build_line(kind, idx, rec):
        topo_sim = "full" if kind == "closed" else "open"
        if kind == "closed":
            half = rec["len_Mm"] / 2.0
            length_arg = half
            sim_len_Mm = rec["len_Mm"]
        else:
            length_arg = min(rec["len_Mm"], args.open_cap_mm)  # capped column
            sim_len_Mm = length_arg
        name = f"{kind[0]}{idx:02d}_{rec['len_Mm']:.0f}Mm"
        dat = os.path.join(args.dat_dir, f"{name}.dat")
        out = os.path.join(args.out_dir, f"{name}.txt")
        cmd = [sys.executable, os.path.join(HERE, "make_loop_dat.py"),
               "--topology", topo_sim, "--half-length-mm", f"{length_arg:.3f}",
               "--fp-B-G", f"{rec['B_G']:.1f}", "--apex-T-MK", f"{args.apex_T_MK}",
               "--output", dat]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  [skip] {name}: make_loop_dat failed\n{r.stderr[-400:]}")
            return
        jobs.append({"name": name, "dat": dat, "out": out,
                     "env": {"FLARE_BEAM_FLUX": args.flare_flux}})
        lines.append({
            "name": name, "kind": kind, "topology": topo_sim,
            "loop_total_arc_Mm": round(rec["len_Mm"], 2),
            "sim_length_Mm": round(sim_len_Mm, 2),
            "footpoint_B_G": round(rec["B_G"], 1),
            "out": out, "dat": dat,
            "path_xyz_rsun": (rec["xyz"] / R_SUN_M).round(5).tolist(),
        })
        print(f"  [{name}] {kind} sim_len={sim_len_Mm:.0f}Mm B={rec['B_G']:.0f}G "
              f"({rec['xyz'].shape[0]} path pts)")

    for i, rec in enumerate(cl_sel):
        build_line("closed", i, rec)
    for i, rec in enumerate(op_sel):
        build_line("open", i, rec)

    manifest = {
        "binary": "build/chromo_main",
        "defaults": {
            "mode": "full", "ioniz": "ionization", "scenario": "pfss_field_line",
            "time_mult": 1.0, "cooling": "cooling",
            "env": {"PFSS_FLARE": "1", "FLARE_DELTA": args.flare_delta,
                    "FLARE_E_CUT": args.flare_ecut},
        },
        "jobs": jobs,
    }
    with open(args.manifest, "w") as f:
        json.dump(manifest, f, indent=2)
    with open(args.lines_json, "w") as f:
        json.dump({"event": "2024-08-01 M8.2 AR13768",
                   "ar_carrington_lon_deg": round(car_lon, 3),
                   "ar_carrington_lat_deg": round(car_lat, 3),
                   "magnetogram": os.path.basename(args.magnetogram),
                   "lines": lines}, f)
    print(f"[ensemble] wrote {len(jobs)} jobs -> {args.manifest}")
    print(f"[ensemble] wrote lines index -> {args.lines_json}")
    print(f"[ensemble] run: python util/run_ensemble.py --manifest {args.manifest}")


if __name__ == "__main__":
    main()
