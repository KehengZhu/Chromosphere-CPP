#!/usr/bin/env python3
"""Resolution-convergence scan for the model_column top-region rho*V ripple.

Answers one question: does the cell-centred mass-flux ripple in 2130-2150 km
vanish as Delta h -> 0, while the *effective* finite-volume face flux stays
smooth and conservative?

All face quantities come from the `<run>.faceflux` sidecar, i.e. a read-only
copy-out of the PRODUCTION `mixture_rhs_explicit` arrays (see
`util/face_flux_diag.py` and `docs/top_ripple_face_flux_diagnosis.md`). Nothing
here re-derives the reconstruction, the limiter or the Riemann solver.

Definitions (index i is the UPPER face i+1/2 of cell i; ds = Delta h):

    M_i        = rho_cell[i] * v_cell[i]              cell-centred mass flux
    f_central  = 0.5[(rho v)_L + (rho v)_R]
    f_diff     = -0.5 a (rho_R - rho_L)
    f_total    = f_central + f_diff                   the production face flux
    f_ref      = f_total at t = 0                     t = 0 baseline flux
    f_eff      = f_total - f_ref                      CHANGE in the flux since t=0
    d_diff     = f_diff - f_diff(ref)

    A_x  = std(x) / |mean(f_eff)|      (for x = d_diff)  or  std(x)/|mean(x)|
                                       (for x = M, f_eff)
    Q_x  = mean|D2 x| / mean|x|        normalised curvature ("roughness")
    R_eff[i] = -(f_eff[i] - f_eff[i-1]) / ds

Physical/numerical mismatch. The two live on different staggered locations, so
one has to be mapped. We map the FACE flux DOWN to cell centres with the
symmetric two-point average

    f_eff_cell[i] = 0.5 ( f_eff[i-1] + f_eff[i] )      (lower + upper face of i)

which is second-order accurate for a smooth field, hence cannot manufacture a
first-order signal. Then

    A_mismatch = RMS[ M_i - f_eff_cell[i] ] / |mean(f_eff)|

evaluated over the physical window only (which is strictly interior to the
capture, so i-1 always exists).

The comparison window is a fixed PHYSICAL height band, not a fixed cell count,
so the three resolutions are compared over the same 20 km of atmosphere.

Usage:
    python util/resolution_scan_diag.py \
        --run 500=outputs/model_column/resconv_ns500_dec_1000s.txt \
        --run 1000=... --run 2000=... \
        --times 500 1000 --window 2130 2150 \
        --json out.json --report out.txt --plot-dir visualization/model_column
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from face_flux_diag import COLUMNS, read_faceflux, relstd, roughness  # noqa: E402
from upper_bc_decouple_diag import COLUMNS as GCOLUMNS, load_gamma_diag  # noqa: E402


# ----------------------------------------------------------------- record pick
def pick_record(records, t_target):
    times = np.array([r[0] for r in records])
    k = int(np.argmin(np.abs(times - t_target)))
    return records[k]


# ------------------------------------------------------------------- metrics
def face_metrics(record, ref, ds_km, win):
    """All ripple/flux metrics for one (run, time) over the physical window."""
    t, step, c = record
    ds = ds_km * 1.0e3
    km_cell = c["cell_km"]
    km_face = c["face_km"]

    M = c["rho_cell"] * c["v_cell"]
    f_tot, f_cen, f_dif = c["f_total"], c["f_central"], c["f_diff"]
    f_eff = f_tot - ref["f_total"]
    d_cen = f_cen - ref["f_central"]
    d_dif = f_dif - ref["f_diff"]

    # face -> cell mapping (see module docstring). index 0 has no lower face.
    f_eff_cell = np.full_like(f_eff, np.nan)
    f_eff_cell[1:] = 0.5 * (f_eff[:-1] + f_eff[1:])

    R_eff = np.full_like(f_eff, np.nan)
    R_eff[1:] = -np.diff(f_eff) / ds
    R_tot = np.full_like(f_tot, np.nan)
    R_tot[1:] = -np.diff(f_tot) / ds

    mc = (km_cell >= win[0]) & (km_cell <= win[1])
    mf = (km_face >= win[0]) & (km_face <= win[1])
    if mc.sum() < 5 or mf.sum() < 5:
        raise SystemExit("window %s too narrow for the capture" % (win,))

    Mw, effw = M[mc], f_eff[mf]
    mean_eff = float(np.mean(effw))
    scale = abs(mean_eff)

    # mismatch on the cell window (strictly interior, so f_eff_cell is defined)
    mm = mc & np.isfinite(f_eff_cell)
    resid = M[mm] - f_eff_cell[mm]

    out = {
        "t": float(t), "step": int(step), "ds_km": ds_km,
        "n_cells_window": int(mc.sum()), "n_faces_window": int(mf.sum()),
        # 1. cell-centred physical mass flux
        "M_mean": float(np.mean(Mw)),
        "A_M": relstd(Mw),
        "Q_M": roughness(Mw),
        # 2. effective finite-volume flux
        "f_eff_mean": mean_eff,
        "A_eff": relstd(effw),
        "Q_eff": roughness(effw),
        "f_total_mean": float(np.mean(f_tot[mf])),
        "A_total": relstd(f_tot[mf]),
        "Q_total": roughness(f_tot[mf]),
        # 3. reference-subtracted diffusive contribution
        "A_diff": float(np.std(d_dif[mf]) / max(scale, 1e-300)),
        "A_central": float(np.std(d_cen[mf]) / max(scale, 1e-300)),
        "corr_dcentral_ddiff": (float(np.corrcoef(d_cen[mf], d_dif[mf])[0, 1])
                                if d_cen[mf].std() > 0 and d_dif[mf].std() > 0
                                else float("nan")),
        "f_diff_mean": float(np.mean(f_dif[mf])),
        "a_over_V_median": float(np.median(
            c["a_face"][mf] / np.maximum(np.abs(c["v_cell"][mf]), 1e-30))),
        # 4. physical / numerical mismatch
        "A_mismatch": float(np.sqrt(np.mean(resid ** 2)) / max(scale, 1e-300)),
        "mismatch_rms": float(np.sqrt(np.mean(resid ** 2))),
        "mismatch_absmax": float(np.max(np.abs(resid))),
        # 5. continuity residual of the effective flux
        "R_eff_rms": float(np.sqrt(np.mean(R_eff[mf][1:] ** 2))) if mf.sum() > 1 else float("nan"),
        "R_eff_absmax": float(np.max(np.abs(R_eff[mf][1:]))) if mf.sum() > 1 else float("nan"),
        "R_eff_rms_norm": float(np.sqrt(np.mean(R_eff[mf][1:] ** 2)) / (scale / ds)),
        "R_total_rms": float(np.sqrt(np.mean(R_tot[mf][1:] ** 2))) if mf.sum() > 1 else float("nan"),
        "R_total_rms_norm": float(np.sqrt(np.mean(R_tot[mf][1:] ** 2)) / (scale / ds)),
        # sanity
        "split_max_relerr": float(np.max(
            np.abs(f_cen + f_dif - f_tot) / np.maximum(np.abs(f_tot), 1e-300))),
    }
    out.update(ripple_geometry(km_cell[mc], Mw, ds_km))
    return out


def ripple_geometry(km, x, ds_km):
    """Extrema positions (km), plus two independent wavelength estimates.

    `lambda_extrema` = 2 x mean spacing between consecutive local extrema (a full
    wavelength spans two extrema). Cheap but biased low when the pattern is close
    to the grid Nyquist limit (2 cells), so it under-reports on coarse meshes.

    `lambda_peak` = period of the largest peak of the periodogram of the
    linearly-detrended window. Independent of the extremum count.
    """
    idx = [i for i in range(1, len(x) - 1)
           if (x[i] - x[i - 1]) * (x[i + 1] - x[i]) < 0.0]
    pos = [float(km[i]) for i in idx]
    if len(pos) >= 2:
        lam_km = float(2.0 * np.mean(np.diff(pos)))
    else:
        lam_km = float("nan")

    # periodogram of the detrended window (uniform mesh, so km spacing is ds_km)
    n = len(x)
    lam_peak = float("nan")
    if n >= 8:
        d = x - np.polyval(np.polyfit(np.arange(n), x, 1), np.arange(n))
        d = d - d.mean()
        p = np.abs(np.fft.rfft(d * np.hanning(n))) ** 2
        k = int(np.argmax(p[1:])) + 1                # skip DC
        lam_peak = float(n * ds_km / k)
    return {
        "extrema_km": pos,
        "n_extrema": len(pos),
        "lambda_km": lam_km,
        "lambda_cells": lam_km / ds_km if np.isfinite(lam_km) else float("nan"),
        "lambda_peak_km": lam_peak,
        "lambda_peak_cells": lam_peak / ds_km if np.isfinite(lam_peak) else float("nan"),
    }


# ---------------------------------------------------------------- global state
def global_metrics(path, t_target, t_wall, win):
    """Column-integral / boundary quantities from the .gamma_diag sidecar."""
    h, frames = load_gamma_diag(path)
    ds_km = float(np.mean(np.diff(h)))
    ds = ds_km * 1.0e3
    iT, iK, iP = GCOLUMNS.index("T"), GCOLUMNS.index("kappa_solver"), GCOLUMNS.index("p")

    t0, a0 = frames[0]
    M0 = float(np.sum(a0[:, 0]) * ds)

    # cumulative conductive input, trapezoid over all frames up to t_target
    ts = np.array([f[0] for f in frames])
    qs = np.array([f[1][-1, iK] * (t_wall - f[1][-1, iT]) / ds for f in frames])
    keep = ts <= t_target + 1e-6
    E_cond = float(np.trapezoid(qs[keep], ts[keep])) if keep.sum() > 1 else 0.0

    k = int(np.argmin(np.abs(ts - t_target)))
    t, arr = frames[k]
    rho, v = arr[:, 0], arr[:, 1]
    band = (h >= win[0]) & (h <= win[1])
    return {
        "t_frame": float(t), "ds_km": ds_km, "ns": len(h),
        "M_column": float(np.sum(rho) * ds),
        "M_rel_change": float(np.sum(rho) * ds / M0 - 1.0),
        "rhoV_window_mean": float(np.mean((rho * v)[band])),
        "rhoV_top": float(rho[-1] * v[-1]),
        "V_top": float(v[-1]),
        "V_Nm2": float(v[-2]),
        "T_top": float(arr[-1, iT]),
        "p_top": float(arr[-1, iP]),
        "q_top": float(arr[-1, iK] * (t_wall - arr[-1, iT]) / ds),
        "E_cond": E_cond,
    }


# ------------------------------------------------------------------ fit helper
def fit_order(ds, y):
    ds, y = np.asarray(ds, float), np.asarray(y, float)
    ok = np.isfinite(y) & (y > 0)
    if ok.sum() < 2:
        return float("nan")
    return float(np.polyfit(np.log(ds[ok]), np.log(y[ok]), 1)[0])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, metavar="NS=RUNPATH",
                    help="ns label = path to the run .txt (sidecars inferred)")
    ap.add_argument("--times", nargs="+", type=float, default=[500.0, 1000.0])
    ap.add_argument("--window", nargs=2, type=float, default=(2130.0, 2150.0))
    ap.add_argument("--t-wall", type=float, default=22000.0)
    ap.add_argument("--p-back", type=float, default=None,
                    help="kOuterPRef [Pa] for the p_top - p_back column")
    ap.add_argument("--json", default=None)
    ap.add_argument("--report", default=None)
    ap.add_argument("--plot-dir", default=None)
    ap.add_argument("--plot-prefix", default="resconv")
    args = ap.parse_args(argv)

    win = tuple(args.window)
    runs = []
    for spec in args.run:
        label, _, path = spec.partition("=")
        ff = path if path.endswith(".faceflux") else path + ".faceflux"
        meta, records = read_faceflux(ff)
        if len(records) < 2:
            raise SystemExit(f"{ff}: need the t=0 reference plus at least one record")
        ds_km = float(meta["ds_km"])
        ref = records[0][2]
        if abs(records[0][0]) > 1e-9:
            raise SystemExit(f"{ff}: first record is t={records[0][0]}, not 0")
        entry = {"label": label, "ns": int(meta["ns"]), "ds_km": ds_km,
                 "path": path, "faceflux": ff, "meta": meta,
                 "n_records": len(records), "times": {}, "series": {}}
        for tt in args.times:
            rec = pick_record(records, tt)
            m = face_metrics(rec, ref, ds_km, win)
            # use the FACE record's own time so the two sidecars describe the
            # same instant as closely as their cadences allow
            g = global_metrics(path, rec[0], args.t_wall, win)
            if args.p_back is not None:
                g["p_top_minus_p_back"] = g["p_top"] - args.p_back
            m["global"] = g
            entry["times"][f"{tt:g}"] = m
            t, step, c = rec
            entry["series"][f"{tt:g}"] = {
                "t": float(t),
                "cell_km": c["cell_km"].tolist(),
                "face_km": c["face_km"].tolist(),
                "rhoV": (c["rho_cell"] * c["v_cell"]).tolist(),
                "f_eff": (c["f_total"] - ref["f_total"]).tolist(),
                "d_diff": (c["f_diff"] - ref["f_diff"]).tolist(),
            }
        runs.append(entry)

    runs.sort(key=lambda r: -r["ds_km"])
    result = {"window": list(win), "times": args.times, "runs": runs,
              "convergence": {}}

    core = ["A_M", "Q_M", "A_eff", "Q_eff", "A_diff", "A_central",
            "A_mismatch", "R_eff_rms", "R_eff_rms_norm", "M_mean", "f_eff_mean"]
    ds = [r["ds_km"] for r in runs]
    for tt in args.times:
        key = f"{tt:g}"
        conv = {}
        for name in core:
            y = [r["times"][key][name] for r in runs]
            conv[name] = {
                "values": y,
                "ratio_coarse_to_mid": (y[0] / y[1]) if len(y) > 1 and y[1] else float("nan"),
                "ratio_mid_to_fine": (y[1] / y[2]) if len(y) > 2 and y[2] else float("nan"),
                "order_p": fit_order(ds, np.abs(y)),
            }
        # what Delta h would be needed for a given ripple amplitude, extrapolated
        # from (a) the 3-point global fit and (b) the LOCAL mid->fine rate
        y_fine = abs(runs[-1]["times"][key]["A_M"])
        dh_fine = runs[-1]["ds_km"]
        p_glob = conv["A_M"]["order_p"]
        r_local = conv["A_M"]["ratio_mid_to_fine"]
        p_loc = (np.log(r_local) / np.log(runs[1]["ds_km"] / runs[2]["ds_km"])
                 if len(runs) > 2 and np.isfinite(r_local) and r_local > 0 else float("nan"))
        req = {"p_global": p_glob, "p_local_mid_to_fine": float(p_loc)}
        DH_KM = 553.0
        for target in (0.05, 0.02):
            for tag, p in (("global", p_glob), ("local", p_loc)):
                if np.isfinite(p) and p > 0:
                    dh = dh_fine * (target / y_fine) ** (1.0 / p)
                    req[f"dh_km_for_{target:g}_{tag}"] = float(dh)
                    req[f"ns_for_{target:g}_{tag}"] = float(DH_KM / dh)
                else:
                    req[f"dh_km_for_{target:g}_{tag}"] = float("nan")
                    req[f"ns_for_{target:g}_{tag}"] = float("nan")
        conv["_A_M_requirement"] = req

        # Richardson: for a MEAN quantity the informative number is the observed
        # order from three successive values and the extrapolated Delta h -> 0 limit.
        # p_obs = log(|y1-y2| / |y2-y3|) / log(r),  r = mesh refinement ratio (2)
        rich = {}
        for name in ("f_eff_mean", "M_mean", "q_top_like"):
            if name == "q_top_like":
                y = [r["times"][key]["global"]["q_top"] for r in runs]
            else:
                y = [r["times"][key][name] for r in runs]
            if len(y) < 3:
                continue
            d1, d2 = y[0] - y[1], y[1] - y[2]
            r_mesh = runs[0]["ds_km"] / runs[1]["ds_km"]
            entry = {"values": y}
            if d2 != 0 and (d1 / d2) > 0:
                p_obs = float(np.log(abs(d1 / d2)) / np.log(r_mesh))
                lim = float(y[2] - d2 / (r_mesh ** p_obs - 1.0))
                entry["p_observed"] = p_obs
                entry["richardson_limit"] = lim
                entry["rel_err_at_each_ns"] = [float(abs(v - lim) / abs(lim)) if lim else float("nan")
                                               for v in y]
            else:
                entry["p_observed"] = float("nan")
                entry["richardson_limit"] = float("nan")
                entry["rel_err_at_each_ns"] = [float("nan")] * len(y)
            rich[name] = entry
        conv["_richardson"] = rich
        result["convergence"][key] = conv

    text = render(result, runs, ds, args)
    print(text)
    if args.report:
        with open(args.report, "w") as fh:
            fh.write(text)
        print(f"wrote {args.report}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(strip_series(result), fh, indent=2)
        print(f"wrote {args.json}")
    if args.plot_dir:
        make_plots(result, runs, args)
    return 0


def strip_series(result):
    out = json.loads(json.dumps(result, default=float))
    for r in out["runs"]:
        r.pop("series", None)
    return out


def fmt(v, w=14):
    if isinstance(v, float):
        if not np.isfinite(v):
            return f"{'nan':>{w}}"
        if v == 0.0:
            return f"{0.0:>{w}.4f}"
        if abs(v) < 1e-3 or abs(v) >= 1e5:
            return f"{v:>{w}.4e}"
        return f"{v:>{w}.5f}"
    return f"{str(v):>{w}}"


def render(result, runs, ds, args):
    L = []
    win = result["window"]
    L.append("=" * 96)
    L.append("model_column top-region ripple: RESOLUTION CONVERGENCE SCAN")
    L.append(f"physical window: {win[0]:g} - {win[1]:g} km   (fixed height band, "
             "not a fixed cell count)")
    L.append("=" * 96)
    L.append("")
    L.append("--- configuration parity (from the .faceflux headers) ---")
    hdr = f"{'field':<22}" + "".join(fmt(f"ns={r['ns']}") for r in runs)
    L.append(hdr)
    L.append("-" * len(hdr))
    for f in ("ns", "ds_km", "uniform_mesh", "mc3", "beta", "first_cell"):
        L.append(f"{f:<22}" + "".join(fmt(r["meta"].get(f, r.get(f, "-"))) for r in runs))
    L.append(f"{'n_faceflux_records':<22}" + "".join(fmt(r["n_records"]) for r in runs))
    L.append(f"{'cells in window':<22}"
             + "".join(fmt(r["times"][f'{args.times[0]:g}']["n_cells_window"]) for r in runs))

    for tt in args.times:
        key = f"{tt:g}"
        L.append("")
        L.append("=" * 96)
        L.append(f"=== t = {tt:g} s  (nearest captured record) ===")
        L.append("=" * 96)
        hdr = f"{'metric':<26}" + "".join(fmt(f"ns={r['ns']}") for r in runs) \
            + fmt("p (fit)") + fmt("coarse/mid") + fmt("mid/fine")
        L.append(hdr)
        L.append("-" * len(hdr))
        rows = [
            ("t_record [s]", "t"),
            ("dh [km]", "ds_km"),
            (None, None),
            ("mean(rho V)", "M_mean"),
            ("A_M = relstd(rho V)", "A_M"),
            ("Q_M = roughness(rho V)", "Q_M"),
            (None, None),
            ("mean(f_eff)", "f_eff_mean"),
            ("A_eff", "A_eff"),
            ("Q_eff", "Q_eff"),
            (None, None),
            ("A_diff", "A_diff"),
            ("A_central", "A_central"),
            ("corr(dcen,ddiff)", "corr_dcentral_ddiff"),
            ("median a/|V|", "a_over_V_median"),
            (None, None),
            ("A_mismatch", "A_mismatch"),
            ("mismatch rms", "mismatch_rms"),
            ("mismatch absmax", "mismatch_absmax"),
            (None, None),
            ("R_eff rms", "R_eff_rms"),
            ("R_eff absmax", "R_eff_absmax"),
            ("R_eff rms / (|F|/dh)", "R_eff_rms_norm"),
            ("R_total rms", "R_total_rms"),
            ("R_total rms / (|F|/dh)", "R_total_rms_norm"),
            (None, None),
            ("mean(f_total)", "f_total_mean"),
            ("A_total", "A_total"),
            ("mean(f_diff)", "f_diff_mean"),
            (None, None),
            ("n extrema in window", "n_extrema"),
            ("lambda_extrema [km]", "lambda_km"),
            ("lambda_extrema [cells]", "lambda_cells"),
            ("lambda_peak [km]", "lambda_peak_km"),
            ("lambda_peak [cells]", "lambda_peak_cells"),
            (None, None),
            ("split_max_relerr", "split_max_relerr"),
        ]
        conv = result["convergence"][key]
        for name, field in rows:
            if name is None:
                L.append("")
                continue
            vals = [r["times"][key][field] for r in runs]
            line = f"{name:<26}" + "".join(fmt(v) for v in vals)
            if field in conv:
                c = conv[field]
                line += fmt(c["order_p"]) + fmt(c["ratio_coarse_to_mid"]) \
                    + fmt(c["ratio_mid_to_fine"])
            L.append(line)

        L.append("")
        L.append("--- global / boundary state (secondary; regime check only) ---")
        g0 = runs[0]["times"][key]["global"]
        hdr = f"{'quantity':<26}" + "".join(fmt(f"ns={r['ns']}") for r in runs)
        L.append(hdr)
        L.append("-" * len(hdr))
        for gf in [k for k in g0 if k != "ds_km"]:
            L.append(f"{gf:<26}"
                     + "".join(fmt(r["times"][key]["global"].get(gf, float('nan')))
                               for r in runs))

        L.append("")
        L.append("--- MEAN-flux convergence (Richardson, r=2) ---")
        for name, rr in conv["_richardson"].items():
            L.append(f"  {name:<12} values = "
                     + " ".join(f"{v:.4e}" for v in rr["values"])
                     + f"   p_obs = {rr['p_observed']:+.3f}"
                     + f"   dh->0 limit = {rr['richardson_limit']:.4e}")
            L.append(f"  {'':<12} rel. error vs that limit: "
                     + "  ".join(f"ns={r['ns']}: {e * 100:.1f}%"
                                 for r, e in zip(runs, rr["rel_err_at_each_ns"])))

        L.append("")
        L.append("--- resolution needed for a given A_M (extrapolated) ---")
        req = conv["_A_M_requirement"]
        L.append(f"  A_M(finest, dh={runs[-1]['ds_km']:.4f} km) = "
                 f"{runs[-1]['times'][key]['A_M']:.5f}")
        L.append(f"  p_global(3-point fit) = {req['p_global']:+.3f}   "
                 f"p_local(mid->fine)   = {req['p_local_mid_to_fine']:+.3f}")
        for target in (0.05, 0.02):
            L.append(f"  A_M < {target:g}: dh = {req[f'dh_km_for_{target:g}_global']:.4f} km "
                     f"(ns ~ {req[f'ns_for_{target:g}_global']:.0f}) [global fit]   |   "
                     f"dh = {req[f'dh_km_for_{target:g}_local']:.4f} km "
                     f"(ns ~ {req[f'ns_for_{target:g}_local']:.0f}) [local rate]")

        L.append("")
        L.append("--- rho V extrema positions [km] ---")
        for r in runs:
            pos = r["times"][key]["extrema_km"]
            L.append(f"  ns={r['ns']:<5d} ({len(pos)}): "
                     + " ".join(f"{p:.2f}" for p in pos))

        L.append("")
        L.append("--- extrema pinned to physical height? (match vs the FINEST mesh) ---")
        fine = np.array(runs[-1]["times"][key]["extrema_km"])
        for r in runs[:-1]:
            pos = np.array(r["times"][key]["extrema_km"])
            if pos.size == 0 or fine.size == 0:
                continue
            d = np.min(np.abs(pos[:, None] - fine[None, :]), axis=1)
            tol = 1.5 * r["ds_km"]
            L.append(f"  ns={r['ns']:<5d}: {int((d <= tol).sum())}/{pos.size} extrema have a "
                     f"ns={runs[-1]['ns']} counterpart within 1.5*dh={tol:.2f} km; "
                     f"max offset {d.max():.2f} km, median {np.median(d):.2f} km")
            L.append(f"           extra extrema on the fine mesh: "
                     f"{fine.size - pos.size:+d}")
    return "\n".join(L)


def make_plots(result, runs, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    win = result["window"]
    os.makedirs(args.plot_dir, exist_ok=True)
    made = []

    # 1-3: rho V, f_eff, d_diff vs physical height, one figure per time
    for tt in args.times:
        key = f"{tt:g}"
        fig, axes = plt.subplots(4, 1, figsize=(9.5, 12.5), sharex=True)
        for k, r in enumerate(runs):
            s = r["series"][key]
            cell = np.array(s["cell_km"]); face = np.array(s["face_km"])
            lo, hi = win[0] - 12.0, 2153.5
            mc = (cell >= lo) & (cell <= hi)
            mf = (face >= lo) & (face <= hi)
            lab = rf"ns={r['ns']}, $\Delta h$={r['ds_km']:.4f} km"
            norm = r["times"][key]["f_eff_mean"]
            axes[0].plot(cell[mc], np.array(s["rhoV"])[mc], ".-", ms=3, lw=1.0,
                         color=f"C{k}", label=lab)
            axes[1].plot(cell[mc], np.array(s["rhoV"])[mc] / norm, ".-", ms=3, lw=1.0,
                         color=f"C{k}",
                         label=lab + rf",  $A_M$={r['times'][key]['A_M']:.3f}")
            axes[2].plot(face[mf], np.array(s["f_eff"])[mf], ".-", ms=3, lw=1.0,
                         color=f"C{k}", label=lab)
            axes[3].plot(face[mf], np.array(s["d_diff"])[mf], ".-", ms=3, lw=1.0,
                         color=f"C{k}", label=lab)
        axes[0].set_ylabel(r"cell $\rho V$  [kg m$^{-2}$ s$^{-1}$]")
        axes[1].set_ylabel(r"$\rho V\,/\,\langle F_{\rm eff}\rangle_{\rm window}$")
        axes[1].axhline(1.0, color="0.4", lw=0.7, ls="--", zorder=1)
        axes[2].set_ylabel(r"$F_{\rm eff}=F_{\rm tot}-F_{\rm ref}$")
        axes[3].set_ylabel(r"$\delta F_{\rm diff}$")
        axes[3].set_xlabel("height [km]")
        for ax in axes:
            ax.axvspan(win[0], win[1], color="0.85", zorder=0)
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
        axes[0].set_title(f"Resolution scan at t = {tt:g} s "
                          f"(shaded = analysis window {win[0]:g}-{win[1]:g} km)")
        p = os.path.join(args.plot_dir, f"{args.plot_prefix}_profiles_{tt:g}s.png")
        fig.tight_layout(); fig.savefig(p, dpi=150); plt.close(fig)
        made.append(p)

    # 4: log-log convergence of the core amplitudes
    ds = np.array([r["ds_km"] for r in runs])
    fig, axes = plt.subplots(1, len(args.times), figsize=(6.0 * len(args.times), 5.0),
                             squeeze=False)
    for j, tt in enumerate(args.times):
        ax = axes[0][j]
        conv = result["convergence"][f"{tt:g}"]
        # A_M and A_mismatch nearly coincide (that IS the result), so draw A_M as a
        # large open marker underneath so both stay visible.
        for name, mk, kw in (("A_M", "o", dict(ms=14, mfc="none", mew=1.6, lw=2.5)),
                             ("Q_M", "s", {}), ("A_diff", "^", {}),
                             ("A_mismatch", "D", dict(ms=6)), ("A_eff", "v", {})):
            y = np.abs(np.array(conv[name]["values"], dtype=float))
            ax.loglog(ds, y, mk + "-",
                      label=f"{name}  (p={conv[name]['order_p']:+.2f})", **kw)
        ref = np.array(conv["A_M"]["values"], dtype=float)[0]
        ax.loglog(ds, ref * (ds / ds[0]), "k--", lw=0.8, alpha=0.6, label=r"$\propto\Delta h$")
        ax.loglog(ds, ref * (ds / ds[0]) ** 2, "k:", lw=0.8, alpha=0.6,
                  label=r"$\propto\Delta h^2$")
        ax.set_xlabel(r"$\Delta h$ [km]"); ax.set_ylabel("amplitude metric")
        ax.set_title(f"t = {tt:g} s"); ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
    p = os.path.join(args.plot_dir, f"{args.plot_prefix}_convergence_loglog.png")
    fig.tight_layout(); fig.savefig(p, dpi=150); plt.close(fig)
    made.append(p)

    # 5: mean fluxes vs resolution
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6))
    for j, tt in enumerate(args.times):
        conv = result["convergence"][f"{tt:g}"]
        axes[0].plot(ds, conv["f_eff_mean"]["values"], "o-", label=f"$F_{{eff}}$ t={tt:g}s")
        axes[0].plot(ds, conv["M_mean"]["values"], "s--", label=rf"$\rho V$ t={tt:g}s")
        axes[1].plot(ds, np.abs(np.array(conv["R_eff_rms_norm"]["values"])), "^-",
                     label=f"$R_{{eff}}$ rms norm t={tt:g}s")
    axes[0].set_xlabel(r"$\Delta h$ [km]"); axes[0].set_ylabel("window mean mass flux")
    axes[0].set_title("mean evaporation mass flux vs resolution")
    axes[1].set_xlabel(r"$\Delta h$ [km]"); axes[1].set_yscale("log")
    axes[1].set_ylabel(r"$R_{\rm eff}$ rms / $(|F|/\Delta h)$")
    axes[1].set_title("effective continuity residual")
    for ax in axes:
        ax.grid(alpha=0.3); ax.legend(fontsize=8)
    p = os.path.join(args.plot_dir, f"{args.plot_prefix}_mean_flux.png")
    fig.tight_layout(); fig.savefig(p, dpi=150); plt.close(fig)
    made.append(p)

    # 6: 500 vs 1000 s convergence comparison
    if len(args.times) >= 2:
        fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.6))
        for k, name in enumerate(("A_M", "A_diff", "A_mismatch")):
            for tt in args.times:
                conv = result["convergence"][f"{tt:g}"]
                axes[k].loglog(ds, np.abs(np.array(conv[name]["values"], float)),
                               "o-", label=f"t={tt:g} s (p={conv[name]['order_p']:+.2f})")
            axes[k].set_xlabel(r"$\Delta h$ [km]"); axes[k].set_ylabel(name)
            axes[k].set_title(name); axes[k].grid(alpha=0.3, which="both")
            axes[k].legend(fontsize=8)
        p = os.path.join(args.plot_dir, f"{args.plot_prefix}_time_compare.png")
        fig.tight_layout(); fig.savefig(p, dpi=150); plt.close(fig)
        made.append(p)

    for p in made:
        print(f"wrote {p}")


if __name__ == "__main__":
    sys.exit(main())
