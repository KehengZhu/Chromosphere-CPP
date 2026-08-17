#!/usr/bin/env python3
"""Finite-volume face mass-flux decomposition for the model_column top region.

Reads the `<run>.faceflux` sidecar written by `chromo_main` with
`CHROMO_FACE_FLUX_DIAG=1`. That sidecar is a read-only copy-out of the PRODUCTION
`mixture_rhs_explicit` arrays, so nothing here re-derives the reconstruction: the
production total-mass face flux, its central/diffusive split, the one-sided
reconstructed face states, the spectral radius and the limiter ratios/values are
all the numbers the solver actually used.

Definitions used below (index i = the UPPER face i+1/2 of cell i):

    rhoV_cell[i] = rho_cell[i] * v_cell[i]                  cell-centred mass flux
    f_central[i] = 0.5 ( (rho v)_L + (rho v)_R )            advective part
    f_diff[i]    = f_total[i] - f_central[i]                dissipative part
                   (= -0.5 a[i] (rho_R - rho_L) only in Rusanov mode; the release
                    Godunov flux and the Roe reference flux are defined as a whole)
    f_total[i]   = f_central[i] + f_diff[i]                 what continuity differences
    R_rho[i]     = -( f_total[i] - f_total[i-1] ) / ds      continuity residual

    f_ref[i]     = f_total[i] at t = 0                      t = 0 baseline flux
    f_eff[i]     = f_total[i] - f_ref[i]                    CHANGE in the flux since t=0
    R_eff[i]     = -( f_eff[i] - f_eff[i-1] ) / ds

`f_eff` is a diagnostic baseline subtraction, not something the solver does: the
release RHS differences `f_total` alone. Removing the (large, nearly hydrostatic)
t = 0 flux is simply the cleanest way to see how the transport has CHANGED since
the initial condition. Both are reported.

Roughness of a series x over a window is `mean|D2 x| / mean|x|` with D2 the second
difference -- the same normalised-curvature definition used in
`util/upper_bc_decouple_diag.py`.
"""

from __future__ import annotations

import argparse
import json
import math
import sys

import numpy as np

COLUMNS = ("cell_km face_km rho_cell v_cell T_cell rho_L rho_R v_L v_R T_L T_R "
           "cs_L cs_R a_face f_central f_diff f_total "
           "r_rho phi_plus_rho r_ip1_rho phi_minus_rho r_v phi_plus_v "
           "r_T phi_plus_T").split()


def read_faceflux(path):
    """-> (meta dict, [(t, step, {column: array}), ...]).

    The column names come from the sidecar's own `# columns=` header line; the
    writer has gained columns over time (p_cell/p_L/p_R), so COLUMNS is only the
    fallback for older files that predate the header. Every name in COLUMNS must
    still be present.
    """
    meta, records = {}, []
    columns = list(COLUMNS)
    cur_t = cur_step = None
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("# columns="):
                    columns = line.split("=", 1)[1].split()
                    missing = [c for c in COLUMNS if c not in columns]
                    if missing:
                        raise SystemExit(f"{path}: sidecar lacks column(s) {missing}")
                elif line.startswith("# t "):
                    if rows:
                        records.append((cur_t, cur_step, rows))
                        rows = []
                    parts = line.split()
                    cur_t, cur_step = float(parts[3]), int(parts[6])
                elif "=" in line:
                    for tok in line[1:].split():
                        if "=" in tok and not tok.startswith("columns"):
                            k, _, v = tok.partition("=")
                            meta.setdefault(k, v)
                continue
            rows.append([float(v) for v in line.split()])
    if rows:
        records.append((cur_t, cur_step, rows))
    out = []
    for t, step, rows in records:
        arr = np.asarray(rows, dtype=float)
        if arr.shape[1] != len(columns):
            raise SystemExit(f"{path}: expected {len(columns)} columns, got {arr.shape[1]}")
        out.append((t, step, {name: arr[:, k] for k, name in enumerate(columns)}))
    return meta, out


def roughness(x):
    x = np.asarray(x, dtype=float)
    if x.size < 3:
        return float("nan")
    return float(np.mean(np.abs(np.diff(x, 2))) / max(np.mean(np.abs(x)), 1e-300))


def relstd(x):
    x = np.asarray(x, dtype=float)
    return float(np.std(x) / max(abs(np.mean(x)), 1e-300))


def window_mask(km, lo, hi):
    return (km >= lo) & (km <= hi)


def ripple_peaks(km, x, n=4):
    """Interior local extrema of x, strongest |curvature| first."""
    x = np.asarray(x, dtype=float)
    idx = []
    for i in range(1, len(x) - 1):
        if (x[i] - x[i - 1]) * (x[i + 1] - x[i]) < 0.0:
            idx.append(i)
    idx.sort(key=lambda i: -abs(x[i + 1] - 2 * x[i] + x[i - 1]))
    return idx[:n]


def analyse(label, meta, record, ref, ds_km, win, top_cells):
    t, step, c = record
    ds = ds_km * 1.0e3
    km = c["cell_km"]
    rhoV = c["rho_cell"] * c["v_cell"]
    f_tot, f_cen, f_dif = c["f_total"], c["f_central"], c["f_diff"]
    f_ref = ref["f_total"]
    f_eff = f_tot - f_ref
    # baseline-corrected split: which HALF of the Rusanov flux carries the
    # structure that survives the t = 0 reference subtraction?
    d_cen = f_cen - ref["f_central"]
    d_dif = f_dif - ref["f_diff"]

    # continuity residuals; index 0 of the window has no lower face inside it, so
    # the residual arrays are defined from index 1 of the captured window onward.
    R_tot = -np.diff(f_tot) / ds
    R_eff = -np.diff(f_eff) / ds
    km_R = km[1:]

    masks = {
        "win_%g_%g" % win: window_mask(km, *win),
        "top%d" % top_cells: np.arange(len(km)) >= len(km) - top_cells,
        "top%d_excl_last" % top_cells:
            (np.arange(len(km)) >= len(km) - top_cells) & (np.arange(len(km)) < len(km) - 1),
    }

    out = {"label": label, "t": t, "step": step, "ds_km": ds_km}
    for name, m in masks.items():
        mR = m[1:]
        scale = max(np.mean(np.abs(f_tot[m])), 1e-300) / ds
        frac = np.abs(f_dif[m]) / np.maximum(np.abs(f_tot[m]), 1e-300)
        out[name] = {
            "n": int(m.sum()),
            "rhoV_mean": float(np.mean(rhoV[m])),
            "rhoV_relstd": relstd(rhoV[m]),
            "rhoV_roughness": roughness(rhoV[m]),
            "f_total_mean": float(np.mean(f_tot[m])),
            "f_total_relstd": relstd(f_tot[m]),
            "f_total_roughness": roughness(f_tot[m]),
            "f_central_mean": float(np.mean(f_cen[m])),
            "f_diff_mean": float(np.mean(f_dif[m])),
            "f_eff_mean": float(np.mean(f_eff[m])),
            "f_eff_relstd": relstd(f_eff[m]),
            "f_eff_roughness": roughness(f_eff[m]),
            "d_central_relstd_of_feff": float(np.std(d_cen[m])
                                              / max(abs(np.mean(f_eff[m])), 1e-300)),
            "d_diff_relstd_of_feff": float(np.std(d_dif[m])
                                           / max(abs(np.mean(f_eff[m])), 1e-300)),
            "corr_dcentral_ddiff": (float(np.corrcoef(d_cen[m], d_dif[m])[0, 1])
                                    if d_cen[m].std() > 0 and d_dif[m].std() > 0
                                    else float("nan")),
            "diff_fraction_mean": float(np.mean(frac)),
            "diff_fraction_median": float(np.median(frac)),
            "diff_fraction_max": float(np.max(frac)),
            "R_total_absmax": float(np.max(np.abs(R_tot[mR]))),
            "R_total_rms": float(np.sqrt(np.mean(R_tot[mR] ** 2))),
            "R_total_rms_over_scale": float(np.sqrt(np.mean(R_tot[mR] ** 2)) / scale),
            "R_eff_absmax": float(np.max(np.abs(R_eff[mR]))),
            "R_eff_rms": float(np.sqrt(np.mean(R_eff[mR] ** 2))),
            "R_eff_rms_over_scale": float(np.sqrt(np.mean(R_eff[mR] ** 2)) / scale),
        }
        # (10) is the cell-centred ripple anti-correlated with the diffusive flux?
        face_rhoV = 0.5 * (rhoV + np.roll(rhoV, -1))
        face_rhoV[-1] = rhoV[-1]
        a, b = face_rhoV[m], f_dif[m]
        if a.std() > 0 and b.std() > 0:
            out[name]["corr_rhoV_fdiff"] = float(np.corrcoef(a, b)[0, 1])
            out[name]["corr_fcentral_fdiff"] = float(np.corrcoef(f_cen[m], b)[0, 1])
        else:
            out[name]["corr_rhoV_fdiff"] = float("nan")
            out[name]["corr_fcentral_fdiff"] = float("nan")

    # (9) full state at the strongest ripple peaks of the cell-centred mass flux
    m = masks["win_%g_%g" % win]
    sel = np.where(m)[0]
    peaks = [sel[i] for i in ripple_peaks(km[m], rhoV[m])]
    out["peaks"] = []
    for i in peaks:
        out["peaks"].append({k: float(c[k][i]) for k in COLUMNS} | {
            "rhoV_cell": float(rhoV[i]),
            "f_eff": float(f_eff[i]),
            "f_ref": float(f_ref[i]),
            "branch_rho_L": mc3_branch(c["r_rho"][i], float(meta.get("beta", 2.0)), True),
            "branch_rho_R": mc3_branch(c["r_ip1_rho"][i], float(meta.get("beta", 2.0)), False),
        })

    # sanity: the captured split must reproduce the captured total
    out["split_max_relerr"] = float(np.max(
        np.abs(f_cen + f_dif - f_tot) / np.maximum(np.abs(f_tot), 1e-300)))
    return out


def mc3_branch(r, beta, plus):
    """Which arm of min(beta r, beta, third_order), floored at 0, is active."""
    if not math.isfinite(r):
        return "nonfinite->0"
    third = (2.0 * r + 1.0) / 3.0 if plus else (r + 2.0) / 3.0
    cands = {"beta_r": beta * r, "beta_cap": beta, "third_order": third}
    name = min(cands, key=cands.get)
    return "zero_floor" if cands[name] <= 0.0 else name


def fmt(v, w=12):
    if isinstance(v, float):
        if v == 0.0:
            return f"{0.0:>{w}.4f}"
        if abs(v) < 1e-3 or abs(v) >= 1e5:
            return f"{v:>{w}.4e}"
        return f"{v:>{w}.5f}"
    return f"{str(v):>{w}}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, metavar="LABEL=PATH",
                    help="a <run>.faceflux sidecar (repeatable)")
    ap.add_argument("--window", nargs=2, type=float, default=(2130.0, 2150.0),
                    metavar=("LO_KM", "HI_KM"))
    ap.add_argument("--top-cells", type=int, default=90)
    ap.add_argument("--json", default=None)
    ap.add_argument("--plot", default=None)
    args = ap.parse_args(argv)

    results = []
    series = {}
    for spec in args.run:
        label, _, path = spec.partition("=")
        meta, records = read_faceflux(path)
        if len(records) < 2:
            raise SystemExit(f"{path}: need the t=0 reference record plus a final record")
        ref = records[0][2]
        ds_km = float(meta.get("ds_km", 0.553))
        res = analyse(label, meta, records[-1], ref, ds_km,
                      tuple(args.window), args.top_cells)
        res["path"] = path
        res["meta"] = meta
        results.append(res)
        c = records[-1][2]
        series[label] = {
            "cell_km": c["cell_km"], "face_km": c["face_km"],
            "rhoV": c["rho_cell"] * c["v_cell"],
            "f_total": c["f_total"], "f_central": c["f_central"],
            "f_diff": c["f_diff"], "f_eff": c["f_total"] - ref["f_total"],
        }

    key = "win_%g_%g" % tuple(args.window)
    for name in (key, "top%d" % args.top_cells,
                 "top%d_excl_last" % args.top_cells):
        print(f"\n=== window: {name} ===")
        rows = ["rhoV_mean", "rhoV_relstd", "rhoV_roughness",
                "f_central_mean", "f_diff_mean",
                "f_total_mean", "f_total_relstd", "f_total_roughness",
                "f_eff_mean", "f_eff_relstd", "f_eff_roughness",
                "d_central_relstd_of_feff", "d_diff_relstd_of_feff",
                "corr_dcentral_ddiff",
                "diff_fraction_mean", "diff_fraction_median", "diff_fraction_max",
                "R_total_absmax", "R_total_rms", "R_total_rms_over_scale",
                "R_eff_absmax", "R_eff_rms", "R_eff_rms_over_scale",
                "corr_rhoV_fdiff", "corr_fcentral_fdiff"]
        hdr = f"{'metric':<24}" + "".join(fmt(r['label'], 14) for r in results)
        print(hdr)
        print("-" * len(hdr))
        for row in rows:
            print(f"{row:<24}" + "".join(fmt(r[name][row], 14) for r in results))

    print("\n=== sanity ===")
    for r in results:
        print(f"{r['label']:<12} t={r['t']:<10.4f} step={r['step']:<7d} "
              f"split_max_relerr={r['split_max_relerr']:.3e}")

    print(f"\n=== ripple peaks in {args.window[0]}-{args.window[1]} km ===")
    for r in results:
        print(f"\n-- {r['label']} --")
        for p in r["peaks"]:
            print(f"  cell {p['cell_km']:.2f} km / face {p['face_km']:.2f} km")
            print(f"    rho_L={p['rho_L']:.6e}  rho_R={p['rho_R']:.6e}  "
                  f"drho/rho={(p['rho_R']-p['rho_L'])/p['rho_L']:+.3e}")
            print(f"    V_L  ={p['v_L']:.4f}      V_R  ={p['v_R']:.4f}")
            print(f"    T_L  ={p['T_L']:.3f}      T_R  ={p['T_R']:.3f}")
            print(f"    cs_L ={p['cs_L']:.3f}      cs_R ={p['cs_R']:.3f}  a={p['a_face']:.3f}")
            print(f"    f_central={p['f_central']:+.6e}  f_diff={p['f_diff']:+.6e}  "
                  f"f_total={p['f_total']:+.6e}")
            print(f"    f_ref    ={p['f_ref']:+.6e}  f_eff ={p['f_eff']:+.6e}  "
                  f"rhoV_cell={p['rhoV_cell']:+.6e}")
            print(f"    limiter rho: r={p['r_rho']:+.5f} phi+={p['phi_plus_rho']:.5f} "
                  f"[{p['branch_rho_L']}]  r_ip1={p['r_ip1_rho']:+.5f} "
                  f"phi-={p['phi_minus_rho']:.5f} [{p['branch_rho_R']}]")
            print(f"    limiter V:   r={p['r_v']:+.5f} phi+={p['phi_plus_v']:.5f}   "
                  f"T: r={p['r_T']:+.5f} phi+={p['phi_plus_T']:.5f}")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"window": list(args.window), "top_cells": args.top_cells,
                       "runs": results}, fh, indent=2)
        print(f"\nwrote {args.json}")

    if args.plot:
        make_plot(series, args.plot, tuple(args.window), args.top_cells)
        print(f"wrote {args.plot}")
    return 0


def make_plot(series, path, win, top_cells):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(4, 1, figsize=(9.5, 11.0), sharex=True)
    colors = {}
    for k, (label, s) in enumerate(series.items()):
        colors[label] = f"C{k}"
        n = len(s["cell_km"])
        sl = slice(max(0, n - top_cells), n)
        axes[0].plot(s["cell_km"][sl], s["rhoV"][sl], ".-", ms=3,
                     color=colors[label], label=label)
        axes[1].plot(s["face_km"][sl], s["f_total"][sl], ".-", ms=3,
                     color=colors[label], label=f"{label} total")
        axes[1].plot(s["face_km"][sl], s["f_central"][sl], "--", lw=0.9,
                     color=colors[label], label=f"{label} central")
        axes[2].plot(s["face_km"][sl], s["f_diff"][sl], ".-", ms=3,
                     color=colors[label], label=label)
        axes[3].plot(s["face_km"][sl], s["f_eff"][sl], ".-", ms=3,
                     color=colors[label], label=label)

    axes[0].set_ylabel(r"cell $\rho V$  [kg m$^{-2}$ s$^{-1}$]")
    axes[1].set_ylabel(r"face $F^\rho$: total / central")
    axes[2].set_ylabel(r"face $F^\rho_{\rm diff}$")
    axes[3].set_ylabel(r"$F^\rho_{\rm total}-F^\rho_{\rm ref}$  (change since $t=0$)")
    axes[3].set_xlabel("height [km]")
    for ax in axes:
        ax.axvspan(win[0], win[1], color="0.85", zorder=0)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, ncol=2)
    axes[0].set_title("Production total-mass face flux vs cell-centred "
                      r"$\rho V$ (top %d cells)" % top_cells)
    fig.tight_layout()
    fig.savefig(path, dpi=150)


if __name__ == "__main__":
    sys.exit(main())
