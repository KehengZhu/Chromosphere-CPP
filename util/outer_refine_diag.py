#!/usr/bin/env python3
"""Does OUTER (upper-TR) static local refinement reduce the top-region rho*V ripple?

Compares three `model_column` runs over the identical 1600-2153 km domain and the
identical coarse-equivalent resolution (ISO_NS = 1000):

    U   uniform
    R2  ISO_REFINE_PROFILE=outer, factor 2  (fine band 500-553 km above the base)
    R4  ISO_REFINE_PROFILE=outer, factor 4

Every metric definition is IMPORTED, not re-typed:

  * ripple / flux metrics (A_M, Q_M, mean(F_eff), A_eff, Q_eff, A_diff, A_mismatch)
    come from `resolution_scan_diag.face_metrics` over the fixed PHYSICAL window
    (default 2130-2150 km) of the `.faceflux` sidecar -- a read-only copy-out of the
    production `rhs_explicit_mixture` arrays;
  * the outer conductive-flux split (q_phys / q_num / q_total) comes from
    `numerical_conduction_diag.read_outercond` (the `.outercond` sidecar).

NON-UNIFORM CAVEAT. `face_metrics` takes one scalar `ds` and uses it only to
normalise the continuity residuals R_eff / R_total and the periodogram wavelength.
On a refined mesh we pass the LOCAL cell width inside the analysis window, measured
from the sidecar's own `face_km` column (median face spacing in the window), which
is exactly Delta h on a uniform mesh. A_M, Q_M, A_eff, Q_eff, A_diff, A_mismatch and
the means are ds-independent and therefore unaffected.

Cell widths for the profile panel are taken from the state file's cell-centre row
(w_i ~ c_{i+1} - c_i, i.e. the average of two adjacent widths inside a graded band);
the widths inside the analysis window come from the exact `face_km` faces.

Usage:
    python util/outer_refine_diag.py \
        --run U=outputs/model_column/outref_U_500s.txt \
        --run R2=... --run R4=... \
        --time 500 --window 2130 2150 \
        --json out.json --report out.txt --plot visualization/model_column/x.png
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from face_flux_diag import read_faceflux  # noqa: E402
from numerical_conduction_diag import pick_outer, read_outercond  # noqa: E402
from resolution_scan_diag import face_metrics, pick_record  # noqa: E402
from upper_bc_decouple_diag import load_gamma_diag, pick as pick_gamma  # noqa: E402


def read_cell_heights(path):
    """Cell-centre heights [km] from the state file's second row."""
    with open(path) as fh:
        fh.readline()                      # "<ns> <ncol>"
        return np.array([float(v) for v in fh.readline().split()])


def mesh_stats(state_path, faces_km, win):
    """Mesh geometry.

    Precise widths come from the `.faceflux` `face_km` column (full float32 precision,
    but only over the captured top region). The state file's cell-centre row is printed
    at 6 significant digits, i.e. 10 m granularity, so it is used ONLY for the
    full-domain width panel of the figure and never for a quoted width or ratio.
    """
    h = read_cell_heights(state_path)
    wf = np.diff(faces_km)                 # exact widths, capture region only
    ratio = np.maximum(wf[1:] / wf[:-1], wf[:-1] / wf[1:])
    inwin = (faces_km[:-1] >= win[0]) & (faces_km[:-1] <= win[1])
    return {
        "ncells": int(len(h)),
        "capture_cells": int(len(wf)),
        "capture_width_min_km": float(wf.min()),
        "capture_width_max_km": float(wf.max()),
        "capture_max_adjacent_ratio": float(ratio.max()),
        "window_width_km": float(np.median(wf[inwin])),
        "top_width_km": float(wf[-1]),
        "cell_km": h.tolist(), "width_km": np.diff(h).tolist(),
    }


def outer_boundary_layer(h_km, T, frac=0.1):
    """Outer thermal boundary layer from the `.gamma_diag` temperature profile.

    ONE definition, applied identically at every resolution (do not switch it
    between meshes):

      G(h)   = |dT/ds| at cell centres (np.gradient on the ACTUAL non-uniform
               cell-centre heights, so a graded mesh is handled correctly);
      G_max  = max G over the whole column;
      edge   = the deepest cell of the CONTIGUOUS run of top cells with
               G >= frac * G_max   (frac = 0.1, i.e. the conventional
               10 %-of-peak-gradient layer edge);
      delta  = (outer cell centre) - (edge cell centre), in km.

    `n_cells` is the number of cells in that contiguous run, i.e. how many cells
    actually resolve the layer. `h_Gmax` says whether the steepest gradient sits
    at the physical outer boundary or inside the resolved domain.
    """
    s = np.asarray(h_km) * 1.0e3
    G = np.abs(np.gradient(np.asarray(T), s))
    gmax = float(G.max())
    thr = frac * gmax
    i = len(G) - 1
    while i > 0 and G[i - 1] >= thr:
        i -= 1
    return {
        "bl_thickness_km": float(h_km[-1] - h_km[i]),
        "bl_n_cells": int(len(G) - i),
        "bl_edge_km": float(h_km[i]),
        "bl_T_edge": float(T[i]),
        "G_wall_profile_max": gmax,
        "h_Gmax_km": float(h_km[int(np.argmax(G))]),
        "G_top_cell": float(G[-1]),
        "Gmax_at_outer_cell": int(int(np.argmax(G)) == len(G) - 1),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, metavar="LABEL=PATH")
    ap.add_argument("--time", type=float, default=500.0)
    ap.add_argument("--window", type=float, nargs=2, default=[2130.0, 2150.0],
                    metavar=("LO_KM", "HI_KM"))
    ap.add_argument("--json")
    ap.add_argument("--report")
    ap.add_argument("--plot")
    args = ap.parse_args(argv)

    win = tuple(args.window)
    out, profiles = {}, {}
    for spec in args.run:
        label, path = spec.split("=", 1)
        meta, recs = read_faceflux(path + ".faceflux")
        ref, rec = recs[0][2], pick_record(recs, args.time)
        col = rec[2]
        ms = mesh_stats(path, col["face_km"], win)
        # local window width for the R_eff normalisation (see module docstring)
        fm = face_metrics(rec, ref, ms["window_width_km"], win)
        oc = read_outercond(path + ".outercond")
        ocr = pick_outer(oc, args.time)
        gh, gframes = load_gamma_diag(path)
        _, ga = pick_gamma(gframes, args.time)
        gm = ga[:, 0] * ga[:, 1]
        grade = (gh >= 2079.0) & (gh <= 2101.0)
        gseg = gm[grade]
        gd = np.diff(gseg)
        gd2 = np.diff(gseg, 2)
        out[label] = {
            "path": path,
            "uniform_mesh": int(meta.get("uniform_mesh", 1)),
            "ns_coarse_ds_km": float(meta["ds_km"]),
            "numerical_diffusivity_per_length": float(
                oc["_meta"].get("numerical_diffusivity_per_length", "nan")),
            "legacy_numerical_diffusivity": float(
                oc["_meta"].get("numerical_diffusivity", "nan")),
            "mesh": {k: v for k, v in ms.items() if k not in ("cell_km", "width_km")},
            "t_faceflux": fm["t"], "step": fm["step"], "t_outercond": ocr["t"],
            "n_cells_window": fm["n_cells_window"],
            "A_M": fm["A_M"], "Q_M": fm["Q_M"],
            "A_eff": fm["A_eff"], "Q_eff": fm["Q_eff"],
            "A_diff": fm["A_diff"], "A_mismatch": fm["A_mismatch"],
            "corr_dcentral_ddiff": fm["corr_dcentral_ddiff"],
            "M_mean": fm["M_mean"], "f_eff_mean": fm["f_eff_mean"],
            "n_extrema": fm["n_extrema"], "lambda_km": fm["lambda_km"],
            "T_top": ocr["T_top"], "q_phys": ocr["q_phys"], "q_num": ocr["q_num"],
            "q_total": ocr["q_total"],
            "chi_num_face": ocr["chi_num_face"],
            "kappa_phys_face": ocr["kappa_phys_face"],
            "kappa_num_face": ocr["kappa_num_face"],
            "q_num_frac": ocr["q_num"] / ocr["q_total"] if ocr["q_total"] else float("nan"),
            "split_max_relerr": fm["split_max_relerr"],
            "grade_monotone": int(np.all(gd >= 0.0) or np.all(gd <= 0.0)),
            "grade_Q_M": float(np.mean(np.abs(gd2)) / np.mean(np.abs(gseg))),
            "ds_face_outer_km": float(ocr["ds_face"]) / 1.0e3,
            "T_wall": ocr["T_wall"],
            "dT_wall": ocr["T_wall"] - ocr["T_top"],
            "G_wall": (ocr["T_wall"] - ocr["T_top"]) / ocr["ds_face"],
        }
        out[label].update(outer_boundary_layer(gh, ga[:, 2]))
        profiles[label] = {
            "cell_km": ms["cell_km"], "width_km": ms["width_km"],
            "win_cell_km": col["cell_km"], "rhoV": (col["rho_cell"] * col["v_cell"]),
        }

    lines = ["Outer (upper-TR) static local refinement, t = %g s, window %g-%g km"
             % (args.time, win[0], win[1]), ""]
    labels = list(out)
    rows = [
        ("cells (total)", "mesh.ncells", "%d"),
        ("capture cells", "mesh.capture_cells", "%d"),
        ("min width, capture [km]", "mesh.capture_width_min_km", "%.4f"),
        ("max width, capture [km]", "mesh.capture_width_max_km", "%.4f"),
        ("max ratio, capture", "mesh.capture_max_adjacent_ratio", "%.4f"),
        ("top cell width [km]", "mesh.top_width_km", "%.4f"),
        ("window cell width [km]", "mesh.window_width_km", "%.4f"),
        ("cells in window", "n_cells_window", "%d"),
        ("uniform_mesh flag", "uniform_mesh", "%d"),
        ("C_num [m/s]", "numerical_diffusivity_per_length", "%.4g"),
        ("chi_num outer [m2/s]", "chi_num_face", "%.4g"),
        ("kappa_phys outer", "kappa_phys_face", "%.4g"),
        ("kappa_num outer", "kappa_num_face", "%.4g"),
        ("t (faceflux) [s]", "t_faceflux", "%.3f"),
        ("steps", "step", "%d"),
        ("A_M  (ripple)", "A_M", "%.5f"),
        ("Q_M  (roughness)", "Q_M", "%.5f"),
        ("A_mismatch", "A_mismatch", "%.5f"),
        ("A_diff", "A_diff", "%.5f"),
        ("A_eff (face flux)", "A_eff", "%.5f"),
        ("Q_eff (face flux)", "Q_eff", "%.4e"),
        ("corr(dcen,ddiff)", "corr_dcentral_ddiff", "%.5f"),
        ("mean(rho V)", "M_mean", "%.5e"),
        ("mean(F_eff)", "f_eff_mean", "%.5e"),
        ("n extrema", "n_extrema", "%d"),
        ("ds_face outer [km]", "ds_face_outer_km", "%.5f"),
        ("T_top [K]", "T_top", "%.1f"),
        ("T_wall - T_top [K]", "dT_wall", "%.1f"),
        ("G_wall [K/m]", "G_wall", "%.5f"),
        ("BL thickness [km]", "bl_thickness_km", "%.3f"),
        ("BL cells", "bl_n_cells", "%d"),
        ("BL edge [km]", "bl_edge_km", "%.3f"),
        ("BL T_edge [K]", "bl_T_edge", "%.1f"),
        ("max|dT/ds| [K/m]", "G_wall_profile_max", "%.4f"),
        ("h(max|dT/ds|) [km]", "h_Gmax_km", "%.3f"),
        ("|dT/ds| top cell [K/m]", "G_top_cell", "%.4f"),
        ("max grad at outer cell", "Gmax_at_outer_cell", "%d"),
        ("q_phys [W/m2]", "q_phys", "%.4f"),
        ("q_num  [W/m2]", "q_num", "%.4f"),
        ("q_total [W/m2]", "q_total", "%.4f"),
        ("q_num/q_total", "q_num_frac", "%.4f"),
        ("split max relerr", "split_max_relerr", "%.2e"),
        ("grade monotone (2079-2101)", "grade_monotone", "%d"),
        ("grade Q_M", "grade_Q_M", "%.4e"),
    ]
    w0 = max(len(r[0]) for r in rows) + 2
    wv = max(14, max(len(l) + 2 for l in labels))
    lines.append("".ljust(w0) + "".join(l.rjust(wv) for l in labels))
    for name, key, fmt in rows:
        vals = []
        for l in labels:
            v = out[l]
            for part in key.split("."):
                v = v[part]
            vals.append((fmt % v).rjust(wv))
        lines.append(name.ljust(w0) + "".join(vals))
    report = "\n".join(lines)
    print(report)
    if args.report:
        with open(args.report, "w") as fh:
            fh.write(report + "\n")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=2, sort_keys=True)

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(9.0, 8.0))
        style = {"U": ("k", "-"), "R2": ("tab:blue", "-"), "R4": ("tab:red", "-"),
                 "R4_mult025": ("tab:green", "--")}
        for l in labels:
            c, ls = style.get(l, ("tab:purple", ":"))
            p = profiles[l]
            cell = np.asarray(p["cell_km"])
            ax[0].step(cell[:-1], np.asarray(p["width_km"]) * 1.0e3, where="post",
                       color=c, ls=ls, lw=1.2, label=l)
            ax[1].plot(p["win_cell_km"], p["rhoV"], color=c, ls=ls, lw=1.0, label=l)
        ax[0].set_xlabel("height [km]"); ax[0].set_ylabel(r"cell width $\Delta s$ [m]")
        ax[0].set_yscale("log"); ax[0].set_xlim(2000.0, 2153.0)
        ax[0].set_title("Static outer refinement: cell width vs height (ISO_NS=1000, "
                        "fine band 2100-2153 km)\nfrom output cell centres "
                        r"($\pm$10 m print granularity)", fontsize=10)
        ax[0].legend(); ax[0].grid(alpha=0.3)
        ax[1].set_xlim(win[0] - 5.0, 2153.0)
        ax[1].axvspan(win[0], win[1], color="0.85", zorder=0)
        ax[1].set_xlabel("height [km]")
        ax[1].set_ylabel(r"$\rho V$ [kg m$^{-2}$ s$^{-1}$]")
        ax[1].set_title(r"cell-centred $\rho V$ at t = %g s (shaded: metric window)"
                        % args.time)
        ax[1].legend(); ax[1].grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=140)
        print("wrote %s" % args.plot)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
