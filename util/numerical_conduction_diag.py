#!/usr/bin/env python3
"""Is `numerical_diffusivity = 2000 * Delta h` the reason the evaporation mass flux
does not converge under mesh refinement?

Reads the 2x2 set of runs (ns = 1000, 2000) x (production numerical conduction,
ISO_NUMERICAL_DIFFUSIVITY_MULT=0) and reports, at each requested time:

  * from the `.outercond` sidecar (read-only capture of the OUTER face of the FINAL
    CONVERGED conduction Newton iteration):
        q_phys  = kappa_phys_face * (T_wall - T_top)/ds_face
        q_num   = kappa_num_face  * (T_wall - T_top)/ds_face,
                  kappa_num_face = numerical_diffusivity * C_V(top)
        q_total = q_phys + q_num                    [W m^-2, into the top cell]
    kappa_phys_face is the arithmetic face average of the top cell's and the wall
    ghost's (kappa_e*TRAC + kappa_n), and ds_face = ds_iph_i(ns-1) = Delta h (top
    cell CENTRE to ghost CENTRE).
  * from the `.faceflux` sidecar, via the SHARED definitions in
    resolution_scan_diag.face_metrics: mean(F_eff) over the physical window, with
    F_eff = F_total - F_total(t=0) the equilibrium-reference-subtracted conserved
    total-mass face flux, plus the ripple amplitude A_M for context.

Nothing here modifies or re-derives a metric: the face-flux readers and metric
definitions are imported from the existing scripts.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from face_flux_diag import COLUMNS as FF_COLUMNS, read_faceflux  # noqa: E402
from resolution_scan_diag import face_metrics, pick_record  # noqa: E402

OUTER_COLS = ["t", "step", "T_top", "T_wall", "kappa_phys_face", "chi_num_face",
              "kappa_num_face", "ds_face", "area_ratio", "q_phys", "q_num", "q_total"]
LEGACY_OUTER_COLS = ["t", "step", "T_top", "T_wall", "kappa_phys_face",
                     "kappa_num_face", "ds_face", "area_ratio", "q_phys", "q_num",
                     "q_total"]


def read_faceflux_tolerant(path):
    """read_faceflux, but tolerant of a sidecar truncated mid-write by a crashed run.

    A run that aborts leaves the final line (and hence the final record) partial.
    Drop incomplete lines and any short record; keep the shared reader and metric
    definitions otherwise untouched.
    """
    import tempfile
    ncol = len(FF_COLUMNS)
    kept, dropped = [], 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                kept.append(line)
            elif len(line.split()) == ncol and line.endswith("\n"):
                kept.append(line)
            else:
                dropped += 1
    with tempfile.NamedTemporaryFile("w", suffix=".faceflux", delete=False) as tf:
        tf.writelines(kept)
        tmp = tf.name
    meta, recs = read_faceflux(tmp)
    os.unlink(tmp)
    if recs:
        n = len(recs[0][2]["f_total"])
        full = [r for r in recs if len(r[2]["f_total"]) == n]
        dropped += len(recs) - len(full)
        recs = full
    if dropped:
        print("note: %s -- dropped %d incomplete line(s)/record(s) (truncated sidecar)"
              % (os.path.basename(path), dropped), file=sys.stderr)
    return meta, recs


def read_outercond(path):
    """Parse a `.outercond` sidecar into a dict of column arrays + the header meta."""
    meta, rows = {}, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                for tok in line[1:].split():
                    if "=" in tok:
                        k, v = tok.split("=", 1)
                        meta[k] = v
                continue
            parts = line.split()
            if len(parts) == len(OUTER_COLS):
                rows.append([float(p) for p in parts])
            elif len(parts) == len(LEGACY_OUTER_COLS):
                old = [float(p) for p in parts]
                rows.append(old[:5] + [float("nan")] + old[5:])
    if not rows:
        raise SystemExit("no data rows in %s" % path)
    arr = np.array(rows)
    out = {name: arr[:, i] for i, name in enumerate(OUTER_COLS)}
    out["_meta"] = meta
    return out


def pick_outer(oc, t_target):
    k = int(np.argmin(np.abs(oc["t"] - t_target)))
    return {name: float(oc[name][k]) for name in OUTER_COLS}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, metavar="LABEL=PATH",
                    help="state-file path; .faceflux/.outercond sidecars are implied")
    ap.add_argument("--times", type=float, nargs="+", default=[500.0, 1000.0])
    ap.add_argument("--window", type=float, nargs=2, default=[2130.0, 2150.0],
                    metavar=("LO_KM", "HI_KM"))
    ap.add_argument("--json")
    ap.add_argument("--report")
    ap.add_argument("--plot")
    args = ap.parse_args(argv)

    runs = []
    for spec in args.run:
        label, path = spec.split("=", 1)
        runs.append((label, path))

    result = {}
    for label, path in runs:
        meta, recs = read_faceflux_tolerant(path + ".faceflux")
        ds_km = float(meta["ds_km"])
        ref = recs[0][2]
        oc = read_outercond(path + ".outercond")
        entry = {"path": path, "ds_km": ds_km,
                 "numerical_diffusivity": float(oc["_meta"]["numerical_diffusivity"]),
                 "ns": int(oc["_meta"]["ns"]), "times": {}}
        for t in args.times:
            rec = pick_record(recs, t)
            fm = face_metrics(rec, ref, ds_km, args.window)
            ocr = pick_outer(oc, t)
            entry["times"]["%g" % t] = {
                "faceflux_t": fm["t"], "outercond_t": ocr["t"],
                "f_eff_mean": fm["f_eff_mean"], "M_mean": fm["M_mean"],
                "A_M": fm["A_M"], "n_faces_window": fm["n_faces_window"],
                "q_phys": ocr["q_phys"], "q_num": ocr["q_num"],
                "q_total": ocr["q_total"],
                "q_num_frac": ocr["q_num"] / ocr["q_total"] if ocr["q_total"] else float("nan"),
                "T_top": ocr["T_top"], "T_wall": ocr["T_wall"],
                "dT": ocr["T_wall"] - ocr["T_top"],
                "dT_over_ds": (ocr["T_wall"] - ocr["T_top"]) / ocr["ds_face"],
                "kappa_phys_face": ocr["kappa_phys_face"],
                "kappa_num_face": ocr["kappa_num_face"],
            }
        result[label] = entry

    lines = []

    def emit(s=""):
        lines.append(s)

    emit("numerical-conduction vs evaporation-flux convergence")
    emit("window = %g-%g km   times = %s s" % (args.window[0], args.window[1],
                                               " ".join("%g" % t for t in args.times)))
    emit("q [W m^-2] at the OUTER face, into the top cell; F_eff [kg m^-2 s^-1]")
    emit()
    hdr = ("%-22s %9s %11s %8s %8s %9s %9s %9s %8s %10s %10s %11s"
           % ("run", "ds_km", "D_num", "t_q", "t_F", "T_top", "dT", "q_phys",
              "q_num", "q_total", "q_num/q_tot", "mean(F_eff)"))
    for t in args.times:
        emit("--- t = %g s " % t + "-" * (len(hdr) - 14))
        emit("(t_q / t_F = the ACTUAL sidecar times used; a large gap means the run "
             "ended early)")
        emit(hdr)
        for label, entry in result.items():
            d = entry["times"]["%g" % t]
            emit("%-22s %9.4f %11.4g %8.1f %8.1f %9.1f %9.1f %9.4f %8.4f %10.4f "
                 "%10.3f %11.4e"
                 % (label, entry["ds_km"], entry["numerical_diffusivity"],
                    d["outercond_t"], d["faceflux_t"],
                    d["T_top"], d["dT"], d["q_phys"], d["q_num"], d["q_total"],
                    d["q_num_frac"], d["f_eff_mean"]))
        emit()

    # 1000 -> 2000 change, per branch. Labels are expected to be <branch>_ns<N>.
    emit("=== ns 1000 -> 2000 change, per branch " + "=" * 30)
    branches = {}
    for label, entry in result.items():
        branches.setdefault(label.rsplit("_ns", 1)[0] if "_ns" in label
                            else label.split("_")[0], {})[entry["ns"]] = entry
    for bname, by_ns in branches.items():
        if 1000 not in by_ns or 2000 not in by_ns:
            emit("%-14s incomplete pair (ns present: %s)"
                 % (bname, sorted(by_ns)))
            continue
        for t in args.times:
            a = by_ns[1000]["times"]["%g" % t]
            b = by_ns[2000]["times"]["%g" % t]
            emit("%-14s t=%-6g q_total %8.4f -> %8.4f (%+7.1f %%)   "
                 "q_phys %8.4f -> %8.4f (%+7.1f %%)   "
                 "mean(F_eff) %10.4e -> %10.4e (%+7.1f %%)"
                 % (bname, t, a["q_total"], b["q_total"],
                    100.0 * (b["q_total"] / a["q_total"] - 1.0),
                    a["q_phys"], b["q_phys"],
                    100.0 * (b["q_phys"] / a["q_phys"] - 1.0),
                    a["f_eff_mean"], b["f_eff_mean"],
                    100.0 * (b["f_eff_mean"] / a["f_eff_mean"] - 1.0)))
        emit()

    text = "\n".join(lines)
    print(text)
    if args.report:
        with open(args.report, "w") as fh:
            fh.write(text + "\n")
        print("wrote %s" % args.report)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(result, fh, indent=2, sort_keys=True)
        print("wrote %s" % args.json)

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
        style = {"prod": ("tab:red", "o", "-", "production $D_{num}=2000\\,\\Delta h$"),
                 "dnum0": ("tab:blue", "s", "--", "$D_{num}=0$")}
        for bname, by_ns in branches.items():
            color, marker, ls, lab = style.get(bname, ("k", "^", ":", bname))
            ns = sorted(by_ns)
            for t, alpha in zip(args.times, (0.45, 1.0)):
                q = [by_ns[n]["times"]["%g" % t]["q_total"] for n in ns]
                f = [by_ns[n]["times"]["%g" % t]["f_eff_mean"] for n in ns]
                suffix = "" if alpha == 1.0 else " (%g s)" % t
                axes[0].plot(ns, q, color=color, marker=marker, ls=ls, alpha=alpha,
                             label=(lab + suffix) if alpha == 1.0 else None)
                axes[1].plot(ns, f, color=color, marker=marker, ls=ls, alpha=alpha,
                             label=(lab + suffix) if alpha == 1.0 else None)
            # numerical share, production branch only
            if bname == "prod":
                qn = [by_ns[n]["times"]["%g" % args.times[-1]]["q_num"] for n in ns]
                axes[0].plot(ns, qn, color=color, marker=marker, ls=":", alpha=0.8,
                             label="$q_{num}$ only")
        axes[0].set_ylabel(r"$q$ at the outer face  [W m$^{-2}$]")
        axes[1].set_ylabel(r"$\langle F_{\rm eff}\rangle_{%g-%g\,\rm km}$"
                           r"  [kg m$^{-2}$ s$^{-1}$]" % tuple(args.window))
        for ax in axes:
            ax.set_xlabel("$N$ (cells over 1600-2153 km)")
            ns_all = sorted({e["ns"] for e in result.values()})
            ax.set_xticks(ns_all)
            ax.set_xlim(ns_all[0] - 0.12*(ns_all[-1]-ns_all[0]),
                        ns_all[-1] + 0.12*(ns_all[-1]-ns_all[0]))
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
        axes[0].set_title("outer conductive flux (faint = %g s, solid = %g s)"
                          % (args.times[0], args.times[-1]), fontsize=9)
        axes[1].set_title("mean evaporation mass flux", fontsize=9)
        # Flag any point whose sidecar time is far from the requested one (crashed run).
        for ax in axes:
            for bname, by_ns in branches.items():
                for n, entry in by_ns.items():
                    d = entry["times"]["%g" % args.times[-1]]
                    if abs(d["outercond_t"] - args.times[-1]) > 20.0:
                        ax.annotate("run ended early\n(%s, N=%d: q at %.0f s)"
                                    % (bname, n, d["outercond_t"]),
                                    xy=(0.03, 0.06), xycoords="axes fraction",
                                    fontsize=7, color="0.35")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print("wrote %s" % args.plot)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
