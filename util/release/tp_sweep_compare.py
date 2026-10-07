#!/usr/bin/env python3
"""Simulated against analytic evaporation mass flux over the (T_wall, p_top) sweep.

Post-processing of ``util/run/run_tp_sweep.sh``. Every run goes through
``twall_relaxation_diag.analyse`` (wave phase, TR relaxation, h_c, J, L_TR), and
the flux is read in a window after the run is ready,
t_ready = max(wave phase within 5 %, TR energy ratio within 0.05 of its end value).

The release has no analytic L_TR, so every analytic flux borrows one simulated
quantity. Two closures are compared with the simulated top-cell flux F_sim:

  energy : F_E = q_wall / (h_top - h_c)     borrows q_wall; L-free, marker-free
  TR     : F_J = J(T_top, p_top, h_c) / L   borrows L_TR (7.5 kK marker)

each with h_c measured in the run ("meas") and with h_c = h_Saha(6.7 kK, p_top),
the nominal C7 plateau a lookup table would use ("nom"). F_sim/F_E(meas) is the
diag's energy ratio ER and F_sim/F_J(meas) is its F L/J. The same quantities over
the last window of the run are stored with an "_end" suffix: the flux keeps
drifting after t_ready, and a TR that relaxes on a timescale comparable to the run
can pass the t_ready test while still relaxing.

The per-run time series (t, F, L_TR, ER, F L/J, J, q_wall and the smoothed
F_base/F_top - 1 on the face captures), and the time mean and standard deviation
of the face flux F(z) from 1000 s to the end, go to <dir>/tp_sweep_series.npz, keyed
"T<kK>_p<Pa>_N<ns>/<name>", so the figures of util/release/plot_tp_sweep_vs_analytic.py
can be redrawn without the raw runs.

Usage:
    python util/release/tp_sweep_compare.py [--dir outputs/model_column/tp_sweep]
        [--window 300] [--out-dir visualization/model_column] [--prefix tp_sweep]
        [--workers 6]
"""
import argparse
import glob
import json
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analytic_evaporation_table import enthalpy
from animate_isentropic import load_face_flux
from twall_relaxation_diag import analyse, j_theory
from twall_sweep_diag import COLUMNS, load_gamma_diag

ROOT = os.path.dirname(os.path.dirname(HERE))
I_T, I_P = COLUMNS.index("T"), COLUMNS.index("p")
T_C_NOM = 6.7e3
PROFILE_T0 = 1000.0            # start of the face-flux profile statistics [s]
SERIES = ("t", "F", "L", "ER", "FL_J", "J", "qw", "ft", "ratio_s")
NAME = re.compile(r"tp_T(\d+)kK_p([0-9.eE+-]+)_N(\d+)\.txt$")


def one(path, window):
    m = NAME.search(os.path.basename(path))
    T_w, p_set, ns = float(m.group(1)) * 1e3, float(m.group(2)), int(m.group(3))
    try:
        r = analyse((T_w, path))
    except Exception as exc:                      # a failed run is a result, not a crash
        return dict(T_wall_K=T_w, p_set_Pa=p_set, N=ns, excluded=f"analyse: {exc}")
    if "excluded" in r:
        return dict(T_wall_K=T_w, p_set_Pa=p_set, N=ns, excluded=r["excluded"])
    s, met = r["series"], r["metrics"]
    _, tg, d = load_gamma_diag(path)
    T_top, p_top = d[:, -1, I_T], d[:, -1, I_P]
    t = s["t"]
    t_wave = met["settle_5pct_s"]
    t_tr = met["t_TR_within_5pct_s"]
    t_ready = np.nanmax([x if x is not None else np.nan for x in (t_wave, t_tr)])
    t_end = float(t[-1])
    ready = bool(np.isfinite(t_ready) and t_ready + window <= t_end)
    t0 = min(t_ready, t_end - window) if np.isfinite(t_ready) else t_end - window

    def read(t0):
        """Window means over [t0, t0 + window] and the closures built from them."""
        w = (t >= t0) & (t <= t0 + window) & np.isfinite(s["J"])
        if not w.any():
            return None
        mean = lambda y: float(np.nanmean(np.asarray(y)[w]))
        Tt = float(np.interp(t0 + 0.5 * window, tg, T_top))
        pt = float(np.interp(t0 + 0.5 * window, tg, p_top))
        h_top = float(enthalpy(Tt, pt)[0])
        hc_nom = float(enthalpy(T_C_NOM, pt)[0])
        F, qw, L = mean(s["F"]), mean(s["qw"]), mean(s["L"])
        ER, FLJ = mean(s["ER"]), mean(s["FL_J"])
        J_nom = j_theory(Tt, pt, hc_nom)
        return dict(
            T_top_K=Tt, p_top_Pa=pt, F_sim=F, q_wall=qw, L_TR_km=L / 1e3, h_top=h_top,
            h_c_meas=h_top - qw * ER / F if F > 0 else np.nan,   # ER = F (h_top - h_c)/q_wall
            h_c_nom=hc_nom,
            F_E_meas=F / ER, F_E_nom=qw / (h_top - hc_nom),
            F_J_meas=F / FLJ, F_J_nom=J_nom / L,
            ratio_E_meas=ER, ratio_E_nom=F * (h_top - hc_nom) / qw,
            ratio_J_meas=FLJ, ratio_J_nom=F * L / J_nom)

    at_ready, at_end = read(t0), read(t_end - window)
    if at_ready is None:
        return dict(T_wall_K=T_w, p_set_Pa=p_set, N=ns, excluded="no h_c in window")
    # Face-flux profile statistics after the wave phase: where in height F is unsteady.
    fh, ft, ff = load_face_flux(path + ".faceflux")
    ff = np.asarray(ff)[np.asarray(ft) >= PROFILE_T0]
    prof = dict(face_km=np.asarray(fh), F_mean=ff.mean(0), F_std=ff.std(0))
    return dict(
        T_wall_K=T_w, p_set_Pa=p_set, N=ns,
        t_wave_s=t_wave, t_TR_s=t_tr, t_ready_s=float(t_ready) if np.isfinite(t_ready) else None,
        read_window_s=[float(t0), float(t0 + window)], ready=ready, t_end_s=t_end,
        **at_ready, **{f"{k}_end": v for k, v in (at_end or {}).items()},
        F_drift_per_500s=float(-met["decay_rate_measured_per_s"] * 500.0)
        if np.isfinite(met["decay_rate_measured_per_s"]) else None,
        godunov_fallbacks=r["godunov_fallbacks"],
        series={**{k: np.asarray(s[k]) for k in SERIES}, **prof})


def grid_of(rows, key, Ts, Ps):
    Z = np.full((len(Ts), len(Ps)), np.nan)
    for r in rows:
        if r.get("excluded") or r["N"] != 500:
            continue
        Z[Ts.index(r["T_wall_K"]), Ps.index(r["p_set_Pa"])] = r[key]
    return Z


def panel(ax, Z, Ts, Ps, title, norm, cmap, fmt, ready=None):
    xe = np.sqrt(np.array(Ps[:-1]) * np.array(Ps[1:]))
    xe = np.concatenate([[Ps[0] ** 2 / xe[0]], xe, [Ps[-1] ** 2 / xe[-1]]])
    ye = np.sqrt(np.array(Ts[:-1]) * np.array(Ts[1:]))
    ye = np.concatenate([[Ts[0] ** 2 / ye[0]], ye, [Ts[-1] ** 2 / ye[-1]]]) / 1e3
    pc = ax.pcolormesh(xe, ye, np.ma.masked_invalid(Z), norm=norm, cmap=cmap)
    for i, T in enumerate(Ts):
        for j, p in enumerate(Ps):
            if np.isfinite(Z[i, j]):
                flag = "" if ready is None or ready[i, j] else "*"
                ax.text(p, T / 1e3, format(Z[i, j], fmt) + flag, ha="center", va="center",
                        fontsize=6.5)
    ax.set(xscale="log", yscale="log", xlabel=r"$p_{top}$ [Pa]", ylabel=r"$T_{wall}$ [kK]",
           title=title)
    ax.set_yticks([T / 1e3 for T in Ts], [f"{T/1e3:.0f}" for T in Ts])
    ax.set_xticks(Ps, [f"{p:g}" for p in Ps], fontsize=7)
    ax.minorticks_off()
    plt.colorbar(pc, ax=ax, fraction=0.046)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", default=os.path.join(ROOT, "outputs/model_column/tp_sweep"))
    ap.add_argument("--window", type=float, default=300.0)
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "visualization/model_column"))
    ap.add_argument("--prefix", default="tp_sweep")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    paths = sorted(glob.glob(os.path.join(a.dir, "tp_T*kK_p*_N*.txt")))
    with ProcessPoolExecutor(a.workers) as ex:
        rows = list(ex.map(one, paths, [a.window] * len(paths)))
    rows.sort(key=lambda r: (r["N"], r["T_wall_K"], r["p_set_Pa"]))
    series = {}
    for r in rows:
        tag = f"T{r['T_wall_K']/1e3:.0f}_p{r['p_set_Pa']:g}_N{r['N']}"
        for k, v in r.pop("series", {}).items():
            series[f"{tag}/{k}"] = v
    np.savez_compressed(os.path.join(a.dir, "tp_sweep_series.npz"), **series)
    json.dump(dict(dir=a.dir, window_s=a.window, T_c_nominal_K=T_C_NOM, runs=rows),
              open(os.path.join(a.dir, "tp_sweep_compare.json"), "w"), indent=1)

    print(f"{'T':>4} {'p':>7} {'N':>5} {'t_ready':>7} {'F_sim':>9} {'E_meas':>6} {'E_nom':>6}"
          f" {'J_meas':>6} {'J_nom':>6} {'L_km':>6} {'drift':>6} {'E_mEnd':>6} {'J_nEnd':>6}")
    for r in rows:
        if r.get("excluded"):
            print(f"{r['T_wall_K']/1e3:4.0f} {r['p_set_Pa']:7.4g} {r['N']:5d}  EXCLUDED: {r['excluded']}")
            continue
        tr = f"{r['t_ready_s']:.0f}" if r["t_ready_s"] else "?"
        dr = f"{r['F_drift_per_500s']:.3f}" if r["F_drift_per_500s"] is not None else "?"
        print(f"{r['T_wall_K']/1e3:4.0f} {r['p_set_Pa']:7.4g} {r['N']:5d} {tr:>7}{'' if r['ready'] else '*'}"
              f" {r['F_sim']:9.3e} {r['ratio_E_meas']:6.3f} {r['ratio_E_nom']:6.3f}"
              f" {r['ratio_J_meas']:6.3f} {r['ratio_J_nom']:6.3f} {r['L_TR_km']:6.1f} {dr:>6}"
              f" {r.get('ratio_E_meas_end', np.nan):6.3f} {r.get('ratio_J_nom_end', np.nan):6.3f}")

    good = [r for r in rows if not r.get("excluded") and r["N"] == 500]
    if not good:
        return
    Ts = sorted({r["T_wall_K"] for r in good})
    Ps = sorted({r["p_set_Pa"] for r in good})
    rdy = grid_of(good, "ready", Ts, Ps).astype(float) > 0.5
    Fs = grid_of(good, "F_sim", Ts, Ps)
    FE = grid_of(good, "F_E_nom", Ts, Ps)
    fl = LogNorm(np.nanmin([Fs, FE]), np.nanmax([Fs, FE]))
    rnorm = TwoSlopeNorm(1.0, 0.6, 1.4)
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    panel(ax[0, 0], Fs, Ts, Ps, r"$F_{sim}$ [kg m$^{-2}$ s$^{-1}$] (top cell)", fl, "viridis", ".2e", rdy)
    panel(ax[0, 1], FE, Ts, Ps, r"$F_E = q_{wall}/(h_{top}-h_c)$, $T_c$ = 6.7 kK", fl, "viridis", ".2e", rdy)
    panel(ax[0, 2], grid_of(good, "ratio_E_nom", Ts, Ps), Ts, Ps,
          r"$F_{sim}/F_E$ (nominal $T_c$)", rnorm, "RdBu_r", ".3f", rdy)
    panel(ax[1, 0], grid_of(good, "ratio_E_meas", Ts, Ps), Ts, Ps,
          r"$F_{sim}/F_E$ (measured $h_c$)", rnorm, "RdBu_r", ".3f", rdy)
    panel(ax[1, 1], grid_of(good, "ratio_J_nom", Ts, Ps), Ts, Ps,
          r"$F_{sim}/F_J = F L_{TR}/J$ (nominal $T_c$)", rnorm, "RdBu_r", ".3f", rdy)
    panel(ax[1, 2], grid_of(good, "t_ready_s", Ts, Ps), Ts, Ps,
          r"$t_{ready}$ [s]  (* = not ready before the end)", None, "magma", ".0f", rdy)
    fig.suptitle("model_column (T_wall, p_top) sweep: simulated vs analytic evaporation mass flux, "
                 f"read over {a.window:.0f} s from t_ready", fontsize=11)
    fig.tight_layout()
    out = os.path.join(a.out_dir, f"{a.prefix}_flux_compare.png")
    os.makedirs(a.out_dir, exist_ok=True)
    fig.savefig(out, dpi=140)
    print("wrote", out)


if __name__ == "__main__":
    main()
