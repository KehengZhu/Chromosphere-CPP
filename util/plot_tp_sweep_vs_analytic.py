#!/usr/bin/env python3
"""Simulated against analytic evaporation over the (T_wall, p_top) sweep, as curves.

Reads only the compact products of ``util/tp_sweep_compare.py``
(``tp_sweep_compare.json`` and ``tp_sweep_series.npz``), so it runs anywhere the
two files are, without the raw runs. Two figures:

``<prefix>_sim_vs_analytic.png``
  (a) J = F L_TR against p_top: simulation (markers) and the steady-TR integral
      J(T_wall, p_top, h_c(6.7 kK)) (lines); (b) parity of F_sim against the two
      table-usable closures, both with the nominal T_c; (c) measured over nominal
      chromospheric enthalpy h_c; (d)-(f) the three flux ratios against p_top, at
      t_ready (circles) and over the last window of the run (crosses).
``<prefix>_time_series.png``
  One column per T_wall: the energy ratio ER(t), the top-cell flux normalised by
  its median after 1000 s, and the relative standard deviation of the face flux
  F(z) after 1000 s, which locates where in height the flux is unsteady.

Open markers are runs that were not ready (t_ready + window > t_end); their
values are read over the last window and are not a steady-state measurement.

Usage:
    python util/plot_tp_sweep_vs_analytic.py [--dir outputs/model_column/tp_sweep]
        [--out-dir visualization/model_column] [--prefix tp_sweep]
"""
import argparse
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analytic_evaporation_table import enthalpy
from twall_relaxation_diag import j_theory
from twall_sweep_diag import wall_colors

ROOT = os.path.dirname(HERE)
T_C_NOM = 6.7e3
T_PROFILE0 = 1000.0            # must match PROFILE_T0 of tp_sweep_compare.py [s]


def tag(r):
    return f"T{r['T_wall_K']/1e3:.0f}_p{r['p_set_Pa']:g}_N{r['N']}"


def marker_kw(r, c):
    face = c if r["ready"] else "white"
    return dict(marker="o", ms=6, mfc=face, mec=c, mew=1.3, ls="none")


def ratio_panel(ax, runs, Ts, cols, key, title):
    for T in Ts:
        rr = sorted((r for r in runs if r["T_wall_K"] == T), key=lambda r: r["p_set_Pa"])
        p = [r["p_set_Pa"] for r in rr]
        ax.plot(p, [r[key] for r in rr], "-", color=cols[T], lw=1.0, alpha=0.6)
        for r in rr:
            ax.plot(r["p_set_Pa"], r[key], **marker_kw(r, cols[T]))
            ax.plot(r["p_set_Pa"], r[key + "_end"], "x", color=cols[T], ms=5, mew=1.2)
    ax.axhspan(0.95, 1.05, color="0.9", zorder=0)
    ax.axhline(1.0, color="k", lw=0.8)
    ax.set(xscale="log", xlabel=r"$p_{top}$ [Pa]", title=title, ylim=(0.80, 1.30))


def fig_compare(runs, ctrl, Ts, Ps, cols, out):
    fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))

    # (a) J against p: simulation and the steady-TR integral with nominal T_c
    a = ax[0, 0]
    pg = np.geomspace(0.8 * min(Ps), 1.25 * max(Ps), 40)
    hc_g = enthalpy(np.full_like(pg, T_C_NOM), pg)[0]
    for T in Ts:
        a.plot(pg, [j_theory(T, p, h) for p, h in zip(pg, hc_g)], "-", color=cols[T], lw=1.2,
               label=f"{T/1e3:.0f} kK")
        for r in runs:
            if r["T_wall_K"] == T:
                a.plot(r["p_set_Pa"], r["F_sim"] * r["L_TR_km"] * 1e3, **marker_kw(r, cols[T]))
    for r in ctrl:
        a.plot(r["p_set_Pa"], r["F_sim"] * r["L_TR_km"] * 1e3, "k*", ms=11, mfc="none",
               label="N=1000")
    a.set(xscale="log", yscale="log", xlabel=r"$p_{top}$ [Pa]",
          ylabel=r"$J$ [kg m$^{-1}$ s$^{-1}$]",
          title=r"(a) $J$: simulated $F L_{TR}$ (markers) vs integral, $T_c$ = 6.7 kK (lines)")
    wall_handles = a.get_legend_handles_labels()

    # (b) parity of the flux against the two table-usable closures
    b = ax[0, 1]
    for r in runs:
        c = cols[r["T_wall_K"]]
        for key, m in (("F_E_nom", "o"), ("F_J_nom", "s")):
            b.plot(r[key], r["F_sim"], marker=m, ms=5.5, mec=c, mew=1.2,
                   mfc=c if r["ready"] else "white", ls="none")
    lo, hi = 1e-10, 1.2e-8
    xx = np.geomspace(lo, hi, 2)
    b.fill_between(xx, 0.9 * xx, 1.1 * xx, color="0.9", zorder=0)
    b.plot(xx, xx, "k-", lw=0.8)
    b.plot([], [], "ko", mfc="none", label=r"$F_E = q_{wall}/(h_{top}-h_c)$")
    b.plot([], [], "ks", mfc="none", label=r"$F_J = J/L_{TR}$")
    b.set(xscale="log", yscale="log", xlim=(lo, hi), ylim=(lo, hi),
          xlabel=r"analytic $F$ [kg m$^{-2}$ s$^{-1}$]  ($T_c$ = 6.7 kK)",
          ylabel=r"$F_{sim}$ [kg m$^{-2}$ s$^{-1}$]",
          title=r"(b) parity, grey band $\pm$10 %")
    b.add_artist(b.legend(fontsize=8, loc="upper left"))
    b.legend(*wall_handles, fontsize=8, ncol=2, loc="lower right", title=r"$T_{wall}$",
             title_fontsize=8)

    # (c) chromospheric enthalpy the run actually has, against the nominal closure
    c_ = ax[0, 2]
    for T in Ts:
        rr = sorted((r for r in runs if r["T_wall_K"] == T), key=lambda r: r["p_set_Pa"])
        c_.plot([r["p_set_Pa"] for r in rr], [r["h_c_meas"] / r["h_c_nom"] for r in rr], "-",
                color=cols[T], lw=1.0, alpha=0.6)
        for r in rr:
            c_.plot(r["p_set_Pa"], r["h_c_meas"] / r["h_c_nom"], **marker_kw(r, cols[T]))
    c_.axhline(1.0, color="k", lw=0.8)
    c_.set(xscale="log", xlabel=r"$p_{top}$ [Pa]",
           title=r"(c) $h_{c,meas}/h_c(6.7\,\mathrm{kK},\,p_{top})$")

    ratio_panel(ax[1, 0], runs, Ts, cols, "ratio_E_meas",
                r"(d) $F_{sim}/F_E$, measured $h_c$ (= ER)")
    ratio_panel(ax[1, 1], runs, Ts, cols, "ratio_E_nom",
                r"(e) $F_{sim}/F_E$, nominal $T_c$")
    ratio_panel(ax[1, 2], runs, Ts, cols, "ratio_J_nom",
                r"(f) $F_{sim}/F_J = F L_{TR}/J$, nominal $T_c$")
    for r in ctrl:
        for axx, key in ((ax[1, 0], "ratio_E_meas"), (ax[1, 1], "ratio_E_nom"),
                         (ax[1, 2], "ratio_J_nom")):
            axx.plot(r["p_set_Pa"], r[key], "k*", ms=11, mfc="none")
    for axx in ax[1]:
        axx.set_xticks(Ps, [f"{p:g}" for p in Ps], fontsize=8)
        axx.minorticks_off()
    fig.suptitle("model_column (T_wall, p_top) sweep, N=500/R4, read over 300 s from t_ready. "
                 "Open markers: not ready before 5000 s (value from the last 300 s). "
                 "x: last 300 s (4700-5000 s). Star: N=1000 control.", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("wrote", out)


def fig_series(runs, z, Ts, Ps, out):
    pcols = dict(zip(Ps, plt.get_cmap("viridis")(np.linspace(0.0, 0.9, len(Ps)))))
    fig, ax = plt.subplots(3, len(Ts), figsize=(4.0 * len(Ts), 10.5), sharey="row")
    for j, T in enumerate(Ts):
        for r in sorted((r for r in runs if r["T_wall_K"] == T), key=lambda r: r["p_set_Pa"]):
            k, c = tag(r), pcols[r["p_set_Pa"]]
            t, F, ER = z[k + "/t"], z[k + "/F"], z[k + "/ER"]
            lab = f"{r['p_set_Pa']:g} Pa" + ("" if r["ready"] else " *")
            ax[0, j].plot(t, ER, color=c, lw=0.9, label=lab)
            if r["t_ready_s"] and r["ready"]:
                kk = np.argmin(np.abs(t - r["t_ready_s"]))
                ax[0, j].plot(t[kk], ER[kk], "o", color=c, ms=5, mec="k", mew=0.6)
            ax[1, j].plot(t, F / np.median(F[t >= T_PROFILE0]), color=c, lw=0.8)
            fz, fm, fs = z[k + "/face_km"], z[k + "/F_mean"], z[k + "/F_std"]
            ax[2, j].plot(fz, fs / np.abs(fm), color=c, lw=1.0)
        ax[0, j].set(title=f"$T_{{wall}}$ = {T/1e3:.0f} kK", ylim=(0.6, 2.2), xlabel="t [s]")
        ax[0, j].axhline(1.0, color="k", lw=0.6)
        ax[0, j].legend(fontsize=7, loc="upper right")
        ax[1, j].set(ylim=(-2.0, 4.0), xlabel="t [s]")
        ax[1, j].axhline(0.0, color="k", lw=0.6)
        ax[2, j].set(yscale="log", ylim=(3e-3, 30), xlabel="face height [km]")
    ax[0, 0].set_ylabel(r"ER = $F(h_{top}-h_c)/q_{wall}$" + "\n(dot: $t_{ready}$)")
    ax[1, 0].set_ylabel(r"$F_{top}(t)$ / median$_{t\geq1000\,s}$")
    ax[2, 0].set_ylabel(r"std$_t F(z)$ / |mean$_t F(z)$|, $t \geq$ 1000 s")
    fig.suptitle("model_column (T_wall, p_top) sweep, N=500/R4: TR energy balance, top-cell "
                 "flux and height-resolved flux unsteadiness (* = not ready before 5000 s)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", default=os.path.join(ROOT, "outputs/model_column/tp_sweep"))
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "visualization/model_column"))
    ap.add_argument("--prefix", default="tp_sweep")
    a = ap.parse_args()
    rows = json.load(open(os.path.join(a.dir, "tp_sweep_compare.json")))["runs"]
    rows = [r for r in rows if not r.get("excluded")]
    z = np.load(os.path.join(a.dir, "tp_sweep_series.npz"))
    runs = [r for r in rows if r["N"] == 500]
    ctrl = [r for r in rows if r["N"] != 500]
    Ts = sorted({r["T_wall_K"] for r in runs})
    Ps = sorted({r["p_set_Pa"] for r in runs})
    cols = wall_colors(Ts)[0]
    os.makedirs(a.out_dir, exist_ok=True)
    fig_compare(runs, ctrl, Ts, Ps, cols, os.path.join(a.out_dir, f"{a.prefix}_sim_vs_analytic.png"))
    fig_series(runs, z, Ts, Ps, os.path.join(a.out_dir, f"{a.prefix}_time_series.png"))


if __name__ == "__main__":
    main()
