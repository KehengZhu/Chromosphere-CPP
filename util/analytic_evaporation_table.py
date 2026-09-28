#!/usr/bin/env python3
"""Analytic steady-evaporation lookup table J(T_top, p_top, T_c) for the release column.

Post-processing only (no new run). The steady, radiation-free, isobaric,
thin transition region of the release equations obeys kappa dT/ds = F (h - h_c),
so the product of the evaporation mass flux F and the TR length L is fixed by
the boundary state alone (twall_mass_flux_relaxation_timescales.md, eq. 1):

    J = F * L_{T1} = int_{T1}^{T_top} kappa(T, p) / (h(T, p) - h_c) dT,

with the chromospheric enthalpy closed here as h_c = h_Saha(T_c, p_top). J is
what the table stores; F needs a length L that TR physics does not set in the
release. The script

  1. tabulates J, h_c, x_c and the local log-derivatives on a (T_top, p_top, T_c)
     grid;
  2. checks the closure end to end against the accepted T_wall sweep runs
     at 4000 s: the table at the run's (T_top, p_top) against the measured F*L
     and against the J the diagnostic formed from the measured h_c;
  3. measures the sensitivity to T_c and to the TR-base marker T1.

Study: docs/studies/boundaries/analytic_evaporation_J_table.md

Usage:
    python util/analytic_evaporation_table.py [--dir outputs/model_column/twall_sweep]
        [--out-dir outputs/model_column/analytic_J_table]
        [--fig visualization/model_column/analytic_J_table.png]
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from conductive_flux import physical_heat_flux
from twall_relaxation_diag import (CHI, MH, T_MARK, geometry, isotherm_height,
                                   j_theory, saha_rho, saha_x)
from twall_sweep_diag import COLUMNS, load_gamma_diag, load_outercond

ROOT = os.path.dirname(HERE)
I_RHO, I_T, I_X, I_NE, I_NHI, I_P = (
    COLUMNS.index(k) for k in ("rho", "T", "x_eq", "n_e", "n_HI", "p"))

# Grid. T_top 15-100 kK is the range the sweep validated (F L/J = 0.89-1.01);
# 100-200 kK and every p_top other than the release 0.0103 Pa are unvalidated.
T_TOP = np.geomspace(15.0e3, 200.0e3, 29)          # [K]
P_TOP = np.geomspace(1.0e-3, 10.0, 33)             # [Pa]
T_C = np.array([6.0e3, 6.35e3, 6.7e3, 7.0e3])      # [K]; 6.7 kK = C7 plateau
T_C_NOMINAL = 6.7e3
P_RELEASE = 0.0102856                              # kOuterPRef of the release IC
T_VALID_MAX = 100.0e3


def enthalpy(T, p):
    """Saha specific enthalpy h = 5/2 p/rho + x chi/m_H [J/kg] at (T, p)."""
    rho = saha_rho(np.asarray(p, float), np.asarray(T, float))
    return 2.5 * p / rho + saha_x(rho, T) * CHI / MH, saha_x(rho, T)


def build_table(t_mark=T_MARK):
    nT, nP, nC = len(T_TOP), len(P_TOP), len(T_C)
    J = np.empty((nT, nP, nC))
    hc = np.empty((nP, nC))
    xc = np.empty((nP, nC))
    for ip, p in enumerate(P_TOP):
        for ic, tc in enumerate(T_C):
            hc[ip, ic], xc[ip, ic] = enthalpy(tc, p)
            for it, tt in enumerate(T_TOP):
                J[it, ip, ic] = j_theory(tt, p, hc[ip, ic], T_lo=t_mark)
    return J, hc, xc


def run_state_4000(path):
    """-> (T_top, p_top, h_c, T at the h_c cell) at the last frame, as the diag forms h_c."""
    h_face, t, d = load_gamma_diag(path)
    oc = load_outercond(path)
    _, h = geometry(h_face)
    T, p, rho, x = d[-1, :, I_T], d[-1, :, I_P], d[-1, :, I_RHO], d[-1, :, I_X]
    qw = np.interp(t[-1], oc["t"], oc["q_face"])
    q = physical_heat_flux(h, T, d[-1, :, I_NE], d[-1, :, I_NHI])
    h75 = isotherm_height(h, T, T_MARK)
    j = j0 = int(np.searchsorted(h, h75))
    while j > 0 and abs(q[j]) > 0.01 * abs(qw) and h[j0] - h[j] < 30.0:
        j -= 1
    hen = 2.5 * p[j] / rho[j] + x[j] * CHI / MH
    return float(T[-1]), float(p[-1]), float(hen), float(T[j]), float(p[j])


def validate(sweep_dir):
    metrics = json.load(open(os.path.join(sweep_dir, "twall_relaxation_metrics.json")))
    rows = []
    for r in metrics["runs"]:
        if not r.get("J_theory_4000"):
            continue
        tw = r["T_wall_K"]
        cand = sorted(glob.glob(os.path.join(sweep_dir, f"twall_{int(tw/1e3)}kK_cfl*.txt")))
        if not cand:
            continue
        T_top, p_top, hc_run, Tc_run, pc_run = run_state_4000(cand[0])
        FL = r["FL_over_J_4000"] * r["J_theory_4000"]
        hc_nom, _ = enthalpy(T_C_NOMINAL, p_top)
        rows.append(dict(
            T_wall_K=tw, T_top_K=T_top, p_top_Pa=p_top,
            T_c_run_K=Tc_run, p_c_run_Pa=pc_run, h_c_run=hc_run, h_c_nominal=float(hc_nom),
            FL_measured=FL, J_measured_hc=r["J_theory_4000"],
            J_table_nominal=j_theory(T_top, p_top, hc_nom),
            J_table_run_Tc=j_theory(T_top, p_top, float(enthalpy(Tc_run, p_top)[0]))))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", default=os.path.join(ROOT, "outputs/model_column/twall_sweep"))
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "outputs/model_column/analytic_J_table"))
    ap.add_argument("--fig", default=os.path.join(ROOT, "visualization/model_column/analytic_J_table.png"))
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    os.makedirs(os.path.dirname(a.fig), exist_ok=True)

    J, hc, xc = build_table()
    J10, _, _ = build_table(t_mark=10.0e3)          # marker sensitivity
    lT, lP = np.log(T_TOP), np.log(P_TOP)
    nT_exp = np.gradient(np.log(J), lT, axis=0)      # d ln J / d ln T_top
    nP_exp = np.gradient(np.log(J), lP, axis=1)      # d ln J / d ln p_top
    ic0 = int(np.argmin(abs(T_C - T_C_NOMINAL)))

    np.savez(os.path.join(a.out_dir, "J_table.npz"), T_top_K=T_TOP, p_top_Pa=P_TOP,
             T_c_K=T_C, J=J, J_marker10kK=J10, h_c=hc, x_c=xc,
             dlnJ_dlnT=nT_exp, dlnJ_dlnp=nP_exp, T_mark_K=T_MARK)
    with open(os.path.join(a.out_dir, "J_table.csv"), "w") as f:
        f.write("# J = F*L_TR [kg m^-1 s^-1], steady radiation-free isobaric TR, "
                f"marker T1={T_MARK:.0f} K, h_c = h_Saha(T_c, p_top). "
                f"Validated only for 15<=T_top<=100 kK at p_top={P_RELEASE} Pa, T_c=6.7 kK.\n")
        f.write("T_top_K,p_top_Pa,T_c_K,J,h_c,x_c,dlnJ_dlnT,dlnJ_dlnp\n")
        for it, tt in enumerate(T_TOP):
            for ip, p in enumerate(P_TOP):
                for ic, tc in enumerate(T_C):
                    f.write(f"{tt:.6g},{p:.6g},{tc:.6g},{J[it,ip,ic]:.6e},{hc[ip,ic]:.6e},"
                            f"{xc[ip,ic]:.6f},{nT_exp[it,ip,ic]:.4f},{nP_exp[it,ip,ic]:.4f}\n")

    rows = validate(a.dir)

    # Summary numbers quoted in the study note.
    summary = dict(validation=rows, T_c_sensitivity={}, marker_sensitivity={},
                   p_exponent={}, T_exponent={})
    for tt in (15e3, 22e3, 50e3, 100e3):
        it = int(np.argmin(abs(T_TOP - tt)))
        for p in (1e-3, P_RELEASE, 0.1, 1.0, 10.0):
            ip = int(np.argmin(abs(P_TOP - p)))
            key = f"T{T_TOP[it]/1e3:.0f}kK_p{P_TOP[ip]:.3g}"
            summary["T_c_sensitivity"][key] = {f"{tc:.0f}": float(J[it, ip, ic] / J[it, ip, ic0])
                                               for ic, tc in enumerate(T_C)}
            summary["marker_sensitivity"][key] = float(J10[it, ip, ic0] / J[it, ip, ic0])
            summary["p_exponent"][key] = float(nP_exp[it, ip, ic0])
            summary["T_exponent"][key] = float(nT_exp[it, ip, ic0])
    json.dump(summary, open(os.path.join(a.out_dir, "J_table_summary.json"), "w"), indent=1)

    print("validation at 4000 s (T_c nominal 6.7 kK, closure h_c = h(T_c, p_top)):")
    print(" T_wall  T_c,run  hc_run/hc_nom  Jtab/FLmeas  Jtab/Jmeas(hc)  Jtab(Tc_run)/Jmeas(hc)")
    for r in rows:
        print(f" {r['T_wall_K']/1e3:5.0f}  {r['T_c_run_K']:7.0f}  {r['h_c_run']/r['h_c_nominal']:12.3f}"
              f"  {r['J_table_nominal']/r['FL_measured']:11.3f}  {r['J_table_nominal']/r['J_measured_hc']:14.3f}"
              f"  {r['J_table_run_Tc']/r['J_measured_hc']:22.3f}")
    for k in summary["T_c_sensitivity"]:
        print(k, "T_c ratios", {c: round(v, 3) for c, v in summary["T_c_sensitivity"][k].items()},
              "marker10/7.5", round(summary["marker_sensitivity"][k], 3),
              "dlnJ/dlnp", round(summary["p_exponent"][k], 3),
              "dlnJ/dlnT", round(summary["T_exponent"][k], 3))

    # ---- figure ---------------------------------------------------------------
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
    cmap = plt.get_cmap("viridis")
    p_show = [1e-3, 1e-2, 0.1, 1.0, 10.0]
    a0 = ax[0, 0]
    for k, p in enumerate(p_show):
        ip = int(np.argmin(abs(P_TOP - p)))
        a0.loglog(T_TOP / 1e3, J[:, ip, ic0], color=cmap(k / (len(p_show) - 1)),
                  label=f"$p_{{top}}$ = {P_TOP[ip]:.3g} Pa")
    if rows:
        a0.loglog([r["T_top_K"] / 1e3 for r in rows], [r["FL_measured"] for r in rows], "k^",
                  ms=6, label=r"sweep $F L_{TR}$ at 4000 s ($p_{top}$=0.0103 Pa)")
    a0.axvspan(T_VALID_MAX / 1e3, T_TOP[-1] / 1e3, color="0.9", zorder=0)
    a0.set(xlabel=r"$T_{top}$ [kK]", ylabel=r"$J = F\,L_{TR}$ [kg m$^{-1}$ s$^{-1}$]",
           title=r"(a) $J(T_{top})$ at $T_c$ = 6.7 kK")
    a0.legend(fontsize=7)

    a1 = ax[0, 1]
    T_show = [15e3, 22e3, 50e3, 100e3, 200e3]
    for k, tt in enumerate(T_show):
        it = int(np.argmin(abs(T_TOP - tt)))
        a1.semilogx(P_TOP, J[it, :, ic0] / J[it, np.argmin(abs(P_TOP - P_RELEASE)), ic0],
                    color=cmap(k / (len(T_show) - 1)), label=f"$T_{{top}}$ = {T_TOP[it]/1e3:.0f} kK")
    a1.axvline(P_RELEASE, color="k", ls=":", lw=1)
    a1.set(xlabel=r"$p_{top}$ [Pa]", ylabel=r"$J / J(p_{top}=0.0103\,$Pa$)$",
           title=r"(b) pressure dependence at $T_c$ = 6.7 kK")
    a1.legend(fontsize=7)

    a2 = ax[1, 0]
    ip0 = int(np.argmin(abs(P_TOP - P_RELEASE)))
    ls = ["-", "--"]
    for jp, ip in enumerate([ip0, int(np.argmin(abs(P_TOP - 1.0)))]):
        for ic, tc in enumerate(T_C):
            if ic == ic0:
                continue
            a2.semilogx(T_TOP / 1e3, J[:, ip, ic] / J[:, ip, ic0], ls[jp],
                        color=plt.get_cmap("coolwarm")(ic / (len(T_C) - 1)),
                        label=f"$T_c$ = {tc/1e3:.2f} kK, $p$ = {P_TOP[ip]:.3g} Pa")
    a2.axhline(1, color="k", lw=0.8)
    a2.set(xlabel=r"$T_{top}$ [kK]", ylabel=r"$J(T_c) / J(6.7\,$kK$)$",
           title=r"(c) sensitivity to the chromospheric temperature $T_c$")
    a2.legend(fontsize=7)

    a3 = ax[1, 1]
    cs = a3.contourf(P_TOP, T_TOP / 1e3, nP_exp[:, :, ic0], levels=21, cmap="RdBu_r",
                     vmin=-max(abs(nP_exp[:, :, ic0].min()), abs(nP_exp[:, :, ic0].max())),
                     vmax=max(abs(nP_exp[:, :, ic0].min()), abs(nP_exp[:, :, ic0].max())))
    fig.colorbar(cs, ax=a3, label=r"$\partial\ln J/\partial\ln p_{top}$")
    a3.contour(P_TOP, T_TOP / 1e3, nP_exp[:, :, ic0], levels=[0], colors="k", linewidths=0.8)
    # Where J depends on the TR-base marker by more than 20 % (T1 = 10 vs 7.5 kK):
    # there the ionization zone sits inside the TR and "L_TR" must name its isotherm.
    a3.contour(P_TOP, T_TOP / 1e3, J10[:, :, ic0] / J[:, :, ic0], levels=[0.8],
               colors="m", linestyles="--", linewidths=1.2)
    a3.plot([], [], "m--", label=r"$J(T_1{=}10\,$kK$)/J(7.5\,$kK$)$ = 0.8")
    a3.plot([P_RELEASE] * 2, [15, 100], "k-", lw=2.5, label="validated by the sweep")
    a3.set(xscale="log", yscale="log", xlabel=r"$p_{top}$ [Pa]", ylabel=r"$T_{top}$ [kK]",
           title=r"(d) local pressure exponent, $T_c$ = 6.7 kK")
    a3.legend(fontsize=7, loc="lower left")
    fig.suptitle(r"Analytic steady evaporation invariant $J=\int\kappa/(h-h_c)\,dT$ "
                 r"(release: no radiation, Saha, $T_1$ = 7.5 kK)", fontsize=11)
    fig.tight_layout()
    fig.savefig(a.fig, dpi=150)
    print("wrote", a.fig, "and", a.out_dir)


if __name__ == "__main__":
    main()
