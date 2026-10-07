#!/usr/bin/env python3
"""Anatomy of the analytic steady-TR invariant J = int_{T1}^{T_top} kappa/(h - h_c) dT.

Post-processing only (no run). Uses exactly the Saha EOS, conductivity and
quadrature of ``twall_relaxation_diag.j_theory`` at the nominal chromospheric
temperature T_c = 6.7 kK, and draws, for three top pressures of the
(T_wall, p_top) sweep:

  (a) the Saha ionization fraction x(T, p) and the enthalpy jump
      (h - h_c)/(h_top - h_c) that the integrand divides by;
  (b) the integrand kappa/(h - h_c), with its log divergence at T_c and its
      Spitzer T^{3/2} tail;
  (c) the cumulative integral J(T)/J(T_top), which in a steady TR is exactly the
      normalised distance above the T1 marker, i.e. the steady TR profile.

It also prints J against the fully ionized Spitzer limit and the share of
h_top - h_c gained between T_c and 10 kK (the "unresolved foot" share quoted in
docs/studies/boundaries/tp_sweep_evaporation_closure_validation.md).

Usage:
    python util/plot_J_integrand.py [--T-top 70000]
        [--out visualization/model_column/analytic_J_integrand.png]
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analytic_evaporation_table import enthalpy
from conductive_flux import physical_conductivity
from twall_relaxation_diag import KB, MH, T_MARK, j_theory, saha_rho, saha_x

ROOT = os.path.dirname(HERE)
T_C = 6.7e3
PRESSURES = (0.003, 0.0103, 0.08)          # [Pa]; ends and C7 line of the sweep
KAPPA0 = 9.2048e-12                        # Spitzer coefficient [W m^-1 K^-7/2]


def integrand(T, p, h_c):
    rho = saha_rho(np.full_like(T, p), T)
    x = saha_x(rho, T)
    nh = rho / MH
    kap = physical_conductivity(x * nh, (1.0 - x) * nh, T)
    h = enthalpy(T, p)[0]
    return kap / (h - h_c), x, h


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--T-top", type=float, default=70.0e3)
    ap.add_argument("--out", default=os.path.join(
        ROOT, "visualization/model_column/analytic_J_integrand.png"))
    a = ap.parse_args()
    colors = dict(zip(PRESSURES, ("#1f77b4", "#2ca02c", "#d62728")))

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    Tw = np.geomspace(6.0e3, a.T_top, 600)
    for p in PRESSURES:
        c = colors[p]
        h_c, x_c = enthalpy(T_C, p)
        h_top = enthalpy(a.T_top, p)[0]
        f, x, h = integrand(Tw, p, h_c)
        above = Tw > T_C
        ax[0].plot(Tw / 1e3, x, color=c, label=f"{p:g} Pa")
        ax[0].plot(Tw[above] / 1e3, (h[above] - h_c) / (h_top - h_c), color=c, ls="--")
        ax[1].plot(Tw[above] / 1e3, f[above], color=c, label=f"{p:g} Pa")
        Tj = np.geomspace(T_MARK, a.T_top, 300)      # the j_theory quadrature grid
        fj = integrand(Tj, p, h_c)[0]
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (fj[1:] + fj[:-1]) * np.diff(Tj))])
        J = j_theory(a.T_top, p, h_c)
        ax[2].plot(cum / cum[-1], Tj / 1e3, color=c, label=f"{p:g} Pa")
        share = (enthalpy(10.0e3, p)[0] - h_c) / (h_top - h_c)
        J_sp = 2.0 * KAPPA0 * MH * a.T_top**2.5 / (25.0 * KB)
        print(f"p={p:g} Pa: x_c={x_c:.3f} h_c={h_c:.3e} J={J:.3e} J/J_Spitzer={J/J_sp:.3f}"
              f" share(6.7-10 kK)={share:.3f}  quad-vs-cum={cum[-1]/J:.6f}")
    T32 = np.geomspace(2.0e4, a.T_top, 50)
    f_ref = integrand(np.array([a.T_top]), PRESSURES[1], enthalpy(T_C, PRESSURES[1])[0])[0][0]
    ax[1].plot(T32 / 1e3, f_ref * (T32 / a.T_top)**1.5, "k:", lw=1,
               label=r"$\propto T^{3/2}$ (Spitzer)")
    for axx in ax[:2]:
        axx.axvline(T_C / 1e3, color="0.5", lw=0.8, ls=":")
        axx.axvline(T_MARK / 1e3, color="0.5", lw=0.8, ls="-.")
        axx.set(xscale="log", xlabel="$T$ [kK]")
    ax[0].set(ylabel="fraction", ylim=(-0.02, 1.05),
              title=r"(a) solid $x$; dashed $(h-h_c)/(h_{top}-h_c)$ ($T_c$ = 6.7 kK)")
    ax[0].legend(fontsize=8, loc="center right")
    ax[1].set(yscale="log", ylabel=r"$\kappa/(h-h_c)$ [kg m$^{-1}$ s$^{-1}$ K$^{-1}$]",
              title=r"(b) integrand; dotted $T_c$, dash-dot $T_1$ = 7.5 kK")
    ax[1].legend(fontsize=8)
    ax[2].set(xlabel=r"$J(T)/J(T_{top})$ = $(s - s_1)/L_{TR}$",
              ylabel="$T$ [kK]",
              title=rf"(c) steady TR profile, $T_{{top}}$ = {a.T_top/1e3:.0f} kK")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=150)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
