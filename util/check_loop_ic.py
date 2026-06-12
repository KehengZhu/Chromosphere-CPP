#!/usr/bin/env python3
"""Sanity-check a full-loop .dat IC: plot T, densities, ionization fraction, and
the hydrostatic-equilibrium residual vs height. Confirms the chromosphere→TR→
corona profile is smooth and near-equilibrium before running the solver.

Usage: python util/check_loop_ic.py scenarios/data/loop_closed_test.dat
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M_I = 1.6726219e-27
K_B = 1.380649e-23
G_SI = 274.0


def parse_cells(path):
    rows = []
    in_cells = False
    for ln in open(path):
        s = ln.strip()
        if s.startswith("[CELLS]"):
            in_cells = True; continue
        if s.startswith("[GHOSTS]"):
            break
        if not in_cells or not s or s.startswith("#"):
            continue
        rows.append([float(x) for x in s.split()])
    a = np.array(rows)
    # i ds B_imh B_iph phi_imh phi_iph ne nn T
    ds = a[:, 1]; phi_imh = a[:, 4]; phi_iph = a[:, 5]
    ne = a[:, 6]; nn = a[:, 7]; T = a[:, 8]
    phi_c = 0.5 * (phi_imh + phi_iph)
    height_m = phi_c / G_SI                    # rise above footpoint
    return ds, height_m, ne, nn, T


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "scenarios/data/loop_closed_test.dat"
    out = sys.argv[2] if len(sys.argv) > 2 else "util/visualization/loop_ic_check.png"
    os.makedirs(os.path.dirname(out), exist_ok=True)

    ds, h, ne, nn, T = parse_cells(path)
    hMm = h / 1e6
    f_ion = ne / (ne + nn)
    # total pressure (ion 2x electron + neutral) and mass density
    p = (2.0 * ne + nn) * K_B * T
    rho = (ne + nn) * M_I
    # field-aligned gravity from d(phi)/ds is not stored; approximate g_par from
    # the local d(height)/ds (= cos angle for a loop) times g.
    s = np.concatenate(([0.0], np.cumsum(ds)))
    s_c = 0.5 * (s[:-1] + s[1:])
    dhds = np.gradient(h, s_c)
    g_par = G_SI * dhds
    # HSE residual: |dp/ds + rho g_par| / |rho g_par|
    dpds = np.gradient(p, s_c)
    resid = np.abs(dpds + rho * g_par) / np.maximum(np.abs(rho * g_par), 1e-30)

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(f"Loop IC check — {os.path.basename(path)}  ({len(ds)} cells)", fontsize=13)

    ax[0, 0].semilogy(hMm, T, "r"); ax[0, 0].set_ylabel("T [K]")
    ax[0, 0].set_xlabel("height above footpoint [Mm]"); ax[0, 0].grid(alpha=0.3)
    ax[0, 0].set_title("temperature (full loop)")
    # zoom inset: chromosphere + TR
    axins = ax[0, 0].inset_axes([0.45, 0.12, 0.5, 0.5])
    m = hMm < 1.5
    axins.semilogy(hMm[m], T[m], "r"); axins.grid(alpha=0.3)
    axins.set_title("chromosphere+TR (<1.5 Mm)", fontsize=8)

    ax[0, 1].semilogy(hMm, ne, "C3", label="$n_e$")
    ax[0, 1].semilogy(hMm, nn, "C0", label="$n_n$")
    ax[0, 1].set_ylabel("number density [m$^{-3}$]"); ax[0, 1].set_xlabel("height [Mm]")
    ax[0, 1].legend(); ax[0, 1].grid(alpha=0.3); ax[0, 1].set_title("densities")

    ax[1, 0].plot(hMm, f_ion, "purple"); ax[1, 0].set_ylabel("ionization fraction")
    ax[1, 0].set_xlabel("height [Mm]"); ax[1, 0].grid(alpha=0.3)
    ax[1, 0].axhline(0.99, color="k", ls="--", lw=0.8, label="F_SLAVE=0.99 (single-fluid)")
    ax[1, 0].legend(); ax[1, 0].set_title("ionization fraction")

    ax[1, 1].semilogy(hMm, np.clip(resid, 1e-4, None), "g")
    ax[1, 1].set_ylabel("|dp/ds + ρg∥| / |ρg∥|"); ax[1, 1].set_xlabel("height [Mm]")
    ax[1, 1].grid(alpha=0.3); ax[1, 1].set_title("HSE residual (corona should be ≪1)")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")
    cor = hMm > 2.0
    print(f"  corona HSE residual: median={np.median(resid[cor]):.2e}, max={resid[cor].max():.2e}")
    print(f"  TR location: T crosses 1e5 K at h≈{np.interp(1e5, T, hMm):.3f} Mm")


if __name__ == "__main__":
    main()
