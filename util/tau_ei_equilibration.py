#!/usr/bin/env python3
"""Measured electron-ion temperature-equilibration time for the model_column release.

Replaces a hand order-of-magnitude estimate with a per-cell measurement taken from
the release run's own EOS diagnostics. The `.gamma_diag` sidecar already carries the
decoded `T`, `x_eq`, `n_e` and `n_HI` (they are derived EOS diagnostics, never stored
state), so nothing here re-derives Saha: the ionization state is exactly the one the
run used.

The electron-ion energy-equilibration rate is the NRL Plasma Formulary expression

    nu_eps^{e\\i} = 3.2e-9 * Z^2 * n_i * lnLambda_ei / (mu * T_e^{3/2})   [s^-1]

with n_i in cm^-3, T_e in eV, mu = m_i/m_p = 1 and Z = 1 for hydrogen, and
tau_eq = 1 / nu_eps. The Coulomb logarithm is evaluated per cell on the NRL branch
that actually applies rather than pinned to a fixed value; the branch index is
reported, and a small or non-positive lnLambda is flagged because it means the
formula is being pushed outside its validity.

Two other relaxation channels are estimated for context, because the lower domain is
only partially ionized and electron-ion collisions are not automatically the fastest
route to a common temperature:

  * electron-neutral energy exchange, rate  (2 m_e/m_H) * n_HI * sigma_en * v_th,e
  * ion-neutral energy exchange, rate       n_HI * sigma_in * v_th,in

Both use a constant momentum-transfer cross-section (`--sigma-en`, `--sigma-in`);
they are deliberately order-of-magnitude, and only their ordering against tau_ei
is used.

Competing timescales are taken from the same run: the solver timestep `dt` (parsed
from the console log), the single-cell sound-crossing time `ds / c_s` with the
release frozen sound speed `sqrt(5/3 p/rho)`, the local conduction diffusion time
`ds^2 / chi` built from the run's own `kappa_physical` (thermal conduction is the
only electron-specific energy channel in the release, so this is the physically
relevant thing that could pull T_e away from T_i), and the evaporation /
mass-loading timescale, which is supplied as a given band (`--evap-band`).

The 1 MK row of the table is an EXTRAPOLATION, not a measurement: this domain tops
out near 22 kK and never reaches coronal temperatures. It is computed from
`--tr-ne` / `--tr-T` and labelled as such everywhere it appears.

Usage:
  .venv/bin/python util/tau_ei_equilibration.py
      [--run outputs/model_column/twall_sweep/twall_22kK_cfl0.50.txt]
      [--times 0 4000] [--heights 1600 1850 2153]
      [--sigma-en 1e-19] [--sigma-in 1e-18]
      [--out visualization/model_column/tau_ei_equilibration.png]
      [--json outputs/model_column/tau_ei_equilibration.json]
"""
import argparse
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
from twall_sweep_diag import (IP, IRHO, IT, IV, IX, console_ok,  # noqa: E402
                              load_gamma_diag)

INE = 4          # .gamma_diag column indices not exported by twall_sweep_diag
INHI = 5
IKAPPA = 8

K_B_EV = 8.617333262e-5      # Boltzmann constant [eV K^-1]
K_B_J = 1.380649e-23         # Boltzmann constant [J K^-1]
M_E = 9.1093837015e-31       # electron mass [kg]
M_H = 1.67262192369e-27      # hydrogen/proton mass [kg]
GAMMA_FROZEN = 5.0 / 3.0     # the release Riemann solver's frozen index


# --------------------------------------------------------------------------
# NRL Plasma Formulary, "Collisions and Transport"
# --------------------------------------------------------------------------
def coulomb_log_ei(n_e_cgs, T_e_eV, T_i_eV, Z=1.0, mu=1.0):
    """NRL electron-ion Coulomb logarithm -> (lnLambda, branch index 1|2|3).

    Branches, in the formulary's own order (n in cm^-3, T in eV):
      1  T_i m_e/m_i < T_e < 10 Z^2      : 23 - ln(n_e^{1/2} Z T_e^{-3/2})
      2  T_i m_e/m_i < 10 Z^2 < T_e      : 24 - ln(n_e^{1/2} / T_e)
      3  T_e < T_i Z m_e/m_i             : 30 - ln(n_i^{1/2} T_i^{-3/2} Z^2 / mu)
    """
    n_e_cgs = np.asarray(n_e_cgs, dtype=float)
    T_e_eV = np.asarray(T_e_eV, dtype=float)
    T_i_eV = np.asarray(T_i_eV, dtype=float)
    ratio = M_E / (mu * M_H)

    b1 = 23.0 - np.log(np.sqrt(n_e_cgs) * Z * T_e_eV ** -1.5)
    b2 = 24.0 - np.log(np.sqrt(n_e_cgs) / T_e_eV)
    b3 = 30.0 - np.log(np.sqrt(n_e_cgs / Z) * T_i_eV ** -1.5 * Z ** 2 / mu)

    branch = np.where(T_e_eV < T_i_eV * Z * ratio, 3,
                      np.where(T_e_eV < 10.0 * Z ** 2, 1, 2))
    lam = np.where(branch == 3, b3, np.where(branch == 1, b1, b2))
    return lam, branch


def tau_ei(n_i_cgs, T_e_eV, lam, Z=1.0, mu=1.0):
    """NRL electron-ion energy-equilibration time 1/nu_eps^{e\\i} [s]."""
    nu = 3.2e-9 * Z ** 2 * np.asarray(n_i_cgs) * lam / (mu * np.asarray(T_e_eV) ** 1.5)
    return 1.0 / nu


def tau_en_energy(n_HI_si, T_K, sigma_en):
    """Electron-neutral ENERGY relaxation time [s].

    Momentum-transfer frequency nu_en = n_HI sigma_en <v_e>, with the mean thermal
    speed <v_e> = sqrt(8 k T / pi m_e); an electron loses a fraction ~2 m_e/m_H of
    its energy per elastic collision, so the energy rate is (2 m_e/m_H) nu_en.
    """
    v_e = np.sqrt(8.0 * K_B_J * np.asarray(T_K) / (np.pi * M_E))
    nu = np.asarray(n_HI_si) * sigma_en * v_e
    return 1.0 / (2.0 * (M_E / M_H) * nu)


def tau_in_energy(n_HI_si, T_K, sigma_in):
    """Ion-neutral ENERGY relaxation time [s], reduced mass m_H/2, O(1) per collision."""
    v = np.sqrt(8.0 * K_B_J * np.asarray(T_K) / (np.pi * (M_H / 2.0)))
    return 1.0 / (np.asarray(n_HI_si) * sigma_in * v)


# --------------------------------------------------------------------------
def cell_geometry(h_km):
    """Face heights (.gamma_diag line 2) -> (centre_km, ds_m).

    Line 2 of the sidecar holds the UPPER-FACE height of each cell (io_formats.md),
    so widths are the forward difference and cell 0 inherits the width of cell 1
    (the mesh is uniform at the base; only the top is refined).
    """
    h = np.asarray(h_km, dtype=float)
    ds_km = np.empty_like(h)
    ds_km[1:] = np.diff(h)
    ds_km[0] = ds_km[1]
    return h - 0.5 * ds_km, ds_km * 1.0e3


def dt_from_log(path):
    """-> (dt_min, dt_max, dt_last) [s] from the run console log, or (None,)*3."""
    log = path.replace(".txt", ".console.log")
    if not os.path.exists(log):
        return None, None, None
    dts = [float(m) for m in re.findall(r"dt = (\S+)", open(log, errors="replace").read())]
    if not dts:
        return None, None, None
    return min(dts), max(dts), dts[-1]


def pick(times, t):
    return int(np.argmin(np.abs(times - t)))


def analyse(frame, h_km, sigma_en, sigma_in):
    """One .gamma_diag frame -> dict of per-cell derived quantities."""
    centre_km, ds_m = cell_geometry(h_km)
    T = frame[:, IT]
    x = frame[:, IX]
    n_e_si = frame[:, INE]
    n_HI_si = frame[:, INHI]
    rho = frame[:, IRHO]
    p = frame[:, IP]

    n_e_cgs = n_e_si * 1.0e-6
    T_eV = T * K_B_EV
    lam, branch = coulomb_log_ei(n_e_cgs, T_eV, T_eV)
    tei = tau_ei(n_e_cgs, T_eV, lam)          # n_i = n_e for Z = 1 hydrogen

    cs = np.sqrt(GAMMA_FROZEN * p / rho)

    # Local conduction diffusion time across one cell, ds^2 / chi, with the thermal
    # diffusivity chi = kappa / ((3/2) n_tot k_B) and n_tot = (1+x) n_H the total
    # particle density the release EOS uses for p_total. Conduction is the only
    # electron-specific energy channel in this model, so this is the fastest
    # PHYSICAL process that could separate T_e from T_i.
    n_tot = n_e_si + n_HI_si + n_e_si          # electrons + neutrals + protons
    chi = frame[:, IKAPPA] / (1.5 * n_tot * K_B_J)

    return dict(centre_km=centre_km, ds_m=ds_m, T=T, x=x, v=frame[:, IV],
                n_e_si=n_e_si, n_e_cgs=n_e_cgs, n_HI_si=n_HI_si, T_eV=T_eV,
                lnL=lam, branch=branch, tau_ei=tei,
                tau_en=tau_en_energy(n_HI_si, T, sigma_en),
                tau_in=tau_in_energy(n_HI_si, T, sigma_in),
                cs=cs, tau_cell=ds_m / cs, tau_cond=ds_m ** 2 / chi)


def row(a, i):
    return dict(h_km=float(a["centre_km"][i]), T=float(a["T"][i]), x=float(a["x"][i]),
                n_e_cgs=float(a["n_e_cgs"][i]), n_e_si=float(a["n_e_si"][i]),
                lnL=float(a["lnL"][i]), branch=int(a["branch"][i]),
                tau_ei=float(a["tau_ei"][i]), tau_en=float(a["tau_en"][i]),
                tau_in=float(a["tau_in"][i]), tau_cell=float(a["tau_cell"][i]),
                tau_cond=float(a["tau_cond"][i]))


# --------------------------------------------------------------------------
def make_figure(a0, a1, t0, t1, dt, evap_band, out_path):
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(8.2, 8.4), sharex=True,
                                 gridspec_kw=dict(height_ratios=[2.1, 1.0]))

    lo, hi = evap_band
    ax.axhspan(lo, hi, color="#9b59b6", alpha=0.13, zorder=0)
    ax.text(a1["centre_km"][-1], np.sqrt(lo * hi), "evaporation / mass loading  ",
            color="#6c3483", fontsize=8, va="center", ha="right")
    if dt is not None:
        ax.axhline(dt, color="#7f8c8d", ls=":", lw=1.4)
        ax.text(a1["centre_km"][-1], dt * 1.18,
                f"solver $\\Delta t$ = {dt*1e3:.2f} ms  ",
                color="#566573", fontsize=8, ha="right")

    ax.plot(a1["centre_km"], a1["tau_cell"], color="#34495e", ls="-.", lw=1.3,
            label=r"cell sound crossing $\Delta s/c_s$")
    ax.plot(a1["centre_km"], a1["tau_cond"], color="#e67e22", ls="-.", lw=1.3,
            label=r"cell conduction $\Delta s^2/\chi$")
    for a, ls, tag in ((a0, "--", f"t = {t0:g} s"), (a1, "-", f"t = {t1:g} s")):
        ax.plot(a["centre_km"], a["tau_ei"], color="#c0392b", ls=ls, lw=1.8,
                label=rf"$\tau_{{eq}}^{{ei}}$ (NRL), {tag}")
        ax.plot(a["centre_km"], a["tau_en"], color="#1f6fb2", ls=ls, lw=1.3,
                label=rf"$\tau^{{en}}_\varepsilon$, {tag}")
        ax.plot(a["centre_km"], a["tau_in"], color="#27ae60", ls=ls, lw=1.3,
                label=rf"$\tau^{{in}}_\varepsilon$, {tag}")
    ax.set_yscale("log")
    # tau_en / tau_in diverge as n_HI -> 0 in the fully ionized top cells; clip so the
    # decade that actually decides the argument stays readable.
    ax.set_ylim(1.0e-5, 1.0e5)
    ax.set_ylabel("relaxation / competing timescale  [s]")
    ax.set_title("model_column release: measured temperature-equilibration times\n"
                 "vs the timescales the solver actually resolves", fontsize=11)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.5, ncol=3, loc="upper left", framealpha=0.92)
    j = int(np.argmin(a1["tau_cell"] / a1["tau_ei"]))
    ax.annotate(f"worst margin\n$\\Delta s/c_s$ / $\\tau_{{eq}}^{{ei}}$ = "
                f"{(a1['tau_cell']/a1['tau_ei'])[j]:.1f}x",
                xy=(a1["centre_km"][j], a1["tau_ei"][j]),
                xytext=(a1["centre_km"][j] - 260, 2.0e-5), fontsize=8, color="#7b241c",
                arrowprops=dict(arrowstyle="->", color="#7b241c", lw=1.0))

    bx.plot(a1["centre_km"], a1["n_e_cgs"], color="#c0392b", lw=1.7,
            label=rf"$n_e$, t = {t1:g} s")
    bx.plot(a0["centre_km"], a0["n_e_cgs"], color="#c0392b", lw=1.3, ls="--",
            label=rf"$n_e$, t = {t0:g} s")
    bx.axhspan(5e9, 1e11, color="#f39c12", alpha=0.16, zorder=0)
    bx.text(a1["centre_km"][5], 2.2e10, "hand-estimate range", color="#9a6209", fontsize=8)
    bx.set_yscale("log")
    bx.set_ylabel(r"$n_e$  [cm$^{-3}$]")
    bx.set_xlabel("height [km]")
    bx.grid(alpha=0.25, which="both")
    cx = bx.twinx()
    cx.plot(a1["centre_km"], a1["x"], color="#1f6fb2", lw=1.4)
    cx.set_ylabel("ionization fraction $x$", color="#1f6fb2")
    cx.tick_params(axis="y", colors="#1f6fb2")
    bx.legend(fontsize=8, loc="upper left")

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    return out_path


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default="outputs/model_column/twall_sweep/twall_22kK_cfl0.50.txt")
    ap.add_argument("--times", nargs=2, type=float, default=[0.0, 4000.0],
                    metavar=("EARLY", "LATE"))
    ap.add_argument("--heights", nargs="*", type=float, default=[1600.0, 1850.0, 2153.0])
    ap.add_argument("--sigma-en", type=float, default=1.0e-19,
                    help="e-H momentum-transfer cross-section [m^2] (default 1e-19)")
    ap.add_argument("--sigma-in", type=float, default=1.0e-18,
                    help="H+-H (charge-exchange dominated) cross-section [m^2]")
    ap.add_argument("--evap-band", nargs=2, type=float, default=[1.0e2, 1.0e3])
    ap.add_argument("--tr-ne", type=float, default=1.0e9,
                    help="EXTRAPOLATION only: upper-TR n_e [cm^-3]")
    ap.add_argument("--tr-T", type=float, default=1.0e6,
                    help="EXTRAPOLATION only: upper-TR T [K]")
    ap.add_argument("--out", default="visualization/model_column/tau_ei_equilibration.png")
    ap.add_argument("--json", default="outputs/model_column/tau_ei_equilibration.json")
    args = ap.parse_args()

    h_km, times, data = load_gamma_diag(args.run)
    term, fallbacks = console_ok(args.run)
    dt_min, dt_max, dt_last = dt_from_log(args.run)
    i0, i1 = pick(times, args.times[0]), pick(times, args.times[1])
    t0, t1 = float(times[i0]), float(times[i1])
    a0 = analyse(data[i0], h_km, args.sigma_en, args.sigma_in)
    a1 = analyse(data[i1], h_km, args.sigma_en, args.sigma_in)

    print(f"run          : {args.run}")
    print(f"cells        : {len(h_km)}   termination={term}  godunov_fallbacks={fallbacks}")
    print(f"frames       : t = {times[0]:g} .. {times[-1]:g} s, n = {len(times)}")
    if dt_last is not None:
        print(f"solver dt    : min {dt_min*1e3:.4f} ms  max {dt_max*1e3:.4f} ms  "
              f"last {dt_last*1e3:.4f} ms")
    print(f"sigma_en     : {args.sigma_en:.1e} m^2    sigma_in: {args.sigma_in:.1e} m^2")
    nb = np.unique(a1["branch"])
    print(f"NRL lnLambda branch(es) in force at t = {t1:g} s: {list(map(int, nb))} "
          "(1 = 23 - ln(n_e^1/2 Z T_e^-3/2), the low-T_e branch)")

    hdr = (f"{'h[km]':>8} {'T[K]':>9} {'x':>10} {'n_e[cm^-3]':>12} {'n_e[m^-3]':>11} "
           f"{'lnL':>6} {'br':>3} {'tau_ei[s]':>11} {'tau_en[s]':>10} {'tau_in[s]':>10} "
           f"{'ds/c_s[s]':>10} {'ds^2/chi[s]':>11}")

    def show(a, tag):
        print(f"\n--- t = {tag} ---")
        print(hdr)
        for hh in args.heights:
            i = int(np.argmin(np.abs(a["centre_km"] - hh)))
            r = row(a, i)
            print(f"{r['h_km']:8.1f} {r['T']:9.1f} {r['x']:10.4g} {r['n_e_cgs']:12.4g} "
                  f"{r['n_e_si']:11.4g} {r['lnL']:6.2f} {r['branch']:3d} "
                  f"{r['tau_ei']:11.4g} {r['tau_en']:10.4g} {r['tau_in']:10.4g} "
                  f"{r['tau_cell']:10.4g} {r['tau_cond']:11.4g}")
        imin, imax = int(np.argmin(a["tau_ei"])), int(np.argmax(a["tau_ei"]))
        print(f"  min tau_ei = {a['tau_ei'][imin]:.4g} s at {a['centre_km'][imin]:.1f} km")
        print(f"  max tau_ei = {a['tau_ei'][imax]:.4g} s at {a['centre_km'][imax]:.1f} km")
        print(f"  n_e range  = {a['n_e_cgs'].min():.4g} .. {a['n_e_cgs'].max():.4g} cm^-3"
              f"   ({a['n_e_si'].min():.4g} .. {a['n_e_si'].max():.4g} m^-3)")
        print(f"  lnLambda   = {a['lnL'].min():.3f} .. {a['lnL'].max():.3f}"
              + ("   *** WARNING: lnLambda <= 2 somewhere; NRL expression is being "
                 "pushed outside its validity ***" if a["lnL"].min() <= 2.0 else ""))
        print(f"  x range    = {a['x'].min():.4g} .. {a['x'].max():.4g}")
        print(f"  T range    = {a['T'].min():.1f} .. {a['T'].max():.1f} K")
        print(f"  ds/c_s     = {a['tau_cell'].min():.4g} .. {a['tau_cell'].max():.4g} s")
        print(f"  ds^2/chi   = {a['tau_cond'].min():.4g} .. {a['tau_cond'].max():.4g} s")
        # Per-cell margins: the worst case is what the T_e = T_i assumption stands on.
        for name, comp in (("ds/c_s", a["tau_cell"]), ("ds^2/chi", a["tau_cond"])):
            ratio = comp / a["tau_ei"]
            j = int(np.argmin(ratio))
            print(f"  min {name}/tau_ei = {ratio[j]:.4g}x at {a['centre_km'][j]:.1f} km")
        if dt_last is not None:
            ratio = dt_last / a["tau_ei"]
            j = int(np.argmin(ratio))
            print(f"  min dt/tau_ei     = {ratio[j]:.4g}x at {a['centre_km'][j]:.1f} km")
        print(f"  min tau_evap/tau_ei = {args.evap_band[0]/a['tau_ei'].max():.4g}x "
              f"(evaporation floor {args.evap_band[0]:g} s)")

    show(a0, f"{t0:g} s")
    show(a1, f"{t1:g} s")

    # EXTRAPOLATION, not measured: this domain never exceeds ~22.6 kK.
    tr_T_eV = args.tr_T * K_B_EV
    tr_lam, tr_br = coulomb_log_ei(args.tr_ne, tr_T_eV, tr_T_eV)
    tr_tau = float(tau_ei(args.tr_ne, tr_T_eV, tr_lam))
    print(f"\n[EXTRAPOLATION - NOT FROM THIS RUN] transition region at T = {args.tr_T:.3g} K, "
          f"n_e = {args.tr_ne:.3g} cm^-3: lnLambda = {float(tr_lam):.2f} "
          f"(NRL branch {int(tr_br)}), tau_ei = {tr_tau:.4g} s")

    print(f"\nWorst-case margins at t = {t1:g} s (competing timescale / tau_ei, "
          "evaluated cell by cell):")
    cands = [("cell sound crossing ds/c_s", a1["tau_cell"] / a1["tau_ei"]),
             ("cell conduction ds^2/chi", a1["tau_cond"] / a1["tau_ei"]),
             ("evaporation floor 1e2 s", args.evap_band[0] / a1["tau_ei"])]
    if dt_last is not None:
        cands.insert(0, ("solver timestep dt", dt_last / a1["tau_ei"]))
    for name, ratio in cands:
        j = int(np.argmin(ratio))
        print(f"  {name:<28s} min {ratio[j]:10.4g}x  at {a1['centre_km'][j]:7.1f} km")

    out = make_figure(a0, a1, t0, t1, dt_last, args.evap_band, args.out)
    print(f"\nwrote {out}")

    if args.json:
        payload = dict(
            run=args.run, termination=term, godunov_fallbacks=fallbacks,
            ncells=len(h_km), sigma_en=args.sigma_en, sigma_in=args.sigma_in,
            dt_min_s=dt_min, dt_max_s=dt_max, dt_last_s=dt_last,
            extrapolated_TR=dict(measured=False, T_K=args.tr_T, n_e_cgs=args.tr_ne,
                                 lnLambda=float(tr_lam), branch=int(tr_br), tau_ei_s=tr_tau),
            frames={})
        for a, tt in ((a0, t0), (a1, t1)):
            payload["frames"][f"{tt:g}"] = dict(
                t_s=tt,
                rows=[row(a, int(np.argmin(np.abs(a["centre_km"] - hh))))
                      for hh in args.heights],
                tau_ei_min_s=float(a["tau_ei"].min()),
                tau_ei_min_h_km=float(a["centre_km"][int(np.argmin(a["tau_ei"]))]),
                tau_ei_max_s=float(a["tau_ei"].max()),
                tau_ei_max_h_km=float(a["centre_km"][int(np.argmax(a["tau_ei"]))]),
                n_e_cgs_min=float(a["n_e_cgs"].min()), n_e_cgs_max=float(a["n_e_cgs"].max()),
                lnL_min=float(a["lnL"].min()), lnL_max=float(a["lnL"].max()),
                x_min=float(a["x"].min()), x_max=float(a["x"].max()),
                T_min=float(a["T"].min()), T_max=float(a["T"].max()),
                tau_cell_min_s=float(a["tau_cell"].min()),
                tau_cond_min_s=float(a["tau_cond"].min()),
                tau_cond_max_s=float(a["tau_cond"].max()),
                margin_min_dt=(None if dt_last is None
                               else float((dt_last / a["tau_ei"]).min())),
                margin_min_sound=float((a["tau_cell"] / a["tau_ei"]).min()),
                margin_min_cond=float((a["tau_cond"] / a["tau_ei"]).min()),
                margin_min_evap=float(args.evap_band[0] / a["tau_ei"].max()),
                tau_en_min_s=float(a["tau_en"].min()), tau_en_max_s=float(a["tau_en"].max()),
                tau_in_min_s=float(a["tau_in"].min()), tau_in_max_s=float(a["tau_in"].max()))
        os.makedirs(os.path.dirname(args.json) or ".", exist_ok=True)
        with open(args.json, "w") as fh:
            json.dump(payload, fh, indent=2)
        print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
