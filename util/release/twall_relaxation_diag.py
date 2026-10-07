#!/usr/bin/env python3
"""Relaxation and decay timescales of the evaporation mass flux in the T_wall sweep.

Post-processing only (no new run). Reads the runs of ``util/run/run_twall_sweep.sh``
(``.gamma_diag``, ``.outercond`` and ``.faceflux`` sidecars) and measures, run by
run, the three stages the conservative mass flux F goes through, each against
its analytic estimate:

  1. wave phase     -- an acoustic pulse bouncing between the two reservoir faces:
                       one-way time t_ac = int ds/c_eq, ringing period, damping
                       time and the time after which F_base/F_top is within 5 %;
  2. TR relaxation  -- the C7 transition region relaxing to the conductive-
                       enthalpy balance, tested by the steady-TR invariant
                       F*L = J(T_top, p_top, h_c) (J computed from first
                       principles: Spitzer+neutral kappa and Saha enthalpy);
  3. secular drift  -- adiabatic contraction of the lifted chromosphere:
                       dL/dPhi = K = int delta(nabla_ad - nabla) g/p ds, so
                       F = J/L falls at the rate r = K F / L.

Two further tests of stage 3 are drawn: the chromospheric temperature change at
fixed column mass against the parcel-lift adiabat, and the Eulerian flux profile
F(z)/F_top = 1 - rho(z) int_base^z delta(nabla_ad - nabla) g/p ds, which is the
~12-13 % "monotone spread" of f_total over 1600-1850 km.

Derivations and results: docs/studies/boundaries/twall_mass_flux_relaxation_timescales.md

Usage:
    python util/release/twall_relaxation_diag.py [--dir outputs/model_column/twall_sweep]
        [--cfl 0.50] [--run T_WALL_K=path] [--release path]
        [--out-dir visualization/model_column] [--prefix twall_relaxation]
        [--json outputs/model_column/twall_sweep/twall_relaxation_metrics.json]
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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))  # util/lib
from animate_isentropic import load_face_flux
from conductive_flux import physical_conductivity, physical_heat_flux
from twall_sweep_diag import (COLUMNS, console_ok, load_gamma_diag,
                              load_outercond, wall_colors)

# Constants exactly as eos.hpp (pure-hydrogen Saha closure).
KB = 1.380649e-23
MH = 1.6726219e-27
ME = 9.1093837015e-31
HP = 6.62607015e-34
CHI = 2.179872361e-18
G_SUN = 273.95                 # model_column gravity [m s^-2]

I_RHO, I_V, I_T, I_X, I_NE, I_NHI, I_P, I_G1 = (
    COLUMNS.index(k) for k in ("rho", "v", "T", "x_eq", "n_e", "n_HI", "p", "Gamma1"))
T_MARK = 7500.0                # TR-base marker isotherm [K]
T_RELEASE = 22.0e3


# --------------------------------------------------------------------------
# Saha pure-hydrogen thermodynamics (reproduces the run's x and p exactly and
# its tabulated Gamma1 to 2e-4)
# --------------------------------------------------------------------------
def saha_x(rho, T):
    """Equilibrium ionization fraction x(rho, T)."""
    a = (2.0 * np.pi * ME * KB * T / HP**2)**1.5 * np.exp(-CHI / (KB * T)) / (rho / MH)
    return 0.5 * (-a + np.sqrt(a * a + 4.0 * a))


def saha_p(rho, T):
    return (1.0 + saha_x(rho, T)) * rho / MH * KB * T


def saha_e(rho, T):
    """Specific internal energy [J/kg], translational + ionization."""
    x = saha_x(rho, T)
    return 1.5 * (1.0 + x) * KB * T / MH + x * CHI / MH


def saha_rho(p, T, iters=60):
    """Density at given pressure and temperature (fixed point, vectorized)."""
    rho = p * MH / (KB * T)
    for _ in range(iters):
        rho = p * MH / ((1.0 + saha_x(rho, T)) * KB * T)
    return rho


def thermo(rho, T, eps=1.0e-6):
    """-> dict(G1, nab_ad, delta, cp) by central differences of the Saha EOS."""
    lr, lT = np.log(rho), np.log(T)
    lp = lambda a, b: np.log(saha_p(np.exp(a), np.exp(b)))
    ee = lambda a, b: saha_e(np.exp(a), np.exp(b))
    chi_rho = (lp(lr + eps, lT) - lp(lr - eps, lT)) / (2 * eps)
    chi_T = (lp(lr, lT + eps) - lp(lr, lT - eps)) / (2 * eps)
    cv = (ee(lr, lT + eps) - ee(lr, lT - eps)) / (2 * eps) / T
    p = saha_p(rho, T)
    g3m1 = p * chi_T / (rho * T * cv)
    g1 = chi_T * g3m1 + chi_rho
    return dict(G1=g1, nab_ad=g3m1 / g1, delta=chi_T / chi_rho, cp=cv * g1 / chi_rho)


def j_theory(T_top, p_top, h_c, T_lo=T_MARK, n=300):
    """Steady-TR invariant J = int_{T_lo}^{T_top} kappa/(h - h_c) dT [kg m^-1 s^-1].

    From the steady, radiation-free, isobaric TR energy balance
    kappa dT/ds = F (h - h_c), so that F * (s_top - s(T_lo)) = J.
    """
    T = np.geomspace(T_lo, T_top, n)
    rho = saha_rho(np.full_like(T, p_top), T)
    x = saha_x(rho, T)
    nh = rho / MH
    kap = physical_conductivity(x * nh, (1.0 - x) * nh, T)
    h = 2.5 * p_top / rho + x * CHI / MH
    return float(np.trapezoid(kap / (h - h_c), T))


# --------------------------------------------------------------------------
def geometry(face_km):
    """-> (cell widths [m], cell-centre heights [km]).

    The ``.gamma_diag`` height column is the cumulative upper-FACE height of each
    cell (it equals the ``.faceflux`` ``face_km`` column), not the cell centre, so
    the centres are rebuilt here from the face ladder.
    """
    f = np.asarray(face_km, dtype=float)
    e = np.concatenate([[f[0] - (f[1] - f[0])], f])
    return np.diff(e) * 1.0e3, 0.5 * (e[:-1] + e[1:])


def isotherm_height(h, T, thr):
    """Height [km] where T first exceeds thr, searching upward."""
    idx = np.where(T > thr)[0]
    if len(idx) == 0 or idx[0] == 0:
        return np.nan
    i = idx[0]
    return float(np.interp(thr, [T[i - 1], T[i]], [h[i - 1], h[i]]))


def smooth(t, y, width):
    """Centred running mean over a window of `width` seconds (irregular t)."""
    out = np.empty_like(y)
    for i, ti in enumerate(t):
        sel = np.abs(t - ti) <= 0.5 * width
        out[i] = y[sel].mean()
    return out


def contraction_weight(rho, T, p, w, mask):
    """int delta (nabla_ad - nabla) g/p ds over the masked cells, cumulative."""
    th = thermo(rho, T)
    nab = np.gradient(np.log(T), np.log(p))
    return np.cumsum(np.where(mask, th["delta"] * (th["nab_ad"] - nab) * G_SUN / p * w, 0.0))


# --------------------------------------------------------------------------
def analyse(entry):
    t_wall, path = entry
    term, fallbacks = console_ok(path)
    if term != "end_time":
        return dict(t_wall=t_wall, path=path, excluded=f"termination={term}")
    h_face, t, d = load_gamma_diag(path)
    oc = load_outercond(path)
    fh, ft, ff = load_face_flux(path + ".faceflux")
    ff = np.asarray(ff)
    w, h = geometry(h_face)                        # widths [m], cell centres [km]
    rho, V, T, p = d[:, :, I_RHO], d[:, :, I_V], d[:, :, I_T], d[:, :, I_P]
    x, ne, nhi, G1 = d[:, :, I_X], d[:, :, I_NE], d[:, :, I_NHI], d[:, :, I_G1]
    nt = len(t)

    # ---- global series ----------------------------------------------------
    F = rho[:, -1] * V[:, -1]                      # top cell; = top face to 0.1 %
    Ftf, Fbf = ff[:, -1], ff[:, 0]                 # top face, first interior face
    phi_c = np.concatenate([[0.0], np.cumsum(0.5 * (Ftf[1:] + Ftf[:-1]) * np.diff(ft))])
    Phi = np.interp(t, ft, phi_c)                  # evaporated mass [kg m^-2]
    M = (rho * w).sum(1)
    h75 = np.array([isotherm_height(h, T[k], T_MARK) for k in range(nt)])
    L = (h[-1] - h75) * 1.0e3                      # TR length, marker -> top cell [m]
    qw = np.interp(t, oc["t"], oc["q_face"]) if oc is not None else np.full(nt, np.nan)
    hen = 2.5 * p / rho + x * CHI / MH

    # ---- wave phase: F_base/F_top - 1 on the face captures -----------------
    sel = ft > 5.0
    tr, rr = ft[sel], Fbf[sel] / Ftf[sel] - 1.0
    rs = smooth(tr, rr, 40.0)                      # suppresses the hot-wall ripple
    ia = np.argmax(np.abs(Fbf[sel]) > 0.1 * np.abs(Ftf[sel]))
    t_arr = float(tr[ia])
    zc = [tr[i] - rs[i] * (tr[i + 1] - tr[i]) / (rs[i + 1] - rs[i])
          for i in range(ia, len(rs) - 1) if rs[i] * rs[i + 1] < 0 and tr[i] < 1000.0]
    period = 2.0 * float(np.mean(np.diff(zc))) if len(zc) >= 3 else np.nan
    peaks = []
    for a, b in zip(zc[:-1], zc[1:]):
        m = (tr >= a) & (tr <= b)
        if m.any():
            j = np.argmax(np.abs(rs[m]))
            peaks.append((tr[m][j], abs(rs[m][j])))
    peaks = [pk for pk in peaks if pk[1] > 0.02]
    tau_d = np.nan
    if len(peaks) >= 3:
        pt, pa = np.array(peaks).T
        tau_d = float(-1.0 / np.polyfit(pt, np.log(pa), 1)[0])

    def settle(eps):
        bad = np.where(np.abs(rs) > eps)[0]
        return float(tr[bad[-1] + 1]) if len(bad) and bad[-1] + 1 < len(tr) else np.nan

    # ---- chromosphere-top enthalpy, energy ratio, steady-TR invariant ------
    upflow = t_wall >= 15.0e3
    hc = np.full(nt, np.nan)
    ER = np.full(nt, np.nan)
    Jth = np.full(nt, np.nan)
    if upflow and oc is not None:
        for k in range(nt):
            if not np.isfinite(h75[k]) or t[k] < 20.0:
                continue
            q = physical_heat_flux(h, T[k], ne[k], nhi[k])
            j = int(np.searchsorted(h, h75[k]))
            j0 = j
            while j > 0 and abs(q[j]) > 0.01 * abs(qw[k]) and h[j0] - h[j] < 30.0:
                j -= 1
            if abs(q[j]) > 0.01 * abs(qw[k]):
                continue
            hc[k] = hen[k, j]
            ER[k] = F[k] * (hen[k, -1] - hc[k]) / qw[k]
            Jth[k] = j_theory(T[k, -1], p[k, -1], hc[k])
    FL_J = F * L / Jth

    # TR relaxation: ER(t) -> ER_inf with an exponential tail after the pulse
    tau_tr, er_inf, er_amp, t_tr5 = np.nan, np.nan, np.nan, np.nan
    if upflow:
        from scipy.optimize import curve_fit
        m = (t >= 300.0) & np.isfinite(ER)
        try:
            (er_inf, er_amp, tau_tr), _ = curve_fit(
                lambda tt, a, b, c: a + b * np.exp(-(tt - 300.0) / c),
                t[m], ER[m], p0=[0.93, 0.4, 500.0], maxfev=20000)
        except RuntimeError:
            pass
        if not (er_amp > 0.05):
            tau_tr = np.nan                           # already relaxed at 300 s
        dev = np.abs(ER - ER[-1]) > 0.05
        late = np.where(dev & (t > 20.0))[0]
        t_tr5 = float(t[late[-1] + 1]) if len(late) and late[-1] + 1 < nt else float(t[1])

    # ---- secular drift: contraction coefficient ----------------------------
    # K_th: linear coefficient of the initial column; K_cur(t): the same weight on
    # the current state, which carries the finite-displacement (nonlinear)
    # correction. mu = M_TR/(rho_f L): the lengthening TR must also absorb
    # chromospheric mass, which amplifies dL/dPhi to K/(1 - mu).
    K_th = float(contraction_weight(rho[0], T[0], p[0], w, h < (h75[0] - 15.0))[-1])
    K_cur = np.full(nt, np.nan)
    mu = np.full(nt, np.nan)
    for k in range(0, nt, 5):
        if not np.isfinite(h75[k]):
            continue
        K_cur[k] = contraction_weight(rho[k], T[k], p[k], w, h < (h75[k] - 5.0))[-1]
        above = h > h75[k]
        i_f = max(int(np.searchsorted(h, h75[k])) - 3, 0)
        mu[k] = (rho[k, above] * w[above]).sum() / (rho[k, i_f] * L[k])
    K_meas = r_meas = r_L = r_FL = r_L_pred = np.nan
    late = t >= t[-1] - 500.0
    kl = np.where(np.isfinite(K_cur) & late)[0]
    K_cur_end = float(np.nanmean(K_cur[kl])) if len(kl) else np.nan
    mu_end = float(np.nanmean(mu[kl])) if len(kl) else np.nan
    if upflow:
        t1 = max(1000.0, t_tr5 if np.isfinite(t_tr5) else 1000.0)
        m = (t >= t1)
        if m.sum() > 5:
            K_meas = float(np.polyfit(Phi[m], L[m], 1)[0])
        # -dlnF/dt = dlnL/dt - dln(F L)/dt over the last 500 s
        r_meas = float(-np.polyfit(t[late], np.log(F[late]), 1)[0])
        r_L = float(np.polyfit(t[late], np.log(L[late]), 1)[0])
        r_FL = float(-np.polyfit(t[late], np.log(F[late] * L[late]), 1)[0])
        r_L_pred = float(K_cur_end * F[-1] / ((1.0 - mu_end) * L[-1]))

    # ---- adiabatic parcel-lift test at the final frame ---------------------
    def m_above(k):
        c = np.cumsum((rho[k] * w)[::-1])[::-1]
        return c - 0.5 * rho[k] * w
    kf = nt - 1
    m0, mk = m_above(0), m_above(kf)
    lnT0, lnp0, lnr0 = np.log(T[0]), np.log(p[0]), np.log(rho[0])
    src = mk + Phi[kf]
    T_src = np.exp(np.interp(src, m0[::-1], lnT0[::-1]))
    p_src = np.exp(np.interp(src, m0[::-1], lnp0[::-1]))
    r_src = np.exp(np.interp(src, m0[::-1], lnr0[::-1]))
    T_pred = T_src * (p[kf] / p_src)**thermo(r_src, T_src)["nab_ad"]
    T_same_m = np.exp(np.interp(mk, m0[::-1], lnT0[::-1]))
    dT_meas, dT_pred = T[kf] - T_same_m, T_pred - T_same_m
    chrom_f = (T[kf] < 7000.0) & (h > 1605.0)
    big = chrom_f & (np.abs(dT_pred) > 1.0) & (h > 1620.0) & (h < h75[kf] - 10.0)
    adiab_err = (float(np.max(np.abs(dT_meas[big] - dT_pred[big]) / np.abs(dT_pred[big])))
                 if big.any() else np.nan)

    # ---- Eulerian flux profile at the final capture ------------------------
    jl = ft >= ft[-1] - 200.0                      # average out the hot-wall ripple
    prof_meas = np.mean(ff[jl] / ff[jl, -1:], axis=0)
    prof_pred = 1.0 - rho[kf] * contraction_weight(rho[kf], T[kf], p[kf], w,
                                                   np.ones_like(h, bool))
    fchrom = (T[kf] < 7000.0) & (fh > 1605.0) & (fh < h75[kf] - 10.0)
    prof_err = float(np.max(np.abs(prof_meas[fchrom] - prof_pred[fchrom])))
    i1, i2 = int(np.argmin(np.abs(fh - 1602.0))), int(np.argmin(np.abs(fh - 1850.0)))

    # ---- acoustic quantities of the initial column -------------------------
    c_eq = np.sqrt(G1[0] * p[0] / rho[0])
    c_fr = np.sqrt(5.0 / 3.0 * p[0] / rho[0])
    late = t >= 2000.0
    D_over_F = float(np.median((M[0] - M[late]) / F[late])) if upflow else np.nan

    return dict(
        t_wall=t_wall, path=path, termination=term, godunov_fallbacks=fallbacks,
        series=dict(t=t, F=F, Phi=Phi, M=M, L=L, h75=h75, ER=ER, J=Jth, FL_J=FL_J,
                    K_cur=K_cur, mu=mu,
                    qw=qw, ft=tr, ratio=rr, ratio_s=rs),
        profile=dict(face_km=fh, meas=prof_meas, pred=prof_pred, mask=fchrom,
                     h=h, dT_meas=dT_meas, dT_pred=dT_pred,
                     chrom=chrom_f & (h < h75[kf] - 10.0)),
        spacetime=dict(ft=ft, fh=fh, ff=ff, h=h, c_eq=c_eq, w=w) if t_wall in
        (T_RELEASE, 100.0e3) else None,
        metrics=dict(
            T_wall_K=t_wall, t_ac_eq_s=float(np.sum(w / c_eq)),
            t_ac_frozen_s=float(np.sum(w / c_fr)), t_arrival_base_s=t_arr,
            ring_period_s=period, ring_efold_s=tau_d,
            settle_10pct_s=settle(0.10), settle_5pct_s=settle(0.05),
            tau_TR_s=float(tau_tr) if np.isfinite(tau_tr) else None,
            ER_inf=float(er_inf) if np.isfinite(er_inf) else None,
            t_TR_within_5pct_s=t_tr5 if np.isfinite(t_tr5) else None,
            FL_over_J_4000=float(FL_J[-1]) if np.isfinite(FL_J[-1]) else None,
            K_theory_m3_per_kg=K_th, K_measured_m3_per_kg=K_meas,
            K_current_end_m3_per_kg=K_cur_end, mu_TR_mass_fraction_end=mu_end,
            K_eff_predicted_m3_per_kg=(K_cur_end / (1.0 - mu_end)
                                       if np.isfinite(mu_end) else np.nan),
            decay_rate_measured_per_s=r_meas, L_growth_rate_measured_per_s=r_L,
            FL_decay_rate_measured_per_s=r_FL, L_growth_rate_predicted_per_s=r_L_pred,
            F_top_4000=float(F[-1]), L_TR_4000_km=float(L[-1] / 1e3),
            J_theory_4000=float(Jth[-1]) if np.isfinite(Jth[-1]) else None,
            Phi_4000=float(Phi[-1]), mass_deficit_over_F_s=D_over_F,
            adiabatic_max_rel_err_4000=adiab_err,
            flux_profile_max_abs_err_4000=prof_err,
            flux_spread_1600_1850_meas=float(prof_meas[i1] - prof_meas[i2]),
            flux_spread_1600_1850_pred=float(prof_pred[i1] - prof_pred[i2])))


def release_profile(path):
    """Measured vs predicted F(z)/F_top at the last capture of one release run."""
    h_face, t, d = load_gamma_diag(path)
    fh, ft, ff = load_face_flux(path + ".faceflux")
    ff = np.asarray(ff)
    w, h = geometry(h_face)
    k = len(t) - 1
    jl = ft >= ft[-1] - 200.0
    rho, T, p = d[k, :, I_RHO], d[k, :, I_T], d[k, :, I_P]
    pred = 1.0 - rho * contraction_weight(rho, T, p, w, np.ones_like(h, bool))
    h75 = isotherm_height(h, T, T_MARK)
    mask = (T < 7000.0) & (fh > 1605.0) & (fh < h75 - 10.0)
    return dict(face_km=fh, meas=np.mean(ff[jl] / ff[jl, -1:], axis=0), pred=pred,
                mask=mask, t=float(t[k]))


# --------------------------------------------------------------------------
def plot_phases(runs, colors, rel, out):
    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5))
    ax = ax.ravel()
    wave_walls = [15e3, 22e3, 30e3, 45e3]
    hot = [r for r in runs if r["t_wall"] >= 15e3]

    # (a) wave phase
    a = ax[0]
    for r in runs:
        if r["t_wall"] not in wave_walls:
            continue
        s = r["series"]
        a.plot(s["ft"], s["ratio_s"], color=colors[r["t_wall"]], lw=1.4,
               label=f"{r['t_wall']/1e3:g} kK")
    t_ac = runs[0]["metrics"]["t_ac_eq_s"]
    a.axhspan(-0.05, 0.05, color="0.85", zorder=0)
    for n in range(1, 5):
        a.axvline(n * t_ac, color="0.6", lw=0.6, ls=":")
    a.set_xlim(0, 1500)
    a.set_ylim(-1.2, 2.0)
    a.set_xlabel("t [s]")
    a.set_ylabel(r"$F_{\rm base}/F_{\rm top}-1$")
    a.set_title(f"(a) wave phase: pulse bounces, $t_{{ac}}$={t_ac:.0f} s (dotted); "
                "grey = ±5 %", fontsize=10)
    a.legend(fontsize=8)

    # (b) TR relaxation: F*L/J -> 1
    a = ax[1]
    for r in hot:
        s = r["series"]
        m = np.isfinite(s["FL_J"])
        a.plot(s["t"][m], s["FL_J"][m], color=colors[r["t_wall"]],
               lw=2.2 if r["t_wall"] == T_RELEASE else 1.2, label=f"{r['t_wall']/1e3:g} kK")
    a.axhspan(0.95, 1.05, color="0.85", zorder=0)
    a.axhline(1.0, color="k", lw=0.5)
    a.set_xscale("log")
    a.set_xlim(40, 4000)
    a.set_ylim(0.8, 2.2)
    a.set_xlabel("t [s]")
    a.set_ylabel(r"$F\,L_{\rm TR}\,/\,J(T_{\rm top},p_{\rm top},h_c)$")
    a.set_title("(b) TR relaxation to the steady conductive-enthalpy balance", fontsize=10)
    a.legend(fontsize=7, ncol=2)

    # (c) drift of F vs the invariant F*L
    a = ax[2]
    for r in hot:
        if r["t_wall"] not in (22e3, 45e3, 100e3):
            continue
        s = r["series"]
        k1 = int(np.argmin(np.abs(s["t"] - 1000.0)))
        c = colors[r["t_wall"]]
        a.plot(s["t"], s["F"] / s["F"][k1], color=c, lw=1.6,
               label=f"{r['t_wall']/1e3:g} kK: $F$")
        FL = s["F"] * s["L"]
        a.plot(s["t"], FL / FL[k1], color=c, lw=1.6, ls="--",
               label=f"{r['t_wall']/1e3:g} kK: $F\\,L_{{\\rm TR}}$")
    a.set_xscale("log")
    a.set_xlim(40, 4000)
    a.set_xlabel("t [s]")
    a.set_ylabel("value / value at 1000 s")
    a.set_title(r"(c) $F$ keeps falling; $F\,L_{\rm TR}$ is invariant", fontsize=10)
    a.legend(fontsize=7, ncol=2)

    # (d) contraction law
    a = ax[3]
    K = np.median([r["metrics"]["K_theory_m3_per_kg"] for r in hot])
    xmax = 0.0
    for r in hot:
        if r["t_wall"] < 22e3:
            continue
        s = r["series"]
        k1 = int(np.argmin(np.abs(s["t"] - 300.0)))
        dP, dL = s["Phi"][k1:] - s["Phi"][k1], (s["L"][k1:] - s["L"][k1]) / 1e3
        a.plot(dP, dL, color=colors[r["t_wall"]], lw=1.6, label=f"{r['t_wall']/1e3:g} kK")
        xmax = max(xmax, dP[-1])
    xx = np.linspace(0, xmax, 10)
    a.plot(xx, K * xx / 1e3, "k--", lw=1.0,
           label=rf"theory $K$={K:.2e} m$^3$/kg")
    a.set_xlabel(r"evaporated mass $\Phi(t)-\Phi(300\,{\rm s})$ [kg m$^{-2}$]")
    a.set_ylabel(r"$L_{\rm TR}(t)-L_{\rm TR}(300\,{\rm s})$ [km]")
    a.set_title(r"(d) TR lengthens as the chromosphere contracts: $dL/d\Phi\approx K$",
                fontsize=10)
    a.legend(fontsize=8)

    # (e) adiabatic parcel-lift test
    a = ax[4]
    for r in hot:
        if r["t_wall"] not in (22e3, 45e3, 100e3):
            continue
        pr = r["profile"]
        m = pr["chrom"]
        c = colors[r["t_wall"]]
        a.plot(pr["h"][m], pr["dT_meas"][m], color=c, lw=1.8,
               label=f"{r['t_wall']/1e3:g} kK measured")
        a.plot(pr["h"][m][::12], pr["dT_pred"][m][::12], "o", ms=4, mfc="none", color=c,
               label=f"{r['t_wall']/1e3:g} kK adiabat")
    a.set_xlabel("height at 4000 s [km]")
    a.set_ylabel(r"$T(m,4000\,{\rm s})-T(m,0)$ at fixed column mass [K]")
    a.set_title("(e) chromosphere cools exactly adiabatically as it is lifted", fontsize=10)
    a.legend(fontsize=7, ncol=2)

    # (f) Eulerian flux profile
    a = ax[5]
    for r in runs:
        if r["t_wall"] not in (15e3, 45e3, 70e3):
            continue
        pr = r["profile"]
        m = pr["mask"]
        c = colors[r["t_wall"]]
        a.plot(pr["face_km"][m], pr["meas"][m], color=c, lw=1.6,
               label=f"{r['t_wall']/1e3:g} kK measured")
        a.plot(pr["face_km"][m], pr["pred"][m], color=c, lw=1.0, ls="--")
    if rel is not None:
        m = rel["mask"]
        a.plot(rel["face_km"][m], rel["meas"][m], color="k", lw=2.2,
               label="release (single T) measured")
        a.plot(rel["face_km"][m], rel["pred"][m], color="0.5", lw=1.2, ls="--",
               label="dashed = contraction prediction")
    a.set_xlabel("face height [km]")
    a.set_ylabel(r"$f_{\rm total}(z)\,/\,f_{\rm top}$, mean of 3800-4000 s")
    a.set_title("(f) the 12-13 % flux spread is the contraction's Eulerian signature",
                fontsize=10)
    a.legend(fontsize=7)
    for a in ax:
        a.grid(alpha=0.25)
    fig.suptitle("model_column T_wall sweep (N=500/R4, 4000 s): the three stages of the "
                 "evaporation mass flux", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=130)
    plt.close(fig)


def plot_spacetime(runs, out):
    sel = [r for r in runs if r["spacetime"] is not None]
    fig, ax = plt.subplots(1, len(sel), figsize=(7.5 * len(sel), 6), squeeze=False)
    for a, r in zip(ax[0], sel):
        st = r["spacetime"]
        m = (st["ft"] > 5.0) & (st["ft"] <= 900.0)
        z = st["ff"][m] / np.abs(st["ff"][m, -1])[:, None]
        im = a.pcolormesh(st["fh"], st["ft"][m], np.clip(z, -1, 3), cmap="RdBu_r",
                          vmin=-1, vmax=3, shading="auto", rasterized=True)
        tt = np.cumsum((st["w"] / st["c_eq"])[::-1])[::-1]   # travel time to the top
        for n in range(8):
            a.plot(st["h"], n * tt[0] + (tt if n % 2 == 0 else tt[0] - tt),
                   "k--", lw=0.6)
        a.set_ylim(0, 900)
        a.set_xlabel("face height [km]")
        a.set_ylabel("t [s]")
        a.set_title(f"{r['t_wall']/1e3:g} kK: $f_{{\\rm total}}(z,t)/f_{{\\rm top}}(t)$; "
                    r"dashed = $\pm c_{\rm eq}$ characteristics from the top at $t=0$",
                    fontsize=10)
        fig.colorbar(im, ax=a, label=r"$f/f_{\rm top}$ (clipped to [-1, 3])")
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


def plot_timescales(runs, colors, out):
    hot = [r for r in runs if r["t_wall"] >= 15e3]
    Tw = np.array([r["t_wall"] for r in hot]) / 1e3
    g = lambda k: np.array([np.nan if r["metrics"][k] is None else r["metrics"][k]
                            for r in hot], dtype=float)
    fig, a = plt.subplots(figsize=(8.5, 7.2))
    t_ac = hot[0]["metrics"]["t_ac_eq_s"]
    a.axhline(t_ac, color="0.5", lw=1.0, ls=":", label=f"one-way acoustic time $t_{{ac}}$ = {t_ac:.0f} s")
    a.plot(Tw, g("settle_5pct_s"), "s-", color="tab:blue", lw=1.2,
           label=r"wave phase over: $|F_{\rm base}/F_{\rm top}-1|<5\,\%$")
    a.plot(Tw, g("t_TR_within_5pct_s"), "o-", color="tab:orange", lw=1.2,
           label="TR relaxed: energy ratio within 0.05 of its 4000 s value")
    a.plot(Tw, g("tau_TR_s"), "o", mfc="none", color="tab:orange", ms=7,
           label=r"TR relaxation e-folding $\tau_{\rm TR}$ (fit)")
    a.plot(Tw, 1.0 / g("decay_rate_measured_per_s"), "D-", color="tab:red", lw=1.2,
           label=r"flux-drift e-folding $1/r$ at 4000 s, measured")
    a.plot(Tw, 1.0 / g("L_growth_rate_predicted_per_s"), "D", mfc="none", color="tab:red",
           ms=7, label=r"contraction only: $(1-\mu)L/(K F)$")
    a.axhline(4000.0, color="k", lw=0.6, ls="--", label="run length 4000 s")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xticks(Tw)
    a.set_xticklabels([f"{v:g}" for v in Tw])
    a.minorticks_off()
    a.set_xlabel(r"$T_{\rm wall}$ [kK]")
    a.set_ylabel("time [s]")
    a.set_title("When the evaporation flux is ready, and how fast it then drifts",
                fontsize=11)
    a.grid(alpha=0.25, which="both")
    a.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2,
             frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="outputs/model_column/twall_sweep")
    ap.add_argument("--cfl", default="0.50")
    ap.add_argument("--run", action="append", default=[],
                    help="explicit T_wall_K=path entry, added to the glob")
    ap.add_argument("--release", default="outputs/model_column/release_singleT_N500_4000s.txt",
                    help="release run whose flux profile is overlaid ('' to skip)")
    ap.add_argument("--out-dir", default="visualization/model_column")
    ap.add_argument("--prefix", default="twall_relaxation")
    ap.add_argument("--json", default="outputs/model_column/twall_sweep/"
                                      "twall_relaxation_metrics.json")
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()

    entries = []
    for path in sorted(glob.glob(os.path.join(args.dir, "twall_*kK_cfl*.txt"))):
        mm = re.search(r"twall_(\d+)kK_cfl([0-9.]+)\.txt$", os.path.basename(path))
        if mm and (args.cfl == "all" or mm.group(2) == args.cfl):
            entries.append((float(mm.group(1)) * 1.0e3, path))
    for spec in args.run:
        tw, path = spec.split("=", 1)
        entries.append((float(tw), path))
    entries.sort()
    with ProcessPoolExecutor(args.workers) as ex:
        results = list(ex.map(analyse, entries))
    runs = []
    for r in results:
        if "excluded" in r:
            print(f"{r['t_wall']/1e3:>5g} kK  EXCLUDED: {r['excluded']}")
            continue
        runs.append(r)
        m = r["metrics"]
        fmt = lambda v, f: "   -   " if v is None or not np.isfinite(v) else format(v, f)
        print(f"{r['t_wall']/1e3:>5g} kK  P={fmt(m['ring_period_s'], '5.0f')} s  "
              f"tau_d={fmt(m['ring_efold_s'], '5.0f')} s  settle5={fmt(m['settle_5pct_s'], '5.0f')} s  "
              f"tau_TR={fmt(m['tau_TR_s'], '5.0f')} s  t_TR5={fmt(m['t_TR_within_5pct_s'], '5.0f')} s  "
              f"FL/J={fmt(m['FL_over_J_4000'], '.3f')}  K_IC={m['K_theory_m3_per_kg']:.2e} "
              f"K_cur={fmt(m['K_current_end_m3_per_kg'], '.2e')} mu={fmt(m['mu_TR_mass_fraction_end'], '.2f')} "
              f"K_eff={fmt(m['K_eff_predicted_m3_per_kg'], '.2e')}/meas {fmt(m['K_measured_m3_per_kg'], '.2e')}  "
              f"r={fmt(m['decay_rate_measured_per_s'], '.2e')} = L {fmt(m['L_growth_rate_measured_per_s'], '.2e')}"
              f"(pred {fmt(m['L_growth_rate_predicted_per_s'], '.2e')}) + FL {fmt(m['FL_decay_rate_measured_per_s'], '.2e')}  "
              f"adiab={fmt(m['adiabatic_max_rel_err_4000'], '.3f')}  "
              f"prof={m['flux_profile_max_abs_err_4000']:.3f}  "
              f"spread={m['flux_spread_1600_1850_meas']:.3f}/{m['flux_spread_1600_1850_pred']:.3f}")
    rel = release_profile(args.release) if args.release else None
    if rel is not None:
        m = rel["mask"]
        i1 = int(np.argmin(np.abs(rel["face_km"] - 1602.0)))
        i2 = int(np.argmin(np.abs(rel["face_km"] - 1850.0)))
        print(f"release {args.release} t={rel['t']:.0f} s: spread 1600-1850 "
              f"{rel['meas'][i1]-rel['meas'][i2]:.3f} measured / "
              f"{rel['pred'][i1]-rel['pred'][i2]:.3f} predicted; max |err| "
              f"{np.max(np.abs(rel['meas'][m]-rel['pred'][m])):.4f}")

    colors, _, _ = wall_colors([r["t_wall"] for r in runs])
    os.makedirs(args.out_dir, exist_ok=True)
    p = lambda n: os.path.join(args.out_dir, f"{args.prefix}_{n}")
    plot_phases(runs, colors, rel, p("phases.png"))
    plot_spacetime(runs, p("spacetime.png"))
    plot_timescales(runs, colors, p("timescales.png"))
    if args.json:
        doc = dict(runs=[r["metrics"] for r in runs])
        if rel is not None:
            m = rel["mask"]
            doc["release_profile"] = dict(
                path=args.release, t=rel["t"],
                max_abs_err=float(np.max(np.abs(rel["meas"][m] - rel["pred"][m]))))
        with open(args.json, "w") as fh:
            json.dump(doc, fh, indent=1, default=float)
        print(f"wrote {args.json}")
    for n in ("phases", "spacetime", "timescales"):
        print(f"wrote {p(n + '.png')}")


if __name__ == "__main__":
    main()
