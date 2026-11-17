#!/usr/bin/env python3
"""Figures and numbers for the 2026-10-27 advisor deck (outer pressure BC; TR pulses).

Part 1 — outer boundary. Compares, at T_wall = 22 kK (4000 s) and 100 kK (first
1500 s), two outer closures that differ only in the ghost pressure:
  fixed p      fixed reservoir back-pressure, hydro ghost T pinned to the wall
               (ISO_HYDRO_T_DECOUPLE=1 ISO_OUTER_T_HYDRO_WALL=1, equivalent to the
               release default ISO_HYDRO_T_DECOUPLE=0)
  floating p   the same plus ISO_OUTER_P_NEUMANN=2: ghost pressure extrapolated
               hydrostatically off the live top cell
Writes outer_bc_rerun_timeseries.pdf: top-cell velocity, top-cell pressure and
column mass against time.

Part 2 — grid-locked pulses in the transition region (TR) at hot walls. From the
face-flux sidecars (cell state every 1780 steps):
  tr_pulses_locking.pdf      100 kK, N=500: TR velocity minus its 60 s running
                             median, with the times at which the 7 kK isotherm
                             (the TR foot) enters a new cell;
  tr_pulses_resolution.pdf   100 kK, N=500 against N=1000 (every cell halved),
                             side by side: T(h) at the foot with the cell centres
                             marked, the TR velocity pulses and the top velocity;
  tr_pulses_face_conductivity.pdf  100 kK, N=500, three face-conductivity rules:
                             default (harmonic mean), ISO_COND_KIRCHHOFF=1 and
                             ISO_COND_ARITH=1;
plus tr_pulses_metrics.json with the pulse amplitude, foot width and mean
top-face mass flux of every run, including the 22/45/70/100 kK sweep runs.

Usage (from the project root; the runs live on the mac-mini):
    .venv/bin/python util/release/plot_outer_bc_and_tr_pulses.py [--out-dir visualization/model_column]
"""
import argparse
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from twall_sweep_diag import COLUMNS, load_gamma_diag

ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = "outputs/model_column"
SWEEP = f"{OUT}/twall_sweep_hse_recon"
SWAP = f"{OUT}/outer_bc_swap_hse_recon"
IV, IR, IP, IT = (COLUMNS.index(c) for c in ("v", "rho", "p", "T"))
T_FOOT = 7000.0          # TR-foot marker isotherm [K]
HP_WINDOW_S = 60.0       # running-median window that removes the slow evolution [s]

BC_RUNS = {
    "22 kK": [("fixed $p$ (today)", f"{SWAP}/tpinned/twall_22kK_cfl0.50.txt", "k", "-"),
              ("floating $p$ (your formulas)", f"{SWAP}/pneumann2/twall_22kK_cfl0.50.txt", "C3", "-")],
    "100 kK": [("fixed $p$ (today)", f"{SWAP}/tpinned/twall_100kK_cfl0.25.txt", "k", "-"),
               ("floating $p$ (your formulas)", f"{SWAP}/pneumann2/twall_100kK_cfl0.25.txt", "C3", "-")],
}
BC_TMAX = {"22 kK": 4000.0, "100 kK": 1500.0}

PULSE_RUNS = {
    "N500": f"{SWEEP}/twall_100kK_cfl0.25.txt",
    "N1000": f"{OUT}/hot_wall_resolution/n1000/twall_100kK_cfl0.25.txt",
    "N500 Kirchhoff": f"{OUT}/cond_kirchhoff/twall_100kK_cfl0.25.txt",
    "N500 arithmetic": f"{OUT}/cond_arith/twall_100kK_cfl0.25.txt",
}
SWEEP_WALLS = {22: "twall_22kK_cfl0.50", 45: "twall_45kK_cfl0.50",
               70: "twall_70kK_cfl0.50", 100: "twall_100kK_cfl0.25"}


def p(rel):
    return os.path.join(ROOT, rel)


# ---------------------------------------------------------------- part 1
def column_mass(h_km, rho):
    """Column mass per unit area [kg m^-2] for each frame (midpoint cell widths)."""
    faces = np.concatenate([[h_km[0] - 0.5 * (h_km[1] - h_km[0])],
                            0.5 * (h_km[1:] + h_km[:-1]),
                            [h_km[-1] + 0.5 * (h_km[-1] - h_km[-2])]])
    return (rho * np.diff(faces) * 1e3).sum(axis=1)


def plot_bc(out_dir):
    fig, ax = plt.subplots(2, 3, figsize=(13, 6.6), constrained_layout=True)
    for r, (wall, runs) in enumerate(BC_RUNS.items()):
        for label, rel, color, ls in runs:
            h, t, d = load_gamma_diag(p(rel))
            h = np.asarray(h)
            keep = t <= BC_TMAX[wall] + 1e-6
            t, d = t[keep], d[keep]
            m = column_mass(h, d[:, :, IR])
            ax[r, 0].plot(t[1:], d[1:, -1, IV], color=color, ls=ls, lw=1.2, label=label)
            ax[r, 1].plot(t, d[:, -1, IP] / d[0, -1, IP], color=color, ls=ls, lw=1.2)
            ax[r, 2].plot(t, 100 * m / m[0], color=color, ls=ls, lw=1.2)
        ax[r, 0].set(yscale="log", ylabel=f"$T_\\mathrm{{wall}}$ = {wall}\n"
                     "top-cell velocity [m/s]")
        ax[r, 1].set(ylabel=r"top-cell pressure / initial")
        ax[r, 2].set(ylabel="column mass [% of initial]")
        for a in ax[r]:
            a.set_xlim(0, BC_TMAX[wall])
    for a in ax[1]:
        a.set_xlabel("t [s]")
    ax[0, 0].legend(fontsize=8, loc="lower left")
    out = os.path.join(out_dir, "outer_bc_rerun_timeseries.pdf")
    fig.savefig(out)
    print("wrote", out)


# ---------------------------------------------------------------- part 2
FACE_COLS = dict(h=0, v=3, T=4, f=16)   # cell_km, v_cell, T_cell, f_total


def load_cells(rel, t0, t1):
    """Cell height, velocity, temperature and upper-face mass flux from a .faceflux."""
    want = sorted(FACE_COLS.values())
    times, frames, cur, rows = [], [], None, []
    with open(p(rel) + ".faceflux") as fh:
        for line in fh:
            if line.startswith("# t ="):
                if cur is not None and rows and t0 <= cur <= t1:
                    times.append(cur)
                    frames.append(np.array(rows))
                cur = float(line.split("=")[1].split()[0])
                rows = []
                if cur > t1:
                    break
                continue
            if line[0] == "#" or not line.strip() or cur is None or not t0 <= cur <= t1:
                continue
            sp = line.split()
            rows.append([float(sp[i]) for i in want])
    if cur is not None and t0 <= cur <= t1 and frames and len(rows) == len(frames[0]):
        times.append(cur)
        frames.append(np.array(rows))
    a = np.array(frames)
    out = {k: a[:, :, want.index(c)] for k, c in FACE_COLS.items()}
    out["h"] = out["h"][0]
    out["t"] = np.array(times)
    return out


def foot_index(T):
    return np.array([np.argmax(row > T_FOOT) for row in T])


def foot_height(h, row, T_mark):
    k = np.argmax(row > T_mark)
    return np.interp(T_mark, row[k - 1:k + 1], h[k - 1:k + 1])


def pulse_metrics(c):
    t, h, V, T, F = c["t"], c["h"], c["v"], c["T"], c["f"]
    dt = float(np.median(np.diff(t)))
    win = max(3, int(round(HP_WINDOW_S / dt)) | 1)
    vhp = V - median_filter(V, size=(win, 1), mode="nearest")
    foot = np.array([foot_height(h, row, T_FOOT) for row in T])
    tr = h[None, :] >= foot[:, None] + 1.0        # TR: more than 1 km above the foot
    rms = float(np.sqrt((vhp ** 2 * tr).sum() / tr.sum()))
    mean_v = float((V * tr).sum() / tr.sum())
    width = [foot_height(h, row, 7500.0) - foot_height(h, row, 7000.0) for row in T]
    i7 = foot_index(T)
    return dict(window_s=[float(t[0]), float(t[-1])], capture_dt_s=dt,
                tr_pulse_rms_m_s=rms, tr_mean_v_m_s=mean_v,
                tr_pulse_rms_rel=rms / abs(mean_v),
                top_v_hp_std_m_s=float(vhp[:, -1].std()),
                mean_top_face_flux=float(F[:, -1].mean()),
                foot_width_km_median=float(np.median(width)),
                cell_km_at_foot=float(np.median(np.diff(h)[i7])),
                foot_cell_crossings=int(np.count_nonzero(np.diff(i7)))), vhp


def plot_locking(out_dir, t0=2000.0, t1=2600.0):
    c = load_cells(PULSE_RUNS["N500"], t0, t1)
    _, vhp = pulse_metrics(c)
    i7 = foot_index(c["T"])
    cross = c["t"][1:][np.diff(i7) != 0]
    fig, ax = plt.subplots(2, 1, figsize=(11, 5.6), sharex=True, constrained_layout=True,
                           gridspec_kw=dict(height_ratios=[1.6, 1]))
    for dz, col in ((2.0, "C1"), (8.0, "C2"), (None, "C3")):
        k = len(c["h"]) - 1 if dz is None else int(np.searchsorted(c["h"], c["h"][i7.max()] + dz))
        lab = "top cell" if dz is None else f"{dz:.0f} km above the 7 kK isotherm"
        ax[0].plot(c["t"], vhp[:, k], color=col, lw=0.9, label=lab)
    for x in cross:
        for a in ax:
            a.axvline(x, color="k", lw=0.5, alpha=0.35)
    ax[0].set(ylabel="velocity minus 60 s\nrunning median [m/s]")
    ax[0].legend(fontsize=8, loc="upper right")
    ax[1].plot(c["t"], c["h"][i7], drawstyle="steps-post", color="k")
    ax[1].set(xlabel="t [s]", ylabel="cell holding the\n7 kK isotherm [km]")
    out = os.path.join(out_dir, "tr_pulses_locking.pdf")
    fig.savefig(out)
    print("wrote", out)


def plot_pair(out_dir, name, labels, t0, t1, metrics, titles=None):
    cases = [(lab, load_cells(PULSE_RUNS[lab], t0, t1)) for lab in labels]
    titles = titles or labels
    n = len(cases)
    fig, ax = plt.subplots(3, n, figsize=(6 * n, 8.2), constrained_layout=True,
                           sharex="row", sharey="row")
    for j, (lab, c) in enumerate(cases):
        m, vhp = pulse_metrics(c)
        metrics[f"100 kK {lab} ({t0:.0f}-{t1:.0f} s)"] = m
        h, T, t = c["h"], c["T"], c["t"]
        foot = foot_height(h, T[-1], T_FOOT)
        sel = (h > foot - 1.5) & (h < foot + 2.5)
        ax[0, j].plot(h[sel], T[-1, sel], "o-", ms=3.5, lw=1.0, color="k")
        ax[0, j].axhspan(7000, 7500, color="C1", alpha=0.25, lw=0)
        ax[0, j].set(yscale="log", xlabel="height [km]",
                     title=f"{titles[j]}\ncells {m['cell_km_at_foot'] * 1e3:.0f} m")
        feet = [foot_height(h, row, T_FOOT) for row in T]       # follows the sinking base
        ks = [int(np.searchsorted(h, f + 3.0)) for f in feet]
        ax[1, j].plot(t, vhp[np.arange(len(ks)), ks], lw=0.8, color="C1")
        ax[1, j].set(title=f"pulse rms over the transition region: {m['tr_pulse_rms_m_s']:.0f} m/s "
                           f"({100 * m['tr_pulse_rms_rel']:.1f} %)")
        ax[2, j].plot(t, c["v"][:, -1], lw=0.8, color="C3")
        ax[2, j].set(xlabel="t [s]",
                     title=f"(mean top mass flux {m['mean_top_face_flux']:.3e} kg m$^{{-2}}$ s$^{{-1}}$)")
    ax[0, 0].set_ylabel(f"T [K] at t = {t1:.0f} s\n(dots = cell centres)")
    ax[1, 0].set_ylabel("pulse 3 km above the\n7 kK isotherm [m/s]")
    ax[2, 0].set_ylabel("top-cell velocity [m/s]")
    out = os.path.join(out_dir, name)
    fig.savefig(out)
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out-dir", default=p("visualization/model_column"))
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    plot_bc(a.out_dir)
    plot_locking(a.out_dir)
    metrics = {}
    for wall, tag in SWEEP_WALLS.items():
        m, _ = pulse_metrics(load_cells(f"{SWEEP}/{tag}.txt", 300.0, 700.0))
        metrics[f"{wall} kK N500 (300-700 s)"] = m
    plot_pair(a.out_dir, "tr_pulses_resolution.pdf", ["N500", "N1000"], 300.0, 500.0, metrics,
              titles=["base mesh (N500)", "every cell halved (N1000)"])
    plot_pair(a.out_dir, "tr_pulses_face_conductivity.pdf",
              ["N500", "N500 Kirchhoff", "N500 arithmetic"], 300.0, 700.0, metrics,
              titles=["harmonic mean (default)", "Kirchhoff, Spitzer part only",
                      "arithmetic mean"])
    out = p(f"{OUT}/hot_wall_resolution/tr_pulses_metrics.json")
    json.dump(metrics, open(out, "w"), indent=1)
    for k, m in metrics.items():
        print(f"{k:28s} rms {m['tr_pulse_rms_m_s']:7.2f} m/s ({100 * m['tr_pulse_rms_rel']:5.2f} %)  "
              f"top hp std {m['top_v_hp_std_m_s']:7.2f}  F_top {m['mean_top_face_flux']:.4e}  "
              f"foot {1e3 * m['foot_width_km_median']:6.0f} m / cell {1e3 * m['cell_km_at_foot']:4.0f} m  "
              f"crossings {m['foot_cell_crossings']}")
    print("wrote", out)


if __name__ == "__main__":
    main()
