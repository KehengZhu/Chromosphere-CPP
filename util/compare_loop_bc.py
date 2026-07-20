#!/usr/bin/env python3
"""Half-loop (reflecting apex) vs full-loop (apex interior) — the apex-BC timing fix.

The half-loop closed model puts a reflecting symmetry plane at the apex: the
evaporation jet from the (single) flaring footpoint bounces straight back. That is
exact ONLY for a symmetric flare. For one-sided heating — the realistic ribbon case
— the jet should stream THROUGH the apex, down the far leg, reflect off footpoint B,
and return only after a full far-leg round trip. This script shows the difference.

Plotted vs ARC LENGTH s (cumulative ds), NOT height: on the full loop height folds
the two legs onto each other, whereas arc length lays the loop out flat
(footpoint A at s=0 | apex at s=L | footpoint B at s=2L).

Usage:
  python util/compare_loop_bc.py            # uses the default half/full outputs
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

m_i = 1.6726219e-27
k_b = 1.380649e-23
g_si = 274.0


def parse_dat(path):
    rows, in_cells = [], False
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
    ds = a[:, 1]
    phi_c = 0.5 * (a[:, 4] + a[:, 5])
    s_ctr = np.cumsum(ds) - 0.5 * ds          # cell-center arc length [m]
    return ds, s_ctr, phi_c, phi_c / g_si     # ds, arc, phi_g, height[m]


def read_frames(path):
    raw = [l for l in open(path).read().splitlines() if l.strip()]
    ns, neq = map(int, raw[0].split())
    i, frames = 2, []
    while i < len(raw):
        t = float(raw[i].split()[3]); i += 1
        d = np.zeros((ns, neq)); ok = True
        for j in range(ns):
            tk = raw[i + j].split()
            if len(tk) != neq:
                ok = False; break
            d[j] = [float(x) for x in tk]
        i += ns
        if ok:
            frames.append((t, d))
    return frames


def stack_V(frames):
    ts, V = [], []
    for t, xn in frames:
        ts.append(t)
        V.append(xn[:, 2] / xn[:, 0])         # ion velocity [m/s], + up the loop
    return np.array(ts), np.array(V)


def main():
    half_dat = "scenarios/data/loop_closed_fine.dat"
    full_dat = "scenarios/data/loop_full_fine.dat"
    half_out = "outputs/loops/loop_half_flare.txt"
    full_out = "outputs/loops/loop_full_flare.txt"
    out = "visualization/loops/loop_bc_compare.png"
    os.makedirs(os.path.dirname(out), exist_ok=True)

    dsh, sh, phih, _ = parse_dat(half_dat)
    dsf, sfu, phif, _ = parse_dat(full_dat)
    th, Vh = stack_V(read_frames(half_out))
    tf, Vf = stack_V(read_frames(full_out))

    shMm, sfMm = sh / 1e6, sfu / 1e6
    apex_h = int(np.argmax(phih))             # half: apex = top cell (= last)
    apex_f = int(np.argmax(phif))             # full: apex = mid cell
    Lh = shMm[apex_h]                         # half-loop apex arc [Mm]
    Lf_apex = sfMm[apex_f]                    # full-loop apex arc [Mm]
    fpB = sfMm[-1]                            # full-loop footpoint B arc [Mm]

    tmax = min(th[-1], tf[-1])
    # leg transit time estimate from the peak upflow (jet) speed near the footpoint
    vjet = np.nanpercentile(np.abs(Vf), 99.5)
    t_leg = (Lf_apex * 1e6) / vjet
    print(f"half: ns={Vh.shape[1]} t<= {th[-1]:.0f}s apex@s={Lh:.1f}Mm(cell {apex_h})")
    print(f"full: ns={Vf.shape[1]} t<= {tf[-1]:.0f}s apex@s={Lf_apex:.1f}Mm(cell {apex_f}) "
          f"fpB@s={fpB:.1f}Mm")
    print(f"jet speed ~{vjet/1e3:.0f} km/s -> one-leg transit ~{t_leg:.0f}s; "
          f"half-loop return is too early by ~2x that for slow flows")

    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], hspace=0.28, wspace=0.22)
    fig.suptitle("Apex boundary condition: half-loop (reflecting apex) vs full-loop "
                 "(apex is interior) — one-sided beam", fontsize=13)

    vlim = np.nanpercentile(np.abs(np.concatenate([Vh.ravel(), Vf.ravel()])), 99) / 1e3
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-vlim, vmax=vlim)

    # --- space-time V, half loop -----------------------------------------
    axh = fig.add_subplot(gs[0, 0])
    im = axh.pcolormesh(th, shMm, Vh.T / 1e3, shading="auto", cmap="RdBu_r", norm=norm)
    axh.axhline(Lh, color="k", ls="--", lw=1.2)
    axh.text(0.98 * tmax, Lh, " reflecting apex (wall)", va="bottom", ha="right",
             fontsize=9, color="k")
    axh.text(0.02 * tmax, 0.3, "footpoint A (beam)", color="k", fontsize=9)
    axh.set_title("HALF loop: jet bounces at the apex wall")
    axh.set_ylabel("arc length s [Mm]"); axh.set_xlabel("t [s]")
    axh.set_xlim(0, tmax); fig.colorbar(im, ax=axh, label="V [km/s] (+ up the loop)")

    # --- space-time V, full loop -----------------------------------------
    axf = fig.add_subplot(gs[0, 1])
    im = axf.pcolormesh(tf, sfMm, Vf.T / 1e3, shading="auto", cmap="RdBu_r", norm=norm)
    axf.axhline(Lf_apex, color="k", ls="--", lw=1.2)
    axf.text(0.98 * tmax, Lf_apex, " apex (interior)", va="bottom", ha="right",
             fontsize=9, color="k")
    axf.text(0.02 * tmax, 0.3, "footpoint A (beam)", color="k", fontsize=9)
    axf.text(0.02 * tmax, fpB - 1.5, "footpoint B (quiet)", color="k", fontsize=9)
    axf.set_title("FULL loop: jet streams through apex to far leg")
    axf.set_ylabel("arc length s [Mm]"); axf.set_xlabel("t [s]")
    axf.set_xlim(0, tmax); fig.colorbar(im, ax=axf, label="V [km/s] (+ up the loop)")

    # --- apex velocity vs time -------------------------------------------
    axv = fig.add_subplot(gs[1, 0])
    axv.plot(th, Vh[:, apex_h] / 1e3, color="0.5", lw=1.6,
             label="half loop (apex = reflecting wall)")
    axv.plot(tf, Vf[:, apex_f] / 1e3, color="C3", lw=1.6,
             label="full loop (apex = interior cell)")
    axv.axhline(0, color="k", lw=0.6)
    axv.set_xlim(0, tmax); axv.set_xlabel("t [s]")
    axv.set_ylabel("V at apex [km/s]"); axv.grid(alpha=0.3)
    axv.set_title("Apex velocity: forced to 0 by the wall vs free jet through-flow")
    axv.legend(fontsize=9, loc="best")

    # --- far-footpoint arrival (full) ------------------------------------
    axb = fig.add_subplot(gs[1, 1])
    # near-footpoint-B downflow (last few % of arc), full loop only
    nb = max(3, int(0.02 * Vf.shape[1]))
    vB = np.nanmean(Vf[:, -nb:], axis=1) / 1e3
    vA = np.nanmean(Vf[:, :nb], axis=1) / 1e3
    axb.plot(tf, vA, color="C0", lw=1.5, label="footpoint A (flaring)")
    axb.plot(tf, vB, color="C1", lw=1.5, label="footpoint B (quiet, far leg)")
    axb.axhline(0, color="k", lw=0.6)
    axb.axvline(t_leg, color="0.6", ls=":", lw=1.2)
    axb.text(t_leg, axb.get_ylim()[1], f" 1 leg transit ~{t_leg:.0f}s",
             va="top", ha="left", fontsize=8, color="0.4")
    axb.set_xlim(0, tmax); axb.set_xlabel("t [s]")
    axb.set_ylabel("mean V near footpoint [km/s]"); axb.grid(alpha=0.3)
    axb.set_title("Full loop predicts a DELAYED disturbance at the far footpoint")
    axb.legend(fontsize=9, loc="best")

    fig.savefig(out, dpi=130, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
