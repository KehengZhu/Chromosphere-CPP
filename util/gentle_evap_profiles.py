#!/usr/bin/env python3
"""Preflare-vs-evaporation field-aligned profiles for a gentle-evaporation run.

Overlays the relaxed PREFLARE profile (just before the heating ramp) and a
DURING/POST-evaporation profile (default: the last frame) of temperature, density,
and field-aligned velocity, on the split x-axis the project standardizes for
chromosphere→corona loops (zoom the sub-Mm chromosphere+TR, compress the extended
corona; region bands marked). Shows, in space, what the EM-rise time series
(gentle_evaporation_diag.py) shows in time: the upper chromosphere heated toward
coronal temperatures, the TR pushed down, the corona filled (n_e up), and a
sustained subsonic upflow. Reuses the split-axis transform from animate_loop_flare.py.

Usage:
  python util/gentle_evap_profiles.py scenarios/data/loop_closed_gentle.dat \
      outputs/gentle_evap.txt [visualization/gentle_evap_profiles.png] [t_pre]
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from animate_loop_flare import (read_meta, parse_dat, read_frames, primitives,
                                build_transform, region_layout)

T_MIN, T_MAX = 3.0e3, 5.0e7
N_MIN, N_MAX = 1.0e11, 1.0e20


def nearest_frame(frames, t_target):
    return min(range(len(frames)), key=lambda k: abs(frames[k][0] - t_target))


def main():
    dat   = sys.argv[1] if len(sys.argv) > 1 else "scenarios/data/loop_closed_gentle.dat"
    evout = sys.argv[2] if len(sys.argv) > 2 else "outputs/gentle_evap.txt"
    out   = sys.argv[3] if len(sys.argv) > 3 else "visualization/gentle_evap_profiles.png"
    t_pre = float(sys.argv[4]) if len(sys.argv) > 4 else 1000.0
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    topology = read_meta(dat).get("topology", "closed")
    ds, s_arc, phi_g, height = parse_dat(dat)
    frames = read_frames(evout)

    k_pre = nearest_frame(frames, t_pre)
    # During-frame: explicit 5th arg, else the post-ramp frame of PEAK TR/lower-
    # corona upflow (the active-evaporation moment, before the loop settles into
    # its new hotter equilibrium and the flows subside).
    if len(sys.argv) > 5:
        k_evp = nearest_frame(frames, float(sys.argv[5]))
    else:
        best, k_evp = -np.inf, len(frames) - 1
        for k in range(k_pre + 1, len(frames)):
            ni_k, _, v_k, Ti_k, _ = primitives(frames[k][1], phi_g)
            m = (Ti_k > 1.0e5) & (Ti_k < 2.0e6)
            vmax = v_k[m].max() if m.any() else -np.inf
            if vmax > best:
                best, k_evp = vmax, k
    t0, xn0 = frames[k_pre]
    t1, xn1 = frames[k_evp]
    ni0, nn0, v0, Ti0, _ = primitives(xn0, phi_g)
    ni1, nn1, v1, Ti1, _ = primitives(xn1, phi_g)

    if topology == "full":
        xMm = s_arc / 1e6
        xmax = float(ds.sum() / 1e6)
        xlabel = "arc length s [Mm]  (footpoint A | apex | footpoint B)"
    else:
        xMm = height / 1e6
        xmax = float(height.max() / 1e6)
        xlabel = "height above footpoint [Mm]  (zoom <2 Mm | corona compressed)"

    kd, kp, xticks = build_transform(xmax, topology)
    seps, labels, apex_x = region_layout(Ti0, xMm, topology, xmax)
    fwd = lambda x: np.interp(x, kd, kp)
    inv = lambda x: np.interp(x, kp, kd)

    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(f"Gentle conduction-driven evaporation — preflare (t={t0:.0f} s) "
                 f"vs evaporation (t={t1:.0f} s)   [{os.path.basename(dat)}]", fontsize=12)
    a_T, a_n, a_V = ax

    a_T.semilogy(xMm, np.clip(Ti0, T_MIN, T_MAX), color="0.55", lw=1.6, ls="--",
                 label=f"preflare (t={t0:.0f} s)")
    a_T.semilogy(xMm, np.clip(Ti1, T_MIN, T_MAX), color="crimson", lw=1.6,
                 label=f"evaporation (t={t1:.0f} s)")
    a_T.set_ylabel(r"$T_i$ [K]"); a_T.set_ylim(T_MIN, T_MAX)
    a_T.set_title("temperature"); a_T.legend(loc="lower right", fontsize=8)

    a_n.semilogy(xMm, np.clip(ni0, N_MIN, None), color="0.55", lw=1.6, ls="--",
                 label="preflare $n_e$")
    a_n.semilogy(xMm, np.clip(ni1, N_MIN, None), color="C3", lw=1.6, label="evaporation $n_e$")
    a_n.set_ylabel(r"$n_e$ [m$^{-3}$]"); a_n.set_ylim(N_MIN, N_MAX)
    a_n.set_title("electron density (corona fills)"); a_n.legend(loc="upper right", fontsize=8)

    a_V.plot(xMm, v0 / 1e3, color="0.55", lw=1.6, ls="--", label="preflare")
    a_V.plot(xMm, v1 / 1e3, color="C0", lw=1.6, label="evaporation")
    a_V.axhline(0, color="k", lw=0.6)
    a_V.set_ylabel(r"$V$ [km/s]  (+ up the loop = upflow)")
    a_V.set_title("velocity (evaporation = sustained +V)"); a_V.legend(loc="upper right", fontsize=8)

    edges = [kd[0]] + list(seps) + [kd[-1]]
    rtypes = (["chromo", "TR", "corona", "TR", "chromo"] if len(seps) == 4
              else ["chromo", "TR", "corona"])
    rcol = {"chromo": "#cfe0f2", "TR": "#d6efd6", "corona": "#fae1c6"}
    for a in ax:
        a.set_xscale("function", functions=(fwd, inv))
        a.set_xlim(kd[0], kd[-1])
        a.set_xticks([tk for tk in xticks if kd[0] <= tk <= kd[-1]])
        a.set_xlabel(xlabel)
        for x0, x1, rt in zip(edges[:-1], edges[1:], rtypes):
            a.axvspan(x0, x1, color=rcol[rt], alpha=0.5, lw=0, zorder=0)
        for sp in seps:
            a.axvline(sp, color="dimgray", ls="--", lw=0.8, alpha=0.6)
        if apex_x is not None:
            a.axvline(apex_x, color="seagreen", ls=":", lw=1.2, alpha=0.65)
        a.grid(True, alpha=0.18)
    trans = a_T.get_xaxis_transform()
    for xpos, lab in labels:
        a_T.text(xpos, 1.02, lab, transform=trans, ha="center", va="bottom",
                 fontsize=8, color="dimgray")

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")
    print(f"  preflare  t={t0:.0f}s: apex T={Ti0.max()/1e6:.2f} MK, "
          f"max +V={v0.max()/1e3:.1f} km/s")
    print(f"  evapor.   t={t1:.0f}s: apex T={Ti1.max()/1e6:.2f} MK, "
          f"max +V={v1.max()/1e3:.1f} km/s")


if __name__ == "__main__":
    main()
