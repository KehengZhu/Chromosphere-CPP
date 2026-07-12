#!/usr/bin/env python3
"""Evolution movie for a GENTLE conduction-driven evaporation run
(docs/gentle_evaporation_plan.md). The gentle counterpart of animate_loop_flare.py:
it reuses that module's topology-aware split x-axis (zoom the sub-Mm chromosphere+TR,
compress the extended corona; region bands marked) and parallel frame rendering, but
relabels the timeline for the heating-ramp driver instead of a beam:

  * RELAXING        — t < t_on: the loop settles to the steady preflare state under
                      the steady ambient coronal heating H(s) = E_H0·exp(−d/s_H).
  * HEATING RAMP ×N — t_on ≤ t ≤ t_on+ramp: H(s) is ramped up ×N.
  * EVAPORATION     — t > t_on+ramp: the enhanced conductive flux drives the (subsonic)
                      upflow that fills the corona.

The footpoint heating profile H(s,t)/E_H0 is overlaid (orange) on the T panel so you
can watch it ramp; the velocity panel is scaled for the gentle (tens of km/s) regime.

Usage:
  python util/animate_gentle.py <loop.dat> <run.txt> <out.mp4> [t_on] [ramp] [enhance]

Examples:
  # relaxation (steady heating, no ramp)
  python util/animate_gentle.py scenarios/data/loop_closed_gentle.dat \
      outputs/gentle_relax.txt visualization/gentle_relax_evolution.mp4 1e9 1 1
  # relaxation + ×3 ramp evaporation
  python util/animate_gentle.py scenarios/data/loop_closed_gentle.dat \
      outputs/gentle_evap.txt visualization/gentle_evap_evolution.mp4 1000 200 3
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _anim_parallel import save_frames_parallel
from animate_loop_flare import (read_meta, parse_dat, read_frames, primitives,
                                build_transform, region_layout,
                                T_MIN, T_MAX, N_MIN, N_MAX, MAX_FRAMES)

V_LIM = float(os.environ.get("GENTLE_VLIM_KMS", "30.0"))   # ± velocity-panel range [km/s]


def enhance_mult(t, t_on, ramp, enhance):
    """Heating amplitude multiplier enhance(t): 1 → enhance over [t_on, t_on+ramp]."""
    if enhance == 1.0:
        return 1.0
    if t <= t_on:
        return 1.0
    if t >= t_on + ramp:
        return enhance
    return 1.0 + (enhance - 1.0) * 0.5 * (1.0 - np.cos(np.pi * (t - t_on) / ramp))


def phase_of(t, t_on, ramp, enhance):
    if enhance == 1.0:
        return "RELAXING → STEADY", "0.30"
    if t < t_on:
        return "RELAXING", "0.40"
    if t <= t_on + ramp:
        return f"HEATING RAMP ×{enhance:.0f}", "orangered"
    return "EVAPORATION", "navy"


def _render_gentle_frame(k, tmpdir, xMm, t, Ti, Te, V, ni, nn, Hshape, hmult,
                         enhance, t_on, ramp, kd, kp, xticks, seps, labels,
                         apex_x, xlabel, vlabel, title_prefix):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    fwd = lambda x: np.interp(x, kd, kp)
    inv = lambda x: np.interp(x, kp, kd)

    fig, ax = plt.subplots(2, 2, figsize=(12.5, 8), sharex=True)
    phase, pcol = phase_of(t, t_on, ramp, enhance)
    fig.suptitle(f"{title_prefix}   t = {t:6.1f} s    [{phase}]",
                 fontsize=13, color=pcol, y=0.98)
    a_T, a_V, a_n, a_f = ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1]

    have_Te = not np.allclose(Te, Ti, rtol=1e-4)
    a_T.semilogy(xMm, np.clip(Te, T_MIN, T_MAX), color="dodgerblue", lw=1.5, ls="--",
                 label=r"$T_e$" if have_Te else r"$T_e=T_i$")
    a_T.semilogy(xMm, np.clip(Ti, T_MIN, T_MAX), color="crimson", lw=1.4, label=r"$T_i$")
    a_T.set_ylabel(r"$T$ [K]"); a_T.set_ylim(T_MIN, T_MAX)
    a_T.legend(loc="upper right", fontsize=8, ncol=2, framealpha=0.7)

    a_V.plot(xMm, V / 1e3, color="C0", lw=1.4); a_V.axhline(0, color="k", lw=0.6)
    a_V.set_ylabel(vlabel); a_V.set_ylim(-V_LIM, V_LIM)

    a_n.semilogy(xMm, np.clip(ni, N_MIN, None), color="C3", lw=1.3, label=r"$n_e$")
    a_n.semilogy(xMm, np.clip(nn, N_MIN, None), color="C0", lw=1.3, label=r"$n_n$")
    a_n.set_ylabel(r"density [m$^{-3}$]"); a_n.set_ylim(N_MIN, N_MAX)
    a_n.set_xlabel(xlabel); a_n.legend(loc="upper right", fontsize=9)

    a_f.plot(xMm, ni / np.maximum(ni + nn, 1.0), color="purple", lw=1.4)
    a_f.set_ylabel("ionization fraction"); a_f.set_ylim(-0.03, 1.05)
    a_f.set_xlabel(xlabel)

    edges = [kd[0]] + list(seps) + [kd[-1]]
    rtypes = (["chromo", "TR", "corona", "TR", "chromo"] if len(seps) == 4
              else ["chromo", "TR", "corona"])
    rcol = {"chromo": "#cfe0f2", "TR": "#d6efd6", "corona": "#fae1c6"}
    for a in (a_T, a_V, a_n, a_f):
        a.set_xscale("function", functions=(fwd, inv))
        a.set_xlim(kd[0], kd[-1])
        a.set_xticks([tk for tk in xticks if kd[0] <= tk <= kd[-1]])
        for x0, x1, rt in zip(edges[:-1], edges[1:], rtypes):
            a.axvspan(x0, x1, color=rcol[rt], alpha=0.55, lw=0, zorder=0)
        for sp in seps:
            a.axvline(sp, color="dimgray", ls="--", lw=0.9, alpha=0.6)
        if apex_x is not None:
            a.axvline(apex_x, color="seagreen", ls=":", lw=1.3, alpha=0.65)
        a.grid(True, alpha=0.18)

    # Footpoint heating overlay H(s,t)/E_H0 on the T panel (orange) — its height
    # tracks the ramp multiplier, so you see the heating grow as evaporation starts.
    aQ = a_T.twinx()
    aQ.fill_between(xMm, 0.0, Hshape * hmult, color="orange", alpha=0.30, lw=0)
    top = max(1.0, enhance) * 4.0          # fill occupies bottom ~1/4 at full ramp
    aQ.set_ylim(0.0, top); aQ.set_yticks([])
    aQ.set_zorder(a_T.get_zorder() - 1); a_T.patch.set_visible(False)

    trans = a_T.get_xaxis_transform()
    for xpos, lab in labels:
        a_T.text(xpos, 1.02, lab, transform=trans, ha="center", va="bottom",
                 fontsize=8, color="dimgray", alpha=0.85)
    if apex_x is not None:
        a_T.text(apex_x, 0.98, "apex", transform=trans, ha="center", va="top",
                 fontsize=8, color="seagreen", alpha=0.9)
    a_T.text(0.015, 0.93, f"heating H(s) = E$_0$·e$^{{-d/s_H}}$  (×{hmult:.2f})",
             transform=a_T.transAxes, fontsize=7.5, color="darkorange")

    fig.tight_layout()
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=110)
    plt.close(fig)


def main():
    dat   = sys.argv[1] if len(sys.argv) > 1 else "scenarios/data/loop_closed_gentle.dat"
    run   = sys.argv[2] if len(sys.argv) > 2 else "outputs/gentle_evap.txt"
    out   = sys.argv[3] if len(sys.argv) > 3 else "visualization/gentle_evap_evolution.mp4"
    t_on    = float(sys.argv[4]) if len(sys.argv) > 4 else float(os.environ.get("GENTLE_T_ON", 1e9))
    ramp    = float(sys.argv[5]) if len(sys.argv) > 5 else float(os.environ.get("GENTLE_RAMP", 30.0))
    enhance = float(sys.argv[6]) if len(sys.argv) > 6 else float(os.environ.get("GENTLE_ENHANCE", 1.0))
    sH_frac = float(os.environ.get("GENTLE_SH_FRAC", 0.4))
    os.makedirs(os.path.dirname(out), exist_ok=True)

    meta = read_meta(dat)
    topology = meta.get("topology", "closed")
    L = float(meta.get("loop_half_length_m", 0.0))
    ds, s_arc, phi_g, height = parse_dat(dat)
    Stot = float(ds.sum())
    if L <= 0.0:
        L = 0.5 * Stot if topology == "full" else Stot

    # Footpoint heating shape H(s)/E_H0 = exp(−d/s_H), peak-normalized for the overlay.
    s_H = sH_frac * L
    if topology == "full":
        d = np.minimum(s_arc, Stot - s_arc)
    else:
        d = s_arc.copy()
    Hshape = np.exp(-np.clip(d, 0.0, None) / s_H)
    Hshape /= Hshape.max()

    if topology == "full":
        xMm = s_arc / 1e6
        xmax = Stot / 1e6
        xlabel = "arc length s [Mm]   (footpoint A | apex | footpoint B)"
        vlabel = r"$V$ [km/s]  (+ along $s$: A$\to$apex$\to$B)"
        title_prefix = f"Gentle evaporation — full loop ({xmax:.0f} Mm)"
    else:
        xMm = height / 1e6
        xmax = float(height.max() / 1e6)
        xlabel = "height above footpoint [Mm]  (zoom <2 Mm | corona compressed)"
        vlabel = r"$V$ [km/s] (+ up the loop = upflow)"
        title_prefix = f"Gentle evaporation — half loop ({xmax:.0f} Mm)"

    frames = read_frames(run)
    if len(frames) > MAX_FRAMES:
        idx = np.linspace(0, len(frames) - 1, MAX_FRAMES).round().astype(int)
        frames = [frames[i] for i in idx]
    print(f"{run}: topology={topology}, {len(xMm)} cells, {len(frames)} frames, "
          f"t<= {frames[-1][0]:.1f} s; ramp: t_on={t_on:g} ramp={ramp:g} enhance={enhance:g}")

    kd, kp, xticks = build_transform(xmax, topology)
    _ni0, _nn0, _v0, Ti0, _Te0 = primitives(frames[0][1], phi_g)
    seps, labels, apex_x = region_layout(Ti0, xMm, topology, xmax)

    args = []
    for k, (t, xn) in enumerate(frames):
        ni, nn, v, Ti, Te = primitives(xn, phi_g)
        hmult = enhance_mult(t, t_on, ramp, enhance)
        args.append((k, xMm, t, Ti, Te, v, ni, nn, Hshape, hmult, enhance, t_on, ramp,
                     kd, kp, xticks, seps, labels, apex_x, xlabel, vlabel, title_prefix))
    fps = int(os.environ.get("ANIM_FPS", "30"))
    save_frames_parallel(_render_gentle_frame, args, out, fps=fps)


if __name__ == "__main__":
    main()
