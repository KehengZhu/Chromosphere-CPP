#!/usr/bin/env python3
"""Overlay synthetic loop Doppler shifts on the OBSERVED IRIS values for a flare.

Built for the data-driven SOL2014-03-29 X1.0 comparison: we configured the loop
(10 Mm quarter-circle, apex 1 MK) and the beam (RHESSI-derived delta=5, E_c=20 keV,
flux 5e10 erg/s/cm^2; Rubio da Costa et al. 2016) to the observed event, ran
model_flare, and here compare the synthesized Fe XXI / Si IV / C II / Mg II Doppler
shifts to what IRIS actually measured (Young et al. 2015 ApJ 799,218; Li et al.
2015 ApJ 811,7).

Reuses synth_doppler.py for the forward model (optically-thin EM-weighted v_Dopp,
disk-center LOS v_z = V*dh/ds). Produces:
  visualization/<stem>_vs_obs.png  — (left) synthetic v_Dopp(t) with observed
     bands overlaid; (right) peak synthetic vs observed-range comparison per line.

Usage:
  python util/compare_doppler_obs.py <loop.dat> <state.txt> [out_stem]
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from synth_doppler import (parse_dat, read_frames, read_meta, synth_timeseries,
                           LINES, T_ON, T_OFF)

# ---------------------------------------------------------------------------
# OBSERVED IRIS Doppler shifts for SOL2014-03-29 X1.0 (km/s; blue<0, red>0).
# Each entry: (lo, hi) observed range, a representative point, and the source.
# Fe XXI: Young 2015 Table 1 per-ribbon -100..-220 (strong ribbons), decaying to
#   ~0; loop-top <10 km/s. Li 2015: explosive pixel 117, profiles exceed 200.
# Si IV : Li 2015 explosive pixel 15 (single-Gaussian) .. 55 (red wing, double).
# C II  : Li 2015 ~36-37 (bisector/moment).
# Mg II : Li 2015 / general flare condensation, tens of km/s redshift.
# ---------------------------------------------------------------------------
OBS = {
    "Fe XXI 1354": dict(lo=-220.0, hi=-100.0, point=-150.0, decays_to_zero=True,
                        src="Young+2015 T1; Li+2015"),
    "Si IV 1394":  dict(lo=+15.0,  hi=+55.0,  point=+35.0,  decays_to_zero=False,
                        src="Li+2015 (1G..2G red wing)"),
    "C II 1334":   dict(lo=+30.0,  hi=+40.0,  point=+36.0,  decays_to_zero=False,
                        src="Li+2015 bisector/moment"),
    "Mg II k 2796": dict(lo=+10.0, hi=+30.0,  point=+20.0,  decays_to_zero=False,
                        src="flare condensation"),
}
COMPARE_LINES = ["Fe XXI 1354", "Si IV 1394", "C II 1334", "Mg II k 2796"]


def _line_index(label):
    return [i for i, L in enumerate(LINES) if L[0] == label][0]


def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    dat_path, txt_path = sys.argv[1], sys.argv[2]
    stem = sys.argv[3] if len(sys.argv) > 3 else \
        os.path.splitext(os.path.basename(txt_path))[0]

    meta = read_meta(dat_path)
    ds, s_arc, phi_g, height = parse_dat(dat_path)
    dhds = np.gradient(height, s_arc)
    frames = read_frames(txt_path)
    times, vDopp, Iline, hform = synth_timeseries(frames, phi_g, ds, dhds)
    print(f"[compare] {stem}: ns={len(ds)}, frames={len(frames)}, "
          f"t=[{times[0]:.1f},{times[-1]:.1f}] s, beam [{T_ON},{T_OFF}] s")

    fig = plt.figure(figsize=(14, 6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.7, 1.0], wspace=0.28)
    ax = fig.add_subplot(gs[0, 0])
    axc = fig.add_subplot(gs[0, 1])

    colors = {"Fe XXI 1354": "#c0392b", "Si IV 1394": "#e67e22",
              "C II 1334": "#2980b9", "Mg II k 2796": "#16a085"}

    # ---- Left: synthetic v_Dopp(t) with observed bands -------------------
    for lab in COMPARE_LINES:
        li = _line_index(lab)
        c = colors[lab]
        ax.plot(times, vDopp[:, li], color=c, lw=2.0, label=f"{lab} (synth)")
        o = OBS[lab]
        ax.axhspan(o["lo"], o["hi"], color=c, alpha=0.13, lw=0)
    ax.axvspan(T_ON, T_OFF, color="gold", alpha=0.15, lw=0, label="beam on")
    ax.axhline(0, color="k", lw=0.8)
    # Zoom the left panel to the dynamic window (spike + decay), not the long cooling tail.
    ax.set_xlim(0, min(float(times[-1]), T_OFF + 60.0))
    ax.set_xlabel("time  [s]")
    ax.set_ylabel("Doppler velocity  [km/s]   (blue<0 up / red>0 down)")
    ax.set_title(f"Synthetic vs observed IRIS Doppler — {stem}\n"
                 "shaded = observed range (SOL2014-03-29; Young+2015, Li+2015)")
    ax.legend(fontsize=8, ncol=2, loc="lower right")
    ax.grid(alpha=0.25)

    # ---- Right: peak synthetic vs observed range -------------------------
    labels = COMPARE_LINES[::-1]
    ypos = np.arange(len(labels))
    for y, lab in zip(ypos, labels):
        li = _line_index(lab)
        c = colors[lab]
        col = vDopp[:, li]
        col = col[np.isfinite(col)]
        if lab.startswith("Fe"):
            synth_peak = np.min(col)          # most-negative (blueshift)
        else:
            synth_peak = np.max(col)          # most-positive (redshift)
        o = OBS[lab]
        # observed range as a horizontal bar
        axc.plot([o["lo"], o["hi"]], [y, y], color=c, lw=8, alpha=0.30,
                 solid_capstyle="butt")
        axc.plot(o["point"], y, "o", color=c, ms=8, label="obs" if y == ypos[-1] else None)
        # synthetic peak as a diamond
        axc.plot(synth_peak, y, "D", color=c, ms=11, markeredgecolor="k",
                 markeredgewidth=1.2, label="synth peak" if y == ypos[-1] else None)
        axc.annotate(f"{synth_peak:+.0f}", (synth_peak, y),
                     textcoords="offset points", xytext=(0, 10), ha="center",
                     fontsize=8, color="k")
    axc.axvline(0, color="k", lw=0.8)
    axc.set_yticks(ypos); axc.set_yticklabels(labels, fontsize=9)
    axc.set_xlabel("peak Doppler velocity  [km/s]")
    axc.set_title("peak synth (◆) vs observed range (bar) + point (●)")
    axc.legend(fontsize=8, loc="lower left")
    axc.margins(y=0.18)
    axc.grid(alpha=0.25, axis="x")

    out = os.path.join("visualization", f"{stem}_vs_obs.png")
    os.makedirs("visualization", exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ---- Console summary table ------------------------------------------
    print("\n  line            synth peak   observed range     point    source")
    for lab in COMPARE_LINES:
        li = _line_index(lab)
        col = vDopp[:, li]; col = col[np.isfinite(col)]
        sp = np.min(col) if lab.startswith("Fe") else np.max(col)
        o = OBS[lab]
        print(f"  {lab:14s} {sp:+8.0f}    [{o['lo']:+.0f}, {o['hi']:+.0f}]"
              f"      {o['point']:+.0f}   {o['src']}")
    print(f"\n[compare] wrote {out}")


if __name__ == "__main__":
    main()
