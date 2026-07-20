#!/usr/bin/env python3
"""Gentle-evaporation signature diagnostic (docs/gentle_evaporation_plan.md Phase 3).

Reads a chromo_main output .txt from a GENTLE run in which the ambient coronal
heating is ramped up (GENTLE_ENHANCE>1) from a relaxed preflare state, and shows
the classic gentle conduction-driven evaporation fingerprint:

  * coronal emission measure EM(t) = Σ n_e² ds  rises as evaporated chromospheric
    plasma fills the corona,
  * a sustained, SUBSONIC upflow (positive V — s increases upward, so +V = upflow)
    develops in the transition region / upper chromosphere after the heating step,
  * the apex / coronal-mean temperature rises to the new, hotter equilibrium.

The heating-ramp window [t_on, t_on+ramp] is shaded so the response is read as the
departure from the preflare baseline. Antiochos & Sturrock (1978): gentle ⇒ the
upflow stays well below the sound speed (Mach number shown).

Usage:
  python util/gentle_evaporation_diag.py outputs/model_column/gentle_evap.txt \
      [visualization/model_column/gentle_evap_signature.png] [t_on] [ramp]
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M_I = 1.6726219e-27
M_N = M_I
K_B = 1.380649e-23
GAMMA = 5.0 / 3.0
CNI, CNN, CNV, CNU, CEI, CEN = 0, 1, 2, 3, 4, 5


def read_frames(path):
    with open(path) as f:
        raw = [ln for ln in f.read().splitlines() if ln.strip()]
    ns, neq = (int(x) for x in raw[0].split())
    arc_km = np.fromstring(raw[1], sep=" ")
    frames, i = [], 2
    while i < len(raw):
        toks = raw[i].split()
        t, stp = float(toks[3]), int(toks[6])
        i += 1
        if i + ns > len(raw):
            break
        data = np.array([np.fromstring(raw[i + j], sep=" ") for j in range(ns)])
        i += ns
        frames.append((t, stp, data))
    return arc_km, frames


def primitives(xn):
    ni = xn[:, CNI] / M_I
    nn = xn[:, CNN] / M_N
    v  = xn[:, CNV] / (M_I * ni)
    p_i = 2.0 / 3.0 * xn[:, CEI] - 1.0 / 3.0 * M_I * ni * v * v
    Ti = p_i / (2.0 * ni * K_B)
    cs = np.sqrt(GAMMA * np.maximum(p_i, 0.0) / xn[:, CNI])   # ion sound speed
    return ni, nn, v, Ti, cs


def main():
    in_path = sys.argv[1] if len(sys.argv) > 1 else "outputs/model_column/gentle_evap.txt"
    out = sys.argv[2] if len(sys.argv) > 2 else "visualization/model_column/gentle_evap_signature.png"
    t_on = float(sys.argv[3]) if len(sys.argv) > 3 else 1000.0
    ramp = float(sys.argv[4]) if len(sys.argv) > 4 else 200.0
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    arc_km, frames = read_frames(in_path)
    # Cell widths from the cumulative-height column. chromo_main writes
    # arc_km[i] = (Σ_{j≤i} ds_j)/km + 1003 (a constant C7 base offset), so the
    # first face sits at 1003 km and ds[0] = arc_km[0] − 1003.
    BASE_KM = 1.003e3
    ds_m = np.empty_like(arc_km)
    ds_m[0] = (arc_km[0] - BASE_KM) * 1.0e3
    ds_m[1:] = np.diff(arc_km) * 1.0e3

    t   = np.array([fr[0] for fr in frames])
    EM  = np.zeros(len(frames))    # coronal column emission measure  [m^-5]
    vup = np.zeros(len(frames))    # peak TR/lower-corona upflow (+V)  [m/s]
    vdn = np.zeros(len(frames))    # peak downflow (−V)                [m/s]
    Mup = np.zeros(len(frames))    # Mach number of the peak upflow
    Tap = np.zeros(len(frames))    # apex (max) T
    Tco = np.zeros(len(frames))    # coronal-mean T (>0.5 MK)

    for k, (_, _, xn) in enumerate(frames):
        ni, nn, v, Ti, cs = primitives(xn)
        cor = Ti > 5.0e5
        EM[k]  = np.sum((ni[cor] ** 2) * ds_m[cor]) if cor.any() else 0.0
        Tap[k] = Ti.max()
        Tco[k] = Ti[cor].mean() if cor.any() else np.nan
        # Evaporation appears in the TR / lower corona (1e5–2e6 K): cool dense
        # chromospheric material heated and pushed up the loop.
        evap = (Ti > 1.0e5) & (Ti < 2.0e6)
        if evap.any():
            vv = v[evap]
            vup[k] = vv.max()
            kmax = np.argmax(vv)
            Mup[k] = vv[kmax] / max(cs[evap][kmax], 1.0)
        vdn[k] = -v.min()

    # Preflare baseline = mean over the ~100 s window just before the ramp.
    pre = (t > t_on - 100) & (t <= t_on)
    EM0 = EM[pre].mean() if pre.any() else EM[0]

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(f"Gentle conduction-driven evaporation — {os.path.basename(in_path)}",
                 fontsize=13)

    def shade(a):
        a.axvspan(t_on, t_on + ramp, color="orange", alpha=0.18,
                  label="heating ramp ×3")
        a.axvline(t_on, color="orange", ls="--", lw=0.8)

    ax[0, 0].plot(t, EM / max(EM0, 1e-30), "C2")
    shade(ax[0, 0]); ax[0, 0].axhline(1.0, color="k", lw=0.6)
    ax[0, 0].set_ylabel("coronal EM / preflare EM"); ax[0, 0].set_xlabel("t [s]")
    ax[0, 0].legend(fontsize=8); ax[0, 0].grid(alpha=0.3)
    ax[0, 0].set_title("emission-measure rise (corona fills with evaporate)")

    ax[0, 1].plot(t, vup / 1e3, "C0", label="peak upflow (+V)")
    ax[0, 1].plot(t, vdn / 1e3, "C3", lw=0.8, label="peak downflow (|−V|)")
    shade(ax[0, 1]); ax[0, 1].set_ylabel("V [km/s]"); ax[0, 1].set_xlabel("t [s]")
    ax[0, 1].legend(fontsize=8); ax[0, 1].grid(alpha=0.3)
    # Focus past the t<~300 s IC-settling transient so the gentle post-ramp upflow
    # (the departure from the relaxed baseline) is legible; clip the residual
    # single-cell TR spikes with a fixed y-range around the bulk evaporation speed.
    ax[0, 1].set_xlim(left=min(300.0, t_on - 200))
    vcap = max(40.0, 1.3 * vup[t > min(600.0, t_on - 100)].max() / 1e3)
    ax[0, 1].set_ylim(-5.0, vcap)
    ax[0, 1].set_title("TR/upper-chromosphere flow (evaporation = sustained +V)")

    ax[1, 0].plot(t, Tap / 1e6, "C2", label="apex T (max)")
    ax[1, 0].plot(t, Tco / 1e6, "C4", label="coronal mean T")
    shade(ax[1, 0]); ax[1, 0].set_ylabel("T [MK]"); ax[1, 0].set_xlabel("t [s]")
    ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=0.3)
    ax[1, 0].set_title("coronal temperature step")

    ax[1, 1].plot(t, Mup, "C1")
    shade(ax[1, 1]); ax[1, 1].axhline(1.0, color="k", ls="--", lw=0.8, label="sonic (M=1)")
    ax[1, 1].set_ylabel("upflow Mach number"); ax[1, 1].set_xlabel("t [s]")
    ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=0.3)
    ax[1, 1].set_xlim(left=min(300.0, t_on - 200))
    ax[1, 1].set_title("gentle regime ⇒ M ≪ 1 (Antiochos & Sturrock 1978)")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")
    post = t > t_on + ramp
    print(f"  preflare coronal EM     = {EM0:.3e} m^-5")
    if post.any():
        print(f"  post-ramp coronal EM    = {EM[post].mean():.3e} m^-5 "
              f"(×{EM[post].mean()/max(EM0,1e-30):.2f})")
        print(f"  peak upflow (post-ramp) = {vup[post].max()/1e3:.2f} km/s "
              f"(Mach {Mup[post].max():.2f})")
        print(f"  apex T: {Tap[pre].mean()/1e6:.2f} → {Tap[post].mean()/1e6:.2f} MK")


if __name__ == "__main__":
    main()
