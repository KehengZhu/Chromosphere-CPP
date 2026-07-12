#!/usr/bin/env python3
"""Writeup summary figure for gentle conduction-driven chromospheric evaporation
(model_gentle). Two panels:
  (a) the driver and response vs time — imposed coronal conductive flux q(t)
      (left axis) and the coronal emission measure EM(t)=∫_{h>2.3Mm} n_e^2 ds
      normalised to its pre-ramp value (right axis);
  (b) electron density vs height, preflare (relaxed) vs evaporation, on a split
      x-axis (zoom the sub-Mm chromosphere+TR, compress the corona) — the corona
      fills and the TR burns down.
Reads outputs/gentle_v3_evap.txt; writes docs/writeup-overleaf/gentle_evaporation.png.

Usage: python util/plot_gentle_summary.py [run.txt] [out.png]
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M_I, K_B = 1.6726219e-27, 1.380649e-23
CNI, CNV, CEI = 0, 2, 4
INK, PRIMARY, PRIMARY_L, ACCENT, GRID = "#1b1b1b", "#1f4e79", "#86a8cc", "#b4341f", "#c9ced6"
Q0, ENH, T_ON, RAMP = 230.0, 3.0, 3000.0, 300.0
ZOOM_END, CHROMO_FRAC = 2.6, 0.5

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.edgecolor": INK, "axes.linewidth": 0.8,
    "axes.spines.top": False, "legend.frameon": False,
})


def read_frames(path):
    raw = [l for l in open(path).read().splitlines() if l.strip()]
    ns, neq = (int(x) for x in raw[0].split()); h = np.fromstring(raw[1], sep=" ")
    fr, i = [], 2
    while i < len(raw):
        t = float(raw[i].split()[3]); i += 1
        if i + ns > len(raw): break
        d = np.array([np.fromstring(raw[i + j], sep=" ") for j in range(ns)]); i += ns
        fr.append((t, d))
    return h, fr


def w_ramp(t):
    x = (t - T_ON) / RAMP
    return np.where(x <= 0, 0.0, np.where(x >= 1, 1.0, 0.5 * (1 - np.cos(np.pi * x))))


def main():
    run = sys.argv[1] if len(sys.argv) > 1 else "outputs/gentle_v3_evap.txt"
    out = sys.argv[2] if len(sys.argv) > 2 else "docs/writeup-overleaf/gentle_evaporation.png"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    h, fr = read_frames(run)
    hMm = h / 1e3
    times = np.array([f[0] for f in fr])
    ne = [f[1][:, CNI] / M_I for f in fr]
    cor = h > 2300.0; dh = np.gradient(h * 1e3)
    EM = np.array([np.sum(n[cor] ** 2 * dh[cor]) for n in ne])
    EM0 = EM[int(np.argmin(np.abs(times - T_ON)))]
    i_pre = int(np.argmin(np.abs(times - T_ON)))
    i_evap = int(np.argmin(np.abs(times - 5000.0)))

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.0, 4.3))

    # ---- (a) driver q(t) + coronal-EM response ----
    tq = np.linspace(0, times[-1], 2000)
    q = Q0 * (1 + (ENH - 1) * w_ramp(tq))
    axA.axvspan(T_ON, T_ON + RAMP, color=ACCENT, alpha=0.12, lw=0)
    axA.plot(tq, q, color=PRIMARY, lw=2.0)
    axA.set_xlabel("time  $t$  [s]")
    axA.set_ylabel(r"imposed flux  $q(t)$  [W m$^{-2}$]", color=PRIMARY)
    axA.tick_params(axis="y", colors=PRIMARY)
    axA.set_ylim(0, ENH * Q0 * 1.15); axA.set_xlim(0, times[-1])
    axA.set_title("(a)  driver and coronal-EM response", loc="left")
    axR = axA.twinx(); axR.spines["top"].set_visible(False)
    axR.plot(times, EM / EM0, color=ACCENT, lw=2.0)
    axR.set_ylabel(r"coronal EM $/$ EM$_0$   ($\int n_e^2\,ds,\ h>2.3$ Mm)", color=ACCENT)
    axR.tick_params(axis="y", colors=ACCENT)
    axR.axhline(1.0, color=GRID, ls="--", lw=0.8)
    axR.set_ylim(0.8, max(2.2, (EM / EM0).max() * 1.08))
    axA.text(T_ON * 0.5, Q0 * 1.1, "relaxation", color=INK, fontsize=8.5, ha="center")
    axA.text((T_ON + times[-1]) / 2, ENH * Q0 * 1.06, "evaporation", color=INK, fontsize=8.5, ha="center")

    # ---- (b) n_e(h) preflare vs evaporation, split x-axis ----
    kd = np.array([hMm[0], ZOOM_END, hMm[-1]]); kp = np.array([0.0, CHROMO_FRAC, 1.0])
    fwd = lambda x: np.interp(x, kd, kp); inv = lambda x: np.interp(x, kp, kd)
    # region boundaries from preflare T
    Ti = (2 / 3 * fr[i_pre][1][:, CEI] - 1 / 3 * M_I * ne[i_pre] *
          (fr[i_pre][1][:, CNV] / fr[i_pre][1][:, CNI]) ** 2) / (2 * ne[i_pre] * K_B)
    def cross(L):
        j = int(np.argmax(Ti >= L)); return float(hMm[j]) if Ti[j] >= L else float(hMm[-1])
    ct, cb = cross(2e4), cross(5e5)
    for x0, x1, c in ((kd[0], ct, "#f1efea"), (ct, cb, "#fbe4cf"), (cb, kd[-1], "#e9eff6")):
        axB.axvspan(x0, x1, color=c, alpha=0.85, lw=0, zorder=0)
    axB.semilogy(hMm, np.clip(ne[i_pre], 1e11, None), color=PRIMARY_L, lw=1.8, ls="--",
                 label=f"preflare ($t={times[i_pre]:.0f}$ s)")
    axB.semilogy(hMm, np.clip(ne[i_evap], 1e11, None), color=PRIMARY, lw=1.8,
                 label=f"evaporation ($t={times[i_evap]:.0f}$ s)")
    axB.set_xscale("function", functions=(fwd, inv)); axB.set_xlim(kd[0], kd[-1])
    axB.set_xticks([1.0, 1.5, 2.0, 2.5, 4.0, 6.0, 8.0, 10.0])
    for sp in (ct, cb):
        axB.axvline(sp, color="dimgray", ls=":", lw=0.8, alpha=0.6)
    axB.set_xlabel("height  $h$  [Mm]"); axB.set_ylabel(r"$n_e$  [m$^{-3}$]")
    axB.set_title("(b)  corona fills, TR burns down", loc="left")
    axB.legend(loc="upper right", fontsize=8.5)
    axB.spines["right"].set_visible(False)

    fig.tight_layout()
    for ext in (".png", ".pdf"):
        fig.savefig(os.path.splitext(out)[0] + ext, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {os.path.splitext(out)[0]}.png/.pdf  EM peak ×{(EM/EM0).max():.2f}")


if __name__ == "__main__":
    main()
