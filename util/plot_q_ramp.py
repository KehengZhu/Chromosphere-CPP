#!/usr/bin/env python3
"""Plot the imposed coronal conductive-flux ramp q(t) used by model_gentle
(docs/q_of_T_boundary_condition.md). q(t) = q0 * [1 + (enhance-1) * w(t)] with the
raised-cosine (Hann) envelope w(t) implemented in scenarios/model_c7.cpp
(model_c7_update_bc). Defaults reproduce the gentle_v3_evolution run.

Usage:
  python util/plot_q_ramp.py [out.png] [q0] [enhance] [t_on] [ramp] [t_max]
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# restrained palette (matches util/animate_gentle_column.py)
INK, PRIMARY, ACCENT, GRID = "#1b1b1b", "#1f4e79", "#b4341f", "#c9ced6"
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 11, "axes.titlesize": 13, "axes.labelsize": 12,
    "axes.edgecolor": INK, "axes.linewidth": 0.8,
    "axes.spines.right": False, "axes.spines.top": False,
    "xtick.color": INK, "ytick.color": INK,
    "legend.frameon": False,
})


def w_ramp(t, t_on, ramp):
    """Raised-cosine 0->1 envelope (model_c7.cpp:489-491)."""
    x = (t - t_on) / ramp
    return np.clip(np.where(x <= 0.0, 0.0,
                   np.where(x >= 1.0, 1.0, 0.5 * (1.0 - np.cos(np.pi * x)))), 0.0, 1.0)


def main():
    out     = sys.argv[1] if len(sys.argv) > 1 else "visualization/model_column/q_ramp.png"
    q0      = float(sys.argv[2]) if len(sys.argv) > 2 else 230.0
    enhance = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
    t_on    = float(sys.argv[4]) if len(sys.argv) > 4 else 3000.0
    ramp    = float(sys.argv[5]) if len(sys.argv) > 5 else 300.0
    t_max   = float(sys.argv[6]) if len(sys.argv) > 6 else 5584.0
    os.makedirs(os.path.dirname(out), exist_ok=True)

    t = np.linspace(0.0, t_max, 4000)
    q = q0 * (1.0 + (enhance - 1.0) * w_ramp(t, t_on, ramp))
    q_hi = enhance * q0

    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    # ramp window highlight
    ax.axvspan(t_on, t_on + ramp, color=ACCENT, alpha=0.12, lw=0)
    # steady-level guides
    for lev, lab in ((q0, rf"$q_0 \approx {q0:.0f}$"), (q_hi, rf"${enhance:.0f}q_0 \approx {q_hi:.0f}$")):
        ax.axhline(lev, color=GRID, ls="--", lw=0.9, zorder=0)
        ax.text(t_max * 0.995, lev, "  " + lab, color=INK, fontsize=10, va="center", ha="left")
    # phase boundaries
    for tb in (t_on, t_on + ramp):
        ax.axvline(tb, color="dimgray", ls=":", lw=0.9, alpha=0.7)
    # the ramp itself
    ax.plot(t, q, color=PRIMARY, lw=2.2)

    # phase labels
    ax.text(t_on * 0.5, q0, "relaxation\n(steady $q_0$)", color=INK, fontsize=9.5,
            ha="center", va="bottom")
    ax.text(t_on + ramp / 2, 0.5 * (q0 + q_hi), "cosine\nramp", color=ACCENT,
            fontsize=9.5, ha="center", va="center", fontweight="bold")
    ax.text((t_on + ramp + t_max) / 2, q_hi, "evaporation", color=INK, fontsize=9.5,
            ha="center", va="bottom")

    ax.set_xlim(0, t_max); ax.set_ylim(0, q_hi * 1.18)
    ax.set_xlabel("time  $t$  [s]")
    ax.set_ylabel(r"imposed conductive flux  $q(t)$  [W m$^{-2}$]")
    ax.set_title(r"Gentle-evaporation driver: $q(t)=q_0\,[\,1+(\mathrm{enhance}-1)\,w(t)\,]$")
    ax.annotate(rf"$t_\mathrm{{on}}={t_on:.0f}$ s", xy=(t_on, q0), xytext=(t_on - 850, q0 * 1.55),
                color="dimgray", fontsize=9,
                arrowprops=dict(arrowstyle="->", color="dimgray", lw=0.8))
    ax.annotate(rf"$\Delta={ramp:.0f}$ s", xy=(t_on + ramp, q_hi), xytext=(t_on + ramp + 500, q_hi * 0.86),
                color="dimgray", fontsize=9,
                arrowprops=dict(arrowstyle="->", color="dimgray", lw=0.8))
    ax.grid(True, axis="y", alpha=0.15)
    fig.tight_layout()
    for ext in (".png", ".pdf"):
        fig.savefig(os.path.splitext(out)[0] + ext, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {os.path.splitext(out)[0]}.png / .pdf  "
          f"(q0={q0:.0f}, x{enhance:.0f} -> {q_hi:.0f} W/m^2, t_on={t_on:.0f}, ramp={ramp:.0f})")


if __name__ == "__main__":
    main()
