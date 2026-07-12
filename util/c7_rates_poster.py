#!/usr/bin/env python3
"""Poster (nature-figure) version of the Model-C7 ionization / recombination
rate figure — enlarged text, short legends, clean styling, sized to read at
poster scale (panel text >= ~5 mm when shown ~300 mm wide).

It reuses the physics and the computed rate arrays from the analysis script
`visualization/plot_c7_ioniz_recomb_rates.py` (importing it also regenerates
that script's own non-poster figure as a harmless side effect) and only changes
the presentation: three stacked panels (ionization, recombination,
photoionization-rate / Route-B closure) on the shared, TR-magnified height axis.

Usage:  python util/c7_rates_poster.py [output.png]
Default output: docs/poster-shine/web/c7_rates.png
"""
import os
import sys
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator

# Pull in the data + physics (module runs its own analysis on import).
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "visualization"))
import plot_c7_ioniz_recomb_rates as C7R
plt.close("all")   # discard the figure the imported module just built

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 16,
    "axes.titlesize": 18,
    "axes.labelsize": 17,
    "xtick.labelsize": 14.5,
    "ytick.labelsize": 14.5,
    "axes.linewidth": 1.1,
    "legend.frameon": True,
    "legend.framealpha": 0.85,
    "legend.fontsize": 13.5,
})

h, T = C7R.h_km, C7R.T
k_ext, h_extrap = C7R.k_ext, C7R.h_extrap
XLIM = C7R.XLIM
TR_STRETCH = C7R.TR_STRETCH
XTICKS = [0, 500, 1000, 1500, 2000, 2150, 2200]

# (rate, short label, colour, linestyle, marker, cool-plasma-fit?)
ion_series = [
    (C7R.R_phot, r"$n_n P_\mathrm{phot}$  (Route B)",  "tab:blue",  "-", "o", False),
    (C7R.R_mlvl, r"$n_e n_n S_\mathrm{CR}$  (CR)",      "tab:green", "-", "s", True),
    (C7R.R_coll, r"$n_e n_n S_i$  (Voronov)",           "tab:red",   "-", "^", False),
]
rec_series = [
    (C7R.R_rad,  r"$n_e^2\alpha_r$  (Hummer)",          "tab:blue",   "-",  "o", False),
    (C7R.R_3b,   r"$n_e^2\alpha_c$  (Hinnov)",          "tab:green",  "-",  "s", True),
    (C7R.R_neut, r"Drawin (neutral-stab.)",             "tab:orange", "--", "v", True),
    (C7R.R_mol,  r"H$_2^+$  (Cretu)",                   "tab:purple", ":",  "d", True),
]


def draw(ax, series):
    for y, lab, c, ls, mk, cool in series:
        if cool:
            ax.plot(h[:k_ext], y[:k_ext], ls, color=c, marker=mk, ms=5, lw=2.4,
                    markevery=3, label=lab)
            ax.plot(h[k_ext - 1:], y[k_ext - 1:], ls, color=c, lw=1.8, alpha=0.30)
        else:
            ax.plot(h, y, ls, color=c, marker=mk, ms=5, lw=2.4, markevery=3, label=lab)


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12.2, 12.6), sharex=True)
draw(ax1, ion_series)
draw(ax2, rec_series)

# panel 3: photoionization rate — Route B closure vs inversion vs Chae FAL-C
ax3.plot(h, C7R.P_phot_routeB, "-", color="0.15", lw=3.0,
         label="Route B (adopted)")
ax3.plot(h, C7R.P_phot_inv_plot, "--", color="tab:red", marker="^", ms=5, lw=2.0,
         markevery=3, label="C7 inversion")
ax3.plot(h, C7R.P_phot_falc, "-", color="tab:blue", marker="o", ms=5, lw=1.9,
         alpha=0.75, markevery=3, label="Chae 2021 FAL-C")

# shared x-axis (chromosphere 1:1, TR magnified) + log y + grid + markers
for ax in (ax1, ax2, ax3):
    ax.set_xscale("function", functions=(C7R._xfwd, C7R._xinv))
    ax.set_xlim(*XLIM)
    ax.set_xticks(XTICKS)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.25)
    ax.axvline(C7R.H_TMIN, color="0.4", ls=(0, (4, 3)), lw=1.2)
    ax.axvline(C7R.H_BREAK, color="0.6", ls=":", lw=1.1)
    ax.axvspan(h_extrap, XLIM[1], color="0.88", alpha=0.55, zorder=0)
    ax.tick_params(length=5)
for ax in (ax1, ax2):
    ax.set_ylim(1.0e6, 1.0e23)
ax3.set_ylim(1.0e-8, 1.0e-1)
ax3.axvline(C7R.H_EQ_NEQ, color="tab:red", ls=(0, (5, 2)), lw=1.4, alpha=0.7)

ax1.set_ylabel(r"rate  [m$^{-3}$ s$^{-1}$]")
ax2.set_ylabel(r"rate  [m$^{-3}$ s$^{-1}$]")
ax3.set_ylabel(r"$P_\mathrm{phot}$  [s$^{-1}$]")
ax1.set_title("(a)  Ionization channels", loc="left", fontweight="bold")
ax2.set_title("(b)  Recombination channels", loc="left", fontweight="bold")
ax3.set_title("(c)  Photoionization: Route B vs FAL-C", loc="left", fontweight="bold")
ax3.set_xlabel(r"height  [km]   ($h{=}0$ at $\tau_{500}{=}1$; TR $\times$%.0f)" % TR_STRETCH)
ax1.legend(loc="lower left")
ax2.legend(loc="lower left", ncol=2)
ax3.legend(loc="upper left")

# region labels across the top of panel (a) + TR band note on (a),(b)
ax1.text(280.0, 3.0e22, "photosphere", va="top", ha="center", fontsize=14, color="0.3")
ax1.text(1300.0, 3.0e22, "chromosphere", va="top", ha="center", fontsize=14, color="0.3")
band_mid = 0.5 * (h_extrap + XLIM[1])
for ax in (ax1, ax2):
    ax.text(band_mid, 3.0e22, "TR", ha="center", va="top", fontsize=14, color="0.3")
ax1.text(C7R.H_TMIN + 26, 6e6, r"$T_\mathrm{min}$", va="bottom", ha="left",
         fontsize=13, color="0.35", rotation=90)

# secondary top axis: temperature at selected heights
axT = ax1.twiny()
axT.set_xscale("function", functions=(C7R._xfwd, C7R._xinv))
axT.set_xlim(*XLIM)
tick_h = [0.0, C7R.H_TMIN, 1003.0, 1520.0, 1989.0, 2160.0, 2208.0]
idx = [int(np.argmin(np.abs(h - hh))) for hh in tick_h]
axT.set_xticks(h[idx])
axT.xaxis.set_minor_locator(NullLocator())
axT.set_xticklabels([(f"{T[i]/1e3:.0f} kK" if T[i] >= 2.0e4 else f"{T[i]:.0f} K")
                     for i in idx], fontsize=12.5)
axT.set_xlabel(r"$T_e$ at selected heights", fontsize=14)

fig.tight_layout()
out = sys.argv[1] if len(sys.argv) > 1 else "docs/poster-shine/web/c7_rates.png"
os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
fig.savefig(out, dpi=300, bbox_inches="tight")
fig.savefig(os.path.splitext(out)[0] + ".pdf", bbox_inches="tight")
plt.close(fig)
print("wrote", out, "+ .pdf")
