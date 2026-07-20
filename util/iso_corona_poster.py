#!/usr/bin/env python3
"""Publication (nature-figure) styled snapshot figure of the resolved-corona,
conduction-driven evaporation testbed (`model_isentropic` with ISO_CORONA=1).

Three panels carry the Antiochos & Sturrock (1978) gentle-evaporation story:

    (a) Temperature  - profiles vs height; corona heats and the TR migrates up
    (b) Velocity     - V~0 chromosphere -> +79 km/s subsonic upflow   [hero]
    (c) Episode      - peak upflow & corona apex T vs time: the upflow builds
                       to t~194 s then the (no-sink) corona drains

Panels (a),(b) overlay six time "frames" colour-graded by time (cividis, navy->
maize, perceptually uniform & CVD-safe). Split x-axis (chromosphere+TR zoomed |
corona compressed), matching the movie and the project convention.

Usage:
    python util/iso_corona_poster.py [input.txt] [output.png]
Defaults: outputs/model_column/iso_corona.txt -> docs/poster-shine/web/iso_corona_evolution.png
"""
import os
import sys
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import animate_isentropic as A  # reuse load(), primitives(), heat_flux()

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    # Poster-sized: text must read at >=5 mm when shown ~299 mm wide on the
    # 914 mm poster. These render larger than the poster's own caption text.
    "font.size": 18,
    "axes.titlesize": 19,
    "axes.labelsize": 18,
    "xtick.labelsize": 16.5,
    "ytick.labelsize": 16.5,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 1.1,
    "legend.frameon": False,
})

# University of Michigan accent palette (matches the poster theme)
NAVY = "#00274C"
ORANGE = "#c05a1f"

# Six snapshots (s): a clean build-up to the peak upflow (t~194 s).
SNAP_TIMES = [0.0, 40.0, 80.0, 120.0, 160.0, 194.0]
T_PEAK_V = 194.0   # peak-upflow time (vertical marker in panel c)


def main():
    in_path = sys.argv[1] if len(sys.argv) > 1 else "outputs/model_column/iso_corona.txt"
    out_path = (sys.argv[2] if len(sys.argv) > 2
                else "docs/poster-shine/web/iso_corona_evolution.png")
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    h, frames = A.load(in_path)
    phi = A.G * (h - h[0]) * 1000.0
    times = np.array([fr[0] for fr in frames])

    # full-run time series for panel (c): peak upflow & corona apex T
    vmax_t = np.array([np.nanmax(A.primitives(fr[1], phi)[0]) / 1e3 for fr in frames])
    tcor_t = np.array([np.nanmax(A.primitives(fr[1], phi)[1]) / 1e6 for fr in frames])

    # six snapshot profiles for panels (a),(b)
    idxs = sorted(set(int(np.argmin(np.abs(times - t))) for t in SNAP_TIMES))
    snaps = []
    for i in idxs:
        V, T, rho, p, ne, nn = A.primitives(frames[i][1], phi)
        snaps.append((times[i], V, T))
    print("snapshots t =", [round(s[0], 1) for s in snaps], " ns =", h.size)

    # --- split x-axis: zoom chromosphere+TR into CHROMO_FRAC of the width --------
    CHROMO_FRAC = 0.55
    T0 = snaps[0][2]
    hot = np.where(T0 > 3.0e4)[0]
    z_split = float(h[hot[0]]) if hot.size else float(h[0] + 0.2 * (h[-1] - h[0]))
    kd = np.array([h[0], z_split, h[-1]])
    kp = np.array([0.0, CHROMO_FRAC, 1.0])
    fwd = lambda x: np.interp(x, kd, kp)
    inv = lambda y: np.interp(y, kp, kd)
    ch_ticks = [t for t in (0, 2000) if h[0] <= t <= z_split]
    co_ticks = [4000, 6000, 8000, 10000]
    xticks = np.array(ch_ticks + co_ticks)
    xlabels = [f"{int(t/1000)}" for t in (ch_ticks + co_ticks)]

    # --- time -> colour (cividis: navy -> maize) --------------------------------
    tmax = max(s[0] for s in snaps)
    norm = Normalize(0, tmax)
    cmap = mpl.colormaps["cividis"]
    colors = [cmap(norm(s[0])) for s in snaps]

    fig = plt.figure(figsize=(13.8, 5.2))
    gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 0.045], wspace=0.32,
                           left=0.055, right=0.94, bottom=0.205, top=0.88)
    axT = fig.add_subplot(gs[0, 0])
    axV = fig.add_subplot(gs[0, 1])
    axE = fig.add_subplot(gs[0, 2])
    cax = fig.add_subplot(gs[0, 3])

    # ---- panel (a) Temperature -------------------------------------------------
    for (t, V, T), c in zip(snaps, colors):
        axT.semilogy(h, T, color=c, lw=2.6)
    axT.set_ylabel("$T$  [K]")
    axT.set_title("(a)  Corona heats", loc="left", fontweight="bold")
    axT.set_ylim(3.5e3, 1.7e6)
    axT.annotate("TR rises", xy=(4300, 1.1e5), xytext=(4500, 9e3),
                 fontsize=16, color=NAVY, fontweight="bold", ha="left",
                 arrowprops=dict(arrowstyle="->", color=NAVY, lw=1.6))
    axT.text(0.05, 0.07, "chromo\n$\\sim$5 kK", transform=axT.transAxes,
             fontsize=15, color="0.4")

    # ---- panel (b) Velocity  (hero) -------------------------------------------
    for (t, V, T), c in zip(snaps, colors):
        axV.plot(h, V / 1e3, color=c, lw=2.6)
    axV.axhline(0, color="0.4", lw=0.9)
    axV.set_ylabel("$V$  [km s$^{-1}$]")
    axV.set_title("(b)  Gentle upflow", loc="left", fontweight="bold")
    axV.set_ylim(-14, 94)
    axV.annotate("+79 km s$^{-1}$\nMach $\\approx$ 0.5", xy=(9717, 79),
                 xytext=(3900, 50), fontsize=16, color=NAVY, fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color=NAVY, lw=1.8))
    axV.text(0.05, 0.12, "$V\\approx0$", transform=axV.transAxes,
             fontsize=15, color="0.4")

    # split-axis cosmetics shared by (a),(b); dashed line + "TR" mark the split
    for a in (axT, axV):
        a.set_xscale("function", functions=(fwd, inv))
        a.set_xlim(h[0], h[-1])
        a.set_xticks(xticks)
        a.set_xticklabels(xlabels)
        a.set_xlabel("height  [Mm]   (split at TR)")
        a.axvline(z_split, color="0.55", ls="--", lw=1.1, alpha=0.85)
        a.grid(True, which="major", alpha=0.22)
        a.tick_params(length=5)

    # ---- panel (c) episode time series ----------------------------------------
    axE.plot(times, vmax_t, color=NAVY, lw=3.0, label="peak upflow")
    axE.set_xlabel("time  [s]")
    axE.set_ylabel("peak $V$  [km s$^{-1}$]", color=NAVY)
    axE.tick_params(axis="y", labelcolor=NAVY)
    axE.set_title("(c)  Rise & drain", loc="left", fontweight="bold")
    axE.set_xlim(0, times[-1])
    axE.set_ylim(0, max(vmax_t) * 1.12)
    axE.grid(True, alpha=0.22)
    # snapshot markers (cividis) tie panel (c) back to the (a),(b) colour key
    for (t, V, Tt), c in zip(snaps, colors):
        j = int(np.argmin(np.abs(times - t)))
        axE.plot(times[j], vmax_t[j], "o", color=c, ms=11, mec="white", mew=1.4,
                 zorder=5)
    axE.axvline(T_PEAK_V, color="0.6", ls=":", lw=1.4)

    axE2 = axE.twinx()
    axE2.spines["top"].set_visible(False)
    axE2.plot(times, tcor_t, color=ORANGE, lw=3.0, ls="--")
    axE2.set_ylabel("apex $T$  [MK]", color=ORANGE)
    axE2.tick_params(axis="y", labelcolor=ORANGE)
    axE2.set_ylim(0.6, 1.35)
    axE.annotate("peaks,\nthen drains\n(no sink)", xy=(T_PEAK_V, max(vmax_t)),
                 xytext=(212, max(vmax_t) * 0.42), fontsize=15, color=NAVY,
                 arrowprops=dict(arrowstyle="->", color="0.4", lw=1.4))
    axE2.annotate("1.25 MK", xy=(times[np.argmax(tcor_t)], max(tcor_t)),
                  xytext=(120, 1.29), fontsize=15.5, color=ORANGE, fontweight="bold",
                  arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.4))

    # ---- shared time colorbar (panels a,b) -------------------------------------
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("snapshot time  [s]", fontsize=16)
    cb.set_ticks([s[0] for s in snaps])
    cb.ax.tick_params(labelsize=14)

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    stem = os.path.splitext(out_path)[0]
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    plt.close(fig)
    print("wrote", out_path, "+ .pdf")


if __name__ == "__main__":
    main()
