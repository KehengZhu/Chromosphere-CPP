#!/usr/bin/env python3
"""Publication (nature-figure) styled INITIAL-vs-FINAL snapshot figure of the
photosphere-to-TR-base, conduction-driven gentle-evaporation testbed
(`model_isentropic`, real Model-C7 IC extended to the 22 kK lower-TR base).

The companion to the `iso_corona` poster figure, but for the thin (2.15 Mm)
chromosphere->TR-base column that is run to a stable quasi-steady state. Two
"shots" -- the hydrostatic V=0 start and the final frame -- are OVERLAID on
each axis so the before/after change reads directly:

    (a) Temperature  - conduction warms the upper chromosphere 6.6 -> ~10 kK
    (b) Velocity     - V~0 -> a gentle +1 km/s subsonic upflow at the TR  [hero]
    (c) Density      - the TR/top fills (~6x) as material evaporates upward
    (d) Mass flux    - rho*V: a net upward (evaporative) mass flow, strongest
                       low down where the dense chromosphere is lifted

Single linear height axis in Mm (the slab is < 3 Mm, so no split axis); the
upper-chromosphere + TR "evaporation region" is shaded. Deeply subsonic upflow
=> gentle evaporation (Antiochos & Sturrock 1978), here spatially resolved
(not imposed) and grid-converged (ns=2000 vs 4000).

Usage:
    python util/iso_t22k_poster.py [input.txt] [output.pdf]
Defaults: outputs/iso_t22k_ns2000.txt -> docs/poster-shine/figs/iso_t22k_snapshots.pdf
"""
import os
import sys
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import animate_isentropic as A  # reuse load(), primitives()

# Adiabatic index of the run being plotted. The default conduction-only testbed
# now uses a reduced γ=1.05 (near-isothermal: c_v=k/(γ−1)→∞) as a polytropic
# stand-in for the radiative sink it omits. The .txt stores no γ, so the EOS
# decode (T, p from energy) MUST use the run's γ — override A's module globals.
GAMMA = float(os.environ.get("ISO_GAMMA", "1.05"))
A.GAMMA = GAMMA
A.GM1 = GAMMA - 1.0

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    # Poster-sized: text must read at >=5 mm when shown ~299 mm wide on the
    # 914 mm poster. Matches util/iso_corona_poster.py so the two Result figures
    # carry identical type and weight.
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

# University of Michigan accent palette (matches the poster theme & iso_corona fig).
# Cool navy = the V=0 "before"; warm orange = the heated/evaporating "after".
NAVY = "#00274C"     # initial  (t = 0)
ORANGE = "#c05a1f"   # final    (evaporating)

# Upper-chromosphere + TR window where conduction heats and the upflow develops.
EVAP_LO_MM = 1.5     # Mm -- below this the final state is unchanged from C7


def main():
    in_path = sys.argv[1] if len(sys.argv) > 1 else "outputs/iso_t22k_ns2000_gamma105.txt"
    out_path = (sys.argv[2] if len(sys.argv) > 2
                else "docs/poster-shine/figs/iso_t22k_snapshots.pdf")
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    h, frames = A.load(in_path)            # h in km
    phi = A.G * (h - h[0]) * 1000.0
    hmm = h / 1000.0                        # km -> Mm

    t0, arr0 = frames[0]
    tf, arrf = frames[-1]
    V0, T0, r0, p0, _, _ = A.primitives(arr0, phi)
    Vf, Tf, rf, pf, _, _ = A.primitives(arrf, phi)
    m0, mf = r0 * V0, rf * Vf      # field-aligned mass flux rho*V [kg/m^2/s]

    # peak upflow + sound-speed Mach (annotations)
    iup = int(np.nanargmax(Vf))
    cs = float(np.sqrt(GAMMA * pf[iup] / rf[iup]))
    mach = float(Vf[iup] / cs)
    # upper-chromosphere T (at ~1.8 Mm) initial vs final — the near-isothermal hold
    i18 = int(np.argmin(np.abs(hmm - 1.8)))
    print(f"gamma={GAMMA}  ns={h.size}  t0={t0:.1f}s  tf={tf:.1f}s  "
          f"peak V={Vf[iup]/1e3:+.2f} km/s @ h={hmm[iup]:.2f} Mm  Mach={mach:.3f}  "
          f"T(1.8Mm) {T0[i18]/1e3:.1f}->{Tf[i18]/1e3:.1f} kK  top rho x{rf[-1]/r0[-1]:.2f}")

    fig = plt.figure(figsize=(13.0, 7.0))
    gs = gridspec.GridSpec(2, 2, wspace=0.26, hspace=0.28,
                           left=0.075, right=0.985, bottom=0.085, top=0.85)
    axT = fig.add_subplot(gs[0, 0])
    axV = fig.add_subplot(gs[0, 1])
    axR = fig.add_subplot(gs[1, 0])
    axM = fig.add_subplot(gs[1, 1])

    LW_I, LW_F = 2.3, 3.2   # final drawn heavier (the "after" state)

    # ---- panel (a) Temperature -------------------------------------------------
    axT.plot(hmm, T0 / 1e3, color=NAVY, lw=LW_I)
    axT.plot(hmm, Tf / 1e3, color=ORANGE, lw=LW_F)
    axT.set_ylabel("$T$  [kK]")
    axT.set_title("(a)  Near-isothermal ($\\gamma{=}1.05$)", loc="left", fontweight="bold")
    axT.set_ylim(3.5, 22.5)
    axT.annotate("upper chromosphere\nheld at $\\sim6.6$ kK\n(heat absorbed at const. $T$)",
                 xy=(1.78, T0[i18] / 1e3), xytext=(0.42, 12.5), fontsize=14.5, color=NAVY,
                 fontweight="bold", ha="left",
                 arrowprops=dict(arrowstyle="->", color=NAVY, lw=1.6))
    axT.annotate("C7 $T$-min", xy=(0.56, 4.4), xytext=(0.95, 4.1), fontsize=14,
                 color="0.45", ha="left",
                 arrowprops=dict(arrowstyle="->", color="0.6", lw=1.2))

    # ---- panel (b) Velocity  (hero) -------------------------------------------
    axV.plot(hmm, V0 / 1e3, color=NAVY, lw=LW_I)
    axV.plot(hmm, Vf / 1e3, color=ORANGE, lw=LW_F)
    axV.axhline(0, color="0.55", lw=0.9, zorder=0)
    axV.set_ylabel("$V$  [km s$^{-1}$]")
    axV.set_title("(b)  Gentle upflow", loc="left", fontweight="bold")
    axV.set_ylim(-0.18, 1.28)
    axV.annotate(f"$+{Vf[iup]/1e3:.1f}$ km s$^{{-1}}$\nMach $\\approx{mach:.2f}$",
                 xy=(hmm[iup], Vf[iup] / 1e3), xytext=(0.85, 0.92), fontsize=15.5,
                 color=ORANGE, fontweight="bold", ha="left",
                 arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.8))
    axV.text(0.55, 0.12, "$V\\approx0$\n(initial)", fontsize=14.5, color=NAVY,
             ha="center")

    # ---- panel (c) Density -----------------------------------------------------
    axR.semilogy(hmm, r0, color=NAVY, lw=LW_I)
    axR.semilogy(hmm, rf, color=ORANGE, lw=LW_F)
    axR.set_ylabel(r"$\rho$  [kg m$^{-3}$]")
    axR.set_title("(c)  TR fills", loc="left", fontweight="bold")
    axR.annotate(f"top $\\rho\\times{rf[-1]/r0[-1]:.1f}$\n(mass evaporates up)",
                 xy=(hmm[-1], rf[-1]), xytext=(0.88, 8e-7), fontsize=14, color=ORANGE,
                 fontweight="bold", ha="left",
                 arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.6))

    # ---- panel (d) Mass flux  (mass flow) -------------------------------------
    # WHOLE PICTURE: full domain, full y-range. The dense base cells dominate ρV
    # (huge ρ × a residual O(Δs) drainage flow — a numerical artifact, not
    # evaporation), so the genuine TR evaporative flux (~10⁻⁶, ~100× smaller) is
    # invisible at this scale -> shown in the h>1.8 Mm INSET below.
    axM.plot(hmm, m0, color=NAVY, lw=LW_I)
    axM.plot(hmm, mf, color=ORANGE, lw=LW_F)
    axM.axhline(0, color="0.55", lw=0.9, zorder=0)
    axM.set_ylabel(r"$\rho V$  [kg m$^{-2}$ s$^{-1}$]")
    axM.set_title("(d)  Mass flow", loc="left", fontweight="bold")
    axM.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    finite = np.isfinite(mf) & np.isfinite(m0)
    allflux = np.concatenate([m0[finite], mf[finite]])
    fhi, flo = float(allflux.max()), float(allflux.min())
    fp = 0.10 * (fhi - flo)
    axM.set_ylim(flo - fp, fhi + fp)
    ibase = int(np.nanargmin(mf))             # the large dense-base spike dominates
    # Label sits in the empty lower-left and points with a SHORT arrow to the spike;
    # kept well clear of the lower-right inset.
    axM.annotate("dense-base flux\n(numerical $O(\\Delta s)$)", xy=(hmm[ibase], flo * 0.92),
                 xytext=(0.16, flo * 0.74), fontsize=12.5, color="0.35",
                 fontweight="bold", ha="left", va="center",
                 arrowprops=dict(arrowstyle="->", color="0.5", lw=1.4))

    # INSET: zoom the h>1.8 Mm TR, where the converged evaporative flux lives.
    HZ = 1.8                                  # Mm
    sel = hmm >= HZ
    axIns = inset_axes(axM, width="48%", height="42%", loc="lower right",
                       borderpad=1.4)
    axIns.plot(hmm[sel], m0[sel], color=NAVY, lw=2.0)
    axIns.plot(hmm[sel], mf[sel], color=ORANGE, lw=2.6)
    axIns.axhline(0, color="0.55", lw=0.8, zorder=0)
    tr = np.concatenate([m0[sel], mf[sel]])
    tlo, thi = float(np.nanmin(tr)), float(np.nanmax(tr))
    tp = 0.12 * (thi - tlo)
    axIns.set_xlim(HZ, hmm[-1]); axIns.set_ylim(tlo - tp, thi + tp)
    axIns.set_title("zoom: $h>1.8$ Mm (TR)", fontsize=12.5, color=ORANGE, pad=2)
    axIns.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axIns.tick_params(labelsize=11.5, length=3)
    axIns.yaxis.get_offset_text().set_size(10.5)
    for s in axIns.spines.values():
        s.set_edgecolor(ORANGE); s.set_linewidth(1.1); s.set_visible(True)
    mark_inset(axM, axIns, loc1=2, loc2=1, fc="none", ec="0.55", lw=0.9)
    # Label sits in the empty band BELOW the (high) orange curve, clear of it.
    axIns.text(HZ + 0.04, thi * 0.30, "evaporative\nupflow",
               fontsize=11.5, color=ORANGE, fontweight="bold",
               ha="left", va="center")

    # ---- shared cosmetics: shade evaporation region, height axis, grid ----------
    for a in (axT, axV, axR, axM):
        a.axvspan(EVAP_LO_MM, hmm[-1], color=NAVY, alpha=0.05, lw=0, zorder=0)
        a.set_xlim(hmm[0], hmm[-1])
        a.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0])
        a.grid(True, which="major", alpha=0.22)
        a.tick_params(length=5)
    for a in (axR, axM):              # bottom row carries the height label
        a.set_xlabel("height  [Mm]")
    for a in (axT, axV):             # top row shares the x-axis -> hide its labels
        a.tick_params(labelbottom=False)
    # one "TR / evaporation region" label, in the temperature panel
    axT.text(1.83, 21.0, "TR", fontsize=15, color="0.4", ha="center",
             fontweight="bold")

    # ---- shared initial/final legend (figure level, top centre) ----------------
    handles = [Line2D([0], [0], color=NAVY, lw=LW_I,
                      label=f"initial  ($t=0$)"),
               Line2D([0], [0], color=ORANGE, lw=LW_F,
                      label=f"final  ($t={tf:.0f}$ s)")]
    fig.legend(handles=handles, ncol=2, loc="upper center",
               bbox_to_anchor=(0.5, 0.975), fontsize=16.5, handlelength=1.6,
               columnspacing=2.2)

    fig.savefig(out_path, bbox_inches="tight")
    stem = os.path.splitext(out_path)[0]
    fig.savefig(stem + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out_path, "+ .png")


if __name__ == "__main__":
    main()
