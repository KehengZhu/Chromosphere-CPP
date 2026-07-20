#!/usr/bin/env python3
"""Plot the effective adiabatic index gamma tabulated by the CRASH statistical-sum
EOS (see util/eos/tabulate_gamma.f90).

Produces, under visualization/eos_gamma/:
  gamma_hydrogen.png     2-panel gamma(T) + ionization(T) for several densities
  gamma_hydrogen_2d.png  gamma(T, Na) heatmap with ionization contours
  gamma_H_vs_HHe.png     pure-H vs H/He-mix comparison at one density

Usage:
  python util/plot_eos_gamma.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors

# table columns written by tabulate_gamma.f90
# 0 T[K] 1 Na 2 Rho 3 Z 4 P 5 Edens 6 Gamma 7 GammaS 8 Gammae 9 GammaSe 10 Cv
COL = dict(T=0, Na=1, Rho=2, Z=3, P=4, E=5, G=6, GS=7, Ge=8, GSe=9, Cv=10)

H_TABLE = "outputs/eos_gamma/gamma_hydrogen.dat"
HE_TABLE = "outputs/eos_gamma/gamma_H90He10.dat"
OUTDIR = "visualization/eos_gamma"
GAMMA_MODEL = 1.05          # the constant currently used in the chromosphere model
GAMMA_IDEAL = 5.0 / 3.0

# densities (Na, m^-3) to draw as lines
NA_LINES = [1e16, 1e18, 1e20, 1e22]


def load_grid(path):
    """Return (Tvals, Navals, dict of 2D arrays[nNa, nT])."""
    d = np.loadtxt(path)
    Na_all = d[:, COL["Na"]]
    Navals = np.unique(Na_all)
    T0 = d[Na_all == Navals[0], COL["T"]]
    Tvals = np.unique(T0)
    nNa, nT = len(Navals), len(Tvals)
    grids = {}
    for key, c in COL.items():
        grids[key] = d[:, c].reshape(nNa, nT)
    return Tvals, Navals, grids


def nearest_na(Navals, target):
    return int(np.argmin(np.abs(np.log10(Navals) - np.log10(target))))


def plot_lines(Tvals, Navals, g, out):
    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(8.2, 8.0), sharex=True,
        gridspec_kw=dict(height_ratios=[2.0, 1.0], hspace=0.08))

    cmap = matplotlib.colormaps["viridis"]
    cnorm = colors.Normalize(vmin=15.5, vmax=22.5)

    for na in NA_LINES:
        i = nearest_na(Navals, na)
        c = cmap(cnorm(np.log10(Navals[i])))
        lab = r"$n_{\rm H}=10^{%.0f}\,{\rm m^{-3}}$" % np.log10(Navals[i])
        ax0.plot(Tvals, g["G"][i], color=c, lw=2.0, label=lab)
        ax0.plot(Tvals, g["GS"][i], color=c, lw=1.3, ls="--")
        ax1.plot(Tvals, g["Z"][i], color=c, lw=2.0)

    # reference levels
    ax0.axhline(GAMMA_IDEAL, color="0.4", ls=":", lw=1.2)
    ax0.text(2.3e3, GAMMA_IDEAL + 0.005, r"$5/3$ (ideal, mono-atomic)",
             color="0.4", fontsize=9, va="bottom")
    ax0.axhline(GAMMA_MODEL, color="crimson", ls=":", lw=1.4)
    ax0.text(2.3e3, GAMMA_MODEL + 0.006, r"model default $\gamma=1.05$",
             color="crimson", fontsize=9, va="bottom")

    # legend: densities (color) + line-style meaning
    leg1 = ax0.legend(loc="center right", fontsize=9, framealpha=0.9,
                      title="density")
    ax0.add_artist(leg1)
    style_handles = [
        plt.Line2D([], [], color="0.2", lw=2.0,
                   label=r"energy $\gamma=1+P/e$  (closes $e=P/(\gamma-1)$)"),
        plt.Line2D([], [], color="0.2", lw=1.3, ls="--",
                   label=r"adiabatic $\Gamma_1=(\partial\ln P/\partial\ln\rho)_S$"),
    ]
    ax0.legend(handles=style_handles, loc="upper right", fontsize=8.5,
               framealpha=0.9)

    ax0.set_ylabel(r"effective adiabatic index $\gamma$")
    ax0.set_ylim(1.0, 1.71)
    ax0.set_title("Effective $\\gamma$ of partially ionized hydrogen  "
                  "(CRASH Saha EOS: ionization + excitation + degeneracy)")

    ax1.set_ylabel(r"ionization $\langle Z\rangle$")
    ax1.set_ylim(-0.03, 1.05)
    ax1.set_xlabel("temperature  T  [K]")
    ax1.set_xscale("log")
    ax1.set_xlim(2e3, 5e6)

    # shade the chromospheric temperature band
    for ax in (ax0, ax1):
        ax.axvspan(4e3, 2.5e4, color="gold", alpha=0.10)
        ax.grid(True, which="both", ls=":", alpha=0.35)
    ax1.text(1.0e4, 0.06, "chromosphere", color="darkgoldenrod",
             fontsize=9, ha="center")

    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def plot_heatmap(Tvals, Navals, g, out):
    fig, ax = plt.subplots(figsize=(8.4, 6.2))
    TT, NN = np.meshgrid(Tvals, Navals)
    pcm = ax.pcolormesh(TT, NN, g["G"], shading="auto",
                        cmap="magma_r", vmin=1.0, vmax=5.0 / 3.0)
    cb = fig.colorbar(pcm, ax=ax, extend="neither")
    cb.set_label(r"energy $\gamma = 1 + P/e$")

    # ionization contours
    cs = ax.contour(TT, NN, g["Z"], levels=[0.1, 0.5, 0.9],
                    colors="cyan", linewidths=1.0)
    ax.clabel(cs, fmt=lambda v: r"$Z=%.1f$" % v, fontsize=8)
    # locus of minimum gamma vs density
    imin = np.argmin(g["G"], axis=1)
    ax.plot(Tvals[imin], Navals, color="white", lw=1.5, ls="--",
            label=r"$\gamma$ minimum")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("temperature  T  [K]")
    ax.set_ylabel(r"hydrogen number density  $n_{\rm H}$  [m$^{-3}$]")
    ax.set_title(r"Effective $\gamma(T, n_{\rm H})$ — partially ionized hydrogen")
    ax.legend(loc="upper right", fontsize=9)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def plot_compare(out):
    Th, Nah, gh = load_grid(H_TABLE)
    The, Nahe, ghe = load_grid(HE_TABLE)
    na = 1e19
    ih, ihe = nearest_na(Nah, na), nearest_na(Nahe, na)

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.plot(Th, gh["G"][ih], color="navy", lw=2.2,
            label="pure H  (energy $\\gamma$)")
    ax.plot(The, ghe["G"][ihe], color="darkorange", lw=2.2,
            label="H 90% / He 10%  (energy $\\gamma$)")
    ax.plot(Th, gh["GS"][ih], color="navy", lw=1.3, ls="--", alpha=0.8)
    ax.plot(The, ghe["GS"][ihe], color="darkorange", lw=1.3, ls="--", alpha=0.8)

    ax.axhline(GAMMA_IDEAL, color="0.4", ls=":", lw=1.2)
    ax.axhline(GAMMA_MODEL, color="crimson", ls=":", lw=1.4)
    ax.text(2.3e3, GAMMA_MODEL + 0.006, r"$\gamma=1.05$", color="crimson",
            fontsize=9, va="bottom")

    ax.set_xscale("log")
    ax.set_xlim(2e3, 5e6)
    ax.set_ylim(1.0, 1.71)
    ax.set_xlabel("temperature  T  [K]")
    ax.set_ylabel(r"effective adiabatic index $\gamma$")
    ax.set_title(r"H vs H/He mix at $n_{\rm H}=10^{19}\,{\rm m^{-3}}$ "
                 "(dashed = adiabatic $\\Gamma_1$)")
    ax.grid(True, which="both", ls=":", alpha=0.35)
    ax.legend(loc="center right", fontsize=9)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    Tvals, Navals, g = load_grid(H_TABLE)
    plot_lines(Tvals, Navals, g, os.path.join(OUTDIR, "gamma_hydrogen.png"))
    plot_heatmap(Tvals, Navals, g, os.path.join(OUTDIR, "gamma_hydrogen_2d.png"))
    plot_compare(os.path.join(OUTDIR, "gamma_H_vs_HHe.png"))


if __name__ == "__main__":
    main()
