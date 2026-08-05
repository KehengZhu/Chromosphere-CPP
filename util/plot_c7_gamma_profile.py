#!/usr/bin/env python3
"""Plot the Model C7 atmosphere (Avrett & Loeser 2008, Table 26) vs height:
temperature, mass density, and the CRASH effective adiabatic index gamma(T,n)
evaluated ALONG the C7 profile.

The C7 (h, n_e, n_HI, T) tables are parsed directly from scenarios/model_c7.cpp
(MODEL_C7_PHOTO + MODEL_C7) so this stays in sync with the code, and the profile
is sampled exactly as c7_full_profile() does (T linear, densities log-linear in h).
The effective gamma = 1 + P/e is bilinearly interpolated from the CRASH statistical-
sum EOS table (util/eos/tabulate_gamma.f90 -> outputs/eos_gamma/gamma_hydrogen.dat)
at each cell's (T, n_H) with n_H = n_HI + n_e (hydrogen nuclei).

Outputs:
  visualization/eos_gamma/c7_gamma_profile.png
  visualization/eos_gamma/c7_gamma_profile_classical_saha.png

Usage:
  python util/plot_c7_gamma_profile.py
  python util/plot_c7_gamma_profile.py --classical-saha
"""
import argparse
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CPP = "scenarios/model_c7.cpp"
GAMMA_TABLE = "outputs/eos_gamma/gamma_hydrogen.dat"
OUTDIR = "visualization/eos_gamma"
OUT = os.path.join(OUTDIR, "c7_gamma_profile.png")
OUT_CLASSICAL = os.path.join(OUTDIR, "c7_gamma_profile_classical_saha.png")

M_H = 1.6726219e-27          # proton mass [kg]
M_E = 9.1093837015e-31
K_B = 1.380649e-23
H_PLANCK = 6.62607015e-34
CHI_H = 2.179872361e-18
GAMMA_IDEAL = 5.0 / 3.0
GAMMA_MODEL = 1.05           # constant index used by the iso_t22k column configuration
H_TOP_KM = 2650.0            # top of the plot (photosphere -> upper TR / corona base)

# EOS-table column layout (see tabulate_gamma.f90)
COL = dict(T=0, Na=1, Rho=2, Z=3, P=4, E=5, G=6, GS=7, Ge=8, GSe=9, Cv=10)


def parse_c7_tables(path):
    """Return MODEL_C7_PHOTO and MODEL_C7 as (N,4) arrays: h[km], n_e, n_HI, T."""
    src = open(path).read()

    def grab(name):
        m = re.search(name + r"\[[^\]]*\]\[4\]\s*=\s*\{(.*?)\};", src, re.S)
        rows = re.findall(r"\{([^}]*)\}", m.group(1))
        clean = [[float(x.strip().rstrip("f")) for x in r.split(",")] for r in rows]
        return np.array(clean)

    return grab("MODEL_C7_PHOTO"), grab("MODEL_C7")


def c7_full_profile(h_km, photo, c7):
    """Vectorized replica of c7_full_profile() in model_c7.cpp.
    Concatenates PHOTO (h<1003) with C7 (h>=1003); T linear, n_e/n_HI log-linear."""
    tab = np.vstack([photo, c7])
    h, ne, nHI, T = tab[:, 0], tab[:, 1], tab[:, 2], tab[:, 3]
    Ti = np.interp(h_km, h, T)
    nei = np.exp(np.interp(h_km, h, np.log(ne)))
    nHIi = np.exp(np.interp(h_km, h, np.log(nHI)))
    return Ti, nei, nHIi


def load_gamma_grid(path, column=COL["G"]):
    """Return (Tvals[nT], Navals[nNa], Gamma[nNa,nT]) from the CRASH EOS table."""
    d = np.loadtxt(path)
    Na_all = d[:, COL["Na"]]
    Navals = np.unique(Na_all)
    Tvals = np.unique(d[Na_all == Navals[0], COL["T"]])
    nNa, nT = len(Navals), len(Tvals)
    G = d[:, column].reshape(nNa, nT)
    return Tvals, Navals, G


def gamma_at(T, Na, Tvals, Navals, G):
    """Bilinear interpolation of gamma on the log(T)-log(Na) grid, clamped to range."""
    lT, lTv = np.log10(np.clip(T, Tvals[0], Tvals[-1])), np.log10(Tvals)
    lN, lNv = np.log10(np.clip(Na, Navals[0], Navals[-1])), np.log10(Navals)

    def frac(x, xv):
        j = np.clip(np.searchsorted(xv, x) - 1, 0, len(xv) - 2)
        return j, (x - xv[j]) / (xv[j + 1] - xv[j])

    jN, fN = frac(lN, lNv)
    jT, fT = frac(lT, lTv)
    g00 = G[jN, jT];     g01 = G[jN, jT + 1]
    g10 = G[jN + 1, jT]; g11 = G[jN + 1, jT + 1]
    return ((1 - fN) * ((1 - fT) * g00 + fT * g01)
            + fN * ((1 - fT) * g10 + fT * g11))


def textbook_saha_gamma_energy(n_h, temperature):
    """Classical pure-H gamma_E=1+p/e with excitation/Fermi/Coulomb off."""
    log_c = 1.5 * np.log(2.0 * np.pi * M_E * K_B / H_PLANCK**2)
    log_a = (log_c + 1.5 * np.log(temperature)
             - CHI_H / (K_B * temperature) - np.log(n_h))
    x = np.empty_like(log_a)
    high = log_a >= 0.0
    x[high] = 2.0 / (1.0 + np.sqrt(1.0 + 4.0 * np.exp(-log_a[high])))
    u = np.exp(0.5 * log_a[~high])
    x[~high] = 2.0 * u / (u + np.sqrt(u * u + 4.0))
    pressure = (1.0 + x) * n_h * K_B * temperature
    energy = 1.5 * pressure + x * n_h * CHI_H
    return 1.0 + pressure / energy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--classical-saha", action="store_true",
        help="plot the signed-off classical pure-H Saha gamma_E and Gamma1")
    parser.add_argument("--output", help="override the output image path")
    args = parser.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    photo, c7 = parse_c7_tables(CPP)
    gamma_column = COL["GS"] if args.classical_saha else COL["G"]
    Tvals, Navals, G = load_gamma_grid(GAMMA_TABLE, gamma_column)

    h = np.linspace(0.0, H_TOP_KM, 1200)
    T, ne, nHI = c7_full_profile(h, photo, c7)
    nH = nHI + ne                       # hydrogen nuclei (n_HI + n_p, n_p ~ n_e)
    rho = M_H * nH                       # pure-H mass density [kg/m^3]
    gam = gamma_at(T, nH, Tvals, Navals, G)
    gamma_energy = textbook_saha_gamma_energy(nH, T) if args.classical_saha else gam

    # coverage note: EOS table starts at Na = 1e15; below that gamma is clamped.
    below = nH < Navals[0]
    if below.any():
        print("note: %d/%d samples have n_H < %.0e (gamma clamped)"
              % (below.sum(), len(h), Navals[0]))

    c_T, c_rho, c_g, c_g1 = "#1f6fb4", "#111111", "#c0392b", "#6f42c1"
    fig, axT = plt.subplots(figsize=(9.2, 5.8))
    fig.subplots_adjust(right=0.80)

    # left axis: temperature
    lT, = axT.plot(h, T, color=c_T, lw=2.2, label=r"temperature $T$")
    axT.set_yscale("log")
    axT.set_xlabel("height above photosphere  $h$  [km]")
    axT.set_ylabel(r"temperature  $T$  [K]", color=c_T)
    axT.tick_params(axis="y", colors=c_T)
    axT.spines["left"].set_color(c_T)
    axT.set_xlim(0, H_TOP_KM)

    # first right axis: mass density
    axR = axT.twinx()
    lR, = axR.plot(h, rho, color=c_rho, lw=2.0, ls="-", label=r"mass density $\rho$")
    axR.set_yscale("log")
    axR.set_ylabel(r"mass density  $\rho$  [kg m$^{-3}$]", color=c_rho)
    axR.tick_params(axis="y", colors=c_rho)
    axR.spines["right"].set_color(c_rho)

    # second right axis (offset): effective gamma
    axG = axT.twinx()
    axG.spines["right"].set_position(("axes", 1.16))
    if args.classical_saha:
        lG, = axG.plot(
            h, gamma_energy, color=c_g, lw=2.4,
            label=r"textbook Saha $\gamma_E=1+P/e$")
        lG1, = axG.plot(
            h, gam, color=c_g1, lw=2.1, ls="--",
            label=r"classical CRASH $\Gamma_1$")
    else:
        lG, = axG.plot(h, gam, color=c_g, lw=2.4,
                      label=r"effective $\gamma=1+P/e$")
        lG1 = None
    axG.set_ylabel(r"effective adiabatic index  $\gamma$", color=c_g)
    axG.tick_params(axis="y", colors=c_g)
    axG.spines["right"].set_color(c_g)
    axG.set_ylim(1.0, 1.71)
    axG.axhline(GAMMA_IDEAL, color=c_g, ls=":", lw=1.1, alpha=0.7)
    axG.text(H_TOP_KM * 0.985, GAMMA_IDEAL - 0.012, r"$5/3$ (ideal)",
             color=c_g, fontsize=8.5, ha="right", va="top", alpha=0.8)
    axG.axhline(GAMMA_MODEL, color=c_g, ls="--", lw=1.0, alpha=0.55)
    axG.text(H_TOP_KM * 0.985, GAMMA_MODEL + 0.012, r"model $\gamma=1.05$",
             color=c_g, fontsize=8.5, ha="right", va="bottom", alpha=0.7)

    # region bands: chromosphere / transition region
    axT.axvspan(0, 2153, color="gold", alpha=0.06)
    axT.axvspan(2153, H_TOP_KM, color="skyblue", alpha=0.06)
    axT.text(1080, axT.get_ylim()[1] * 0.45, "chromosphere",
             color="darkgoldenrod", fontsize=9, ha="center")
    axT.text(2400, axT.get_ylim()[1] * 0.45, "TR", color="steelblue",
             fontsize=9, ha="center")

    handles = [lT, lR, lG] + ([lG1] if lG1 is not None else [])
    axT.legend(handles=handles, loc="center left", fontsize=9.5,
               framealpha=0.9)
    if args.classical_saha:
        axT.set_title(
            "Model C7 with classical pure-H Saha EOS vs height\n"
            "excitation, Fermi, and Coulomb corrections off")
    else:
        axT.set_title("Model C7 atmosphere and CRASH effective $\\gamma(T,n_{\\rm H})$ "
                      "vs height\n(Avrett & Loeser 2008, Table 26)")
    axT.grid(True, which="both", ls=":", alpha=0.3)

    output = args.output or (OUT_CLASSICAL if args.classical_saha else OUT)
    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    fig.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(fig)
    if args.classical_saha:
        i_e = int(np.argmin(gamma_energy))
        i_1 = int(np.argmin(gam))
        print("textbook gamma_E min=%.6f at h=%.1f km" %
              (gamma_energy[i_e], h[i_e]))
        print("classical Gamma1 min=%.6f at h=%.1f km" % (gam[i_1], h[i_1]))
    print("wrote", output)


if __name__ == "__main__":
    main()
