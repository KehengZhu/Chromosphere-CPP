#!/usr/bin/env python3
"""Plot C7, classical-Saha, and CRASH hydrogen ionization degree vs height.

The atmospheric profile is parsed from scenarios/model_c7.cpp and sampled with
the same interpolation used by c7_full_profile().  With no arguments, this
script writes the two individual profiles and their direct comparison:

  visualization/eos_gamma/c7_ionization_profile.png
  visualization/eos_gamma/c7_ionization_profile_classical_saha.png
  visualization/eos_gamma/c7_ionization_profile_compare.png
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Rectangle
import numpy as np

from plot_c7_gamma_profile import (
    CPP,
    H_TOP_KM,
    K_B,
    M_E,
    M_H,
    H_PLANCK,
    CHI_H,
    COL,
    c7_full_profile,
    gamma_at,
    load_gamma_grid,
    parse_c7_tables,
)


OUTDIR = "visualization/eos_gamma"
CRASH_IONIZATION_TABLE = "outputs/eos_gamma/gamma_hydrogen_excitation.dat"
OUTPUTS = {
    "c7": os.path.join(OUTDIR, "c7_ionization_profile.png"),
    "classical-saha": os.path.join(
        OUTDIR, "c7_ionization_profile_classical_saha.png"),
    "compare": os.path.join(OUTDIR, "c7_ionization_profile_compare.png"),
}


def saha_ionization_fraction(n_h, temperature):
    """Stable textbook pure-H Saha ionization degree x=n_p/n_H."""
    log_c = 1.5 * np.log(2.0 * np.pi * M_E * K_B / H_PLANCK**2)
    log_a = (log_c + 1.5 * np.log(temperature)
             - CHI_H / (K_B * temperature) - np.log(n_h))
    x = np.empty_like(log_a)
    high = log_a >= 0.0
    x[high] = 2.0 / (1.0 + np.sqrt(1.0 + 4.0 * np.exp(-log_a[high])))
    u = np.exp(0.5 * log_a[~high])
    x[~high] = 2.0 * u / (u + np.sqrt(u * u + 4.0))
    return x


def make_plot(mode, h, temperature, rho, x_c7, x_saha, x_crash, output):
    c_t, c_rho = "#1f6fb4", "#111111"
    c_c7, c_saha, c_crash = "#c0392b", "#6f42c1", "#008b8b"

    fig, ax_t = plt.subplots(figsize=(9.2, 5.8))
    fig.subplots_adjust(right=0.80)

    line_t, = ax_t.plot(
        h, temperature, color=c_t, lw=2.2, label=r"temperature $T$")
    ax_t.set_yscale("log")
    ax_t.set_xlabel("height above photosphere  $h$  [km]")
    ax_t.set_ylabel(r"temperature  $T$  [K]", color=c_t)
    ax_t.tick_params(axis="y", colors=c_t)
    ax_t.spines["left"].set_color(c_t)
    ax_t.set_xlim(0.0, H_TOP_KM)

    ax_rho = ax_t.twinx()
    line_rho, = ax_rho.plot(
        h, rho, color=c_rho, lw=2.0, label=r"mass density $\rho$")
    ax_rho.set_yscale("log")
    ax_rho.set_ylabel(r"mass density  $\rho$  [kg m$^{-3}$]", color=c_rho)
    ax_rho.tick_params(axis="y", colors=c_rho)
    ax_rho.spines["right"].set_color(c_rho)

    ax_x = ax_t.twinx()
    ax_x.spines["right"].set_position(("axes", 1.16))
    ion_lines = []
    if mode in ("c7", "compare"):
        line_c7, = ax_x.plot(
            h, x_c7, color=c_c7, lw=2.4,
            label=r"C7 $x=n_e/(n_{\rm HI}+n_e)$")
        ion_lines.append(line_c7)
    if mode in ("classical-saha", "compare"):
        line_saha, = ax_x.plot(
            h, x_saha, color=c_saha, lw=2.2,
            ls="--" if mode == "compare" else "-",
            label=r"classical Saha $x=n_p/n_{\rm H}$")
        ion_lines.append(line_saha)
    if mode == "compare":
        line_crash, = ax_x.plot(
            h, x_crash, color=c_crash, lw=1.8, ls=":",
            marker="o", markevery=0.05, ms=3.0, mfc="white", mew=0.9,
            label=r"CRASH excitation-on $\bar{Z}$")
        ion_lines.append(line_crash)
    ion_axis_color = c_saha if mode == "classical-saha" else c_c7
    ax_x.set_yscale("log")
    ax_x.set_ylim(1.0e-6, 1.2)
    ax_x.set_ylabel(r"hydrogen ionization degree  $x$", color=ion_axis_color)
    ax_x.tick_params(axis="y", colors=ion_axis_color)
    ax_x.spines["right"].set_color(ion_axis_color)
    for level in (0.01, 0.5, 0.9):
        ax_x.axhline(level, color="#777777", ls=":", lw=0.8, alpha=0.35)

    if mode == "compare":
        zoom_lo, zoom_hi = 2119.0, 2147.0
        ax_zoom = ax_t.inset_axes([0.57, 0.45, 0.31, 0.24])
        ax_zoom.plot(h, x_c7, color=c_c7, lw=2.2, label="C7")
        ax_zoom.plot(
            h, x_saha, color=c_saha, lw=2.0, ls="--",
            label="classical Saha")
        ax_zoom.plot(
            h, x_crash, color=c_crash, lw=1.8, ls=":",
            label="CRASH excitation-on")
        ax_zoom.set_xlim(zoom_lo, zoom_hi)
        zoom_mask = (h >= zoom_lo) & (h <= zoom_hi)
        zoom_values = np.concatenate((
            x_c7[zoom_mask], x_saha[zoom_mask], x_crash[zoom_mask]))
        zoom_pad = 0.05 * (np.max(zoom_values) - np.min(zoom_values))
        ax_zoom.set_ylim(
            np.min(zoom_values) - zoom_pad,
            min(1.01, np.max(zoom_values) + zoom_pad))
        ax_zoom.set_title("Zoom: 2119–2147 km", fontsize=9)
        ax_zoom.set_xlabel(r"$h$ [km]", fontsize=8)
        ax_zoom.set_ylabel("ionization degree", fontsize=8)
        ax_zoom.tick_params(labelsize=7)
        ax_zoom.grid(True, ls=":", alpha=0.35)
        ax_zoom.legend(
            loc="upper left", fontsize=6.8, framealpha=0.9,
            handlelength=2.2)
        box_lo, box_hi = zoom_lo - 20.0, zoom_hi + 20.0
        box_y_lo, box_y_hi = 0.40, 1.16
        zoom_box = Rectangle(
            (box_lo, box_y_lo), box_hi - box_lo, box_y_hi - box_y_lo,
            fill=False, edgecolor="#e74c3c", linewidth=1.5, zorder=8)
        ax_x.add_patch(zoom_box)
        zoom_arrow = ConnectionPatch(
            xyA=((box_lo + box_hi) / 2.0, box_y_lo),
            coordsA=ax_x.transData,
            xyB=(0.52, 1.11), coordsB=ax_zoom.transAxes,
            arrowstyle="->", color="#e74c3c", linewidth=1.4,
            mutation_scale=11, shrinkA=1.5, shrinkB=2.0, zorder=9)
        fig.add_artist(zoom_arrow)

    ax_t.axvspan(0.0, 2153.0, color="gold", alpha=0.06)
    ax_t.axvspan(2153.0, H_TOP_KM, color="skyblue", alpha=0.06)
    ax_t.text(
        1080.0, ax_t.get_ylim()[1] * 0.45, "chromosphere",
        color="darkgoldenrod", fontsize=9, ha="center")
    ax_t.text(
        2400.0, ax_t.get_ylim()[1] * 0.45, "TR",
        color="steelblue", fontsize=9, ha="center")

    titles = {
        "c7": (
            "Model C7 hydrogen ionization degree vs height\n"
            "(Avrett & Loeser 2008, Table 26)"),
        "classical-saha": (
            "Classical pure-H Saha ionization along the Model C7 atmosphere\n"
            "using the local C7 temperature and hydrogen density"),
        "compare": (
            "Model C7, classical Saha, and excitation-enabled CRASH ionization\n"
            "evaluated along the same C7 temperature-density profile"),
    }
    ax_t.set_title(titles[mode])
    ax_t.legend(
        handles=[line_t, line_rho] + ion_lines,
        loc="center left", fontsize=9.2, framealpha=0.9)
    ax_t.grid(True, which="both", ls=":", alpha=0.3)

    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    fig.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("wrote", output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode", choices=("all", "c7", "classical-saha", "compare"),
        default="all", help="figure to generate (default: all three)")
    args = parser.parse_args()

    photo, c7 = parse_c7_tables(CPP)
    h = np.unique(np.concatenate((
        np.linspace(0.0, H_TOP_KM, 1200),
        np.linspace(2119.0, 2147.0, 400),
    )))
    temperature, n_e, n_hi = c7_full_profile(h, photo, c7)
    n_h = n_hi + n_e
    rho = M_H * n_h
    x_c7 = n_e / n_h
    x_saha = saha_ionization_fraction(n_h, temperature)
    t_grid, n_grid, z_grid = load_gamma_grid(
        CRASH_IONIZATION_TABLE, COL["Z"])
    x_crash = gamma_at(
        temperature, n_h, t_grid, n_grid, z_grid)

    modes = ("c7", "classical-saha", "compare") if args.mode == "all" else (args.mode,)
    for mode in modes:
        make_plot(
            mode, h, temperature, rho, x_c7, x_saha, x_crash, OUTPUTS[mode])

    ratio = x_c7 / x_saha
    i_max = int(np.argmax(ratio))
    i_c7_half = int(np.argmin(np.abs(x_c7 - 0.5)))
    i_saha_half = int(np.argmin(np.abs(x_saha - 0.5)))
    print("largest C7/Saha ratio = %.3e at h = %.1f km" %
          (ratio[i_max], h[i_max]))
    print("nearest x=0.5 heights: C7 %.1f km, Saha %.1f km" %
          (h[i_c7_half], h[i_saha_half]))
    crash_abs = np.abs(x_crash - x_saha)
    crash_rel = crash_abs / np.maximum(x_saha, 1.0e-300)
    print("CRASH excitation-on vs analytic Saha: max abs = %.3e, max rel = %.3e" %
          (np.max(crash_abs), np.max(crash_rel)))


if __name__ == "__main__":
    main()
