#!/usr/bin/env python3
"""Compare full CRASH ionization tables with excitation off/on.

The complete Model C7 (T, n_H) path is overlaid on every panel and labeled at
representative atmospheric heights.  The third panel shows the fractional
reduction in ionization caused by bound excitation.
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import numpy as np

from plot_c7_gamma_profile import (
    COL,
    CPP,
    H_TOP_KM,
    c7_full_profile,
    parse_c7_tables,
)


TABLE_OFF = "outputs/eos_gamma/gamma_hydrogen.dat"
TABLE_ON = "outputs/eos_gamma/gamma_hydrogen_excitation.dat"
OUTPUT = "visualization/eos_gamma/crash_excitation_full_table_c7.png"


def load_ionization_grid(path):
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] != 11 or not np.all(np.isfinite(data)):
        raise RuntimeError(f"{path}: expected a finite 11-column CRASH table")
    temperatures = np.unique(data[:, COL["T"]])
    densities = np.unique(data[:, COL["Na"]])
    expected = len(temperatures) * len(densities)
    if len(data) != expected:
        raise RuntimeError(f"{path}: table is not rectangular")
    if not (
        np.array_equal(data[:, COL["T"]], np.tile(temperatures, len(densities)))
        and np.array_equal(data[:, COL["Na"]], np.repeat(densities, len(temperatures)))
    ):
        raise RuntimeError(f"{path}: unexpected table ordering")
    ionization = data[:, COL["Z"]].reshape(len(densities), len(temperatures))
    if np.any(ionization < 0.0) or np.any(ionization > 1.0):
        raise RuntimeError(f"{path}: ionization degree lies outside [0,1]")
    return temperatures, densities, ionization


def add_c7_track(ax, temperature, density, label=True):
    line, = ax.plot(
        temperature, density, color="#111111", lw=2.4, zorder=6,
        label="Model C7 track" if label else None)
    line.set_path_effects([
        path_effects.Stroke(linewidth=4.4, foreground="white"),
        path_effects.Normal(),
    ])
    ax.scatter(
        temperature[[0, -1]], density[[0, -1]], s=28,
        color="#111111", edgecolor="white", linewidth=0.8, zorder=7)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=OUTPUT, help="output PNG path")
    args = parser.parse_args()

    t_off, n_off, x_off = load_ionization_grid(TABLE_OFF)
    t_on, n_on, x_on = load_ionization_grid(TABLE_ON)
    if not (np.array_equal(t_off, t_on) and np.array_equal(n_off, n_on)):
        raise RuntimeError("excitation-off/on CRASH tables use different grids")

    suppression = 100.0 * (1.0 - x_on / np.maximum(x_off, 1.0e-300))
    suppression = np.clip(suppression, 0.0, None)

    photo, c7 = parse_c7_tables(CPP)
    h = np.linspace(0.0, H_TOP_KM, 1200)
    t_c7, n_e, n_hi = c7_full_profile(h, photo, c7)
    n_c7 = n_e + n_hi

    fig, axes = plt.subplots(
        1, 3, figsize=(16.2, 5.7), sharex=True, sharey=True,
        constrained_layout=True)
    ion_norm = LogNorm(vmin=1.0e-10, vmax=1.0)
    meshes = [
        axes[0].pcolormesh(
            t_off, n_off, x_off, shading="auto", cmap="viridis",
            norm=ion_norm, rasterized=True),
        axes[1].pcolormesh(
            t_on, n_on, x_on, shading="auto", cmap="viridis",
            norm=ion_norm, rasterized=True),
        axes[2].pcolormesh(
            t_on, n_on, suppression, shading="auto", cmap="magma",
            norm=Normalize(vmin=0.0, vmax=45.0), rasterized=True),
    ]

    titles = [
        "CRASH ionization: excitation OFF",
        "CRASH ionization: excitation ON",
        "Ionization reduction from excitation",
    ]
    for ax, title in zip(axes, titles):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(t_off[0], t_off[-1])
        ax.set_ylim(n_off[0], n_off[-1])
        ax.set_title(title, fontsize=12.5)
        ax.set_xlabel(r"temperature  $T$  [K]")
        ax.grid(True, which="major", color="white", ls=":", alpha=0.18)
        add_c7_track(ax, t_c7, n_c7)
    axes[0].set_ylabel(r"hydrogen nuclei density  $n_{\rm H}$  [m$^{-3}$]")

    cbar0 = fig.colorbar(
        meshes[0], ax=axes[:2], location="bottom", shrink=0.72, pad=0.08,
        aspect=36)
    cbar0.set_label(r"CRASH mean ionization degree  $\bar{Z}$")
    cbar1 = fig.colorbar(
        meshes[2], ax=axes[2], location="bottom", shrink=0.90, pad=0.08,
        aspect=18)
    cbar1.set_label(r"reduction  $100(1-\bar{Z}_{\rm on}/\bar{Z}_{\rm off})$  [%]")

    tt, nn = np.meshgrid(t_on, n_on)
    contours = axes[2].contour(
        tt, nn, suppression, levels=(1.0, 10.0, 30.0),
        colors="white", linewidths=(0.8, 1.0, 1.2), alpha=0.9)
    axes[2].clabel(contours, fmt=lambda value: f"{value:.0f}%", fontsize=8)

    landmarks = [
        ("photosphere", 0.0, (1.35, 1.35)),
        (r"$T_{\min}$", 560.0, (1.45, 0.62)),
        ("upper chrom.", 1800.0, (1.55, 1.60)),
        ("TR base", 2153.0, (1.65, 0.55)),
        ("upper TR", H_TOP_KM, (0.34, 1.80)),
    ]
    for label, target_h, offset in landmarks:
        i = int(np.argmin(np.abs(h - target_h)))
        axes[2].annotate(
            label, xy=(t_c7[i], n_c7[i]),
            xytext=(t_c7[i] * offset[0], n_c7[i] * offset[1]),
            fontsize=8.2, color="#111111", ha="left", va="center",
            arrowprops=dict(arrowstyle="-", color="#222222", lw=0.8),
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.82),
            zorder=8)

    max_suppression = float(np.max(suppression))
    i_max = np.unravel_index(np.argmax(suppression), suppression.shape)

    # Bilinear interpolation on log coordinates is sufficient for reporting
    # how far the C7 track lies from the excitation-sensitive region.
    from plot_c7_gamma_profile import gamma_at
    x_off_c7 = gamma_at(t_c7, n_c7, t_off, n_off, x_off)
    x_on_c7 = gamma_at(t_c7, n_c7, t_on, n_on, x_on)
    c7_suppression = 100.0 * (
        1.0 - x_on_c7 / np.maximum(x_off_c7, 1.0e-300))
    c7_max = float(np.max(np.abs(c7_suppression)))
    axes[2].text(
        0.03, 0.035,
        "Model C7: maximum excitation effect\n"
        rf"$|\Delta \bar{{Z}}/\bar{{Z}}|={c7_max:.2e}\%$",
        transform=axes[2].transAxes, fontsize=9.0, color="#111111",
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#555555", alpha=0.92),
        zorder=9)
    axes[0].legend(loc="lower right", framealpha=0.92, fontsize=9)

    fig.suptitle(
        "Full CRASH pure-hydrogen EOS table: effect of bound excitation\n"
        "Model C7 occupies a low-density path outside the excitation-sensitive regime",
        fontsize=15)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=180, bbox_inches="tight")
    plt.close(fig)

    print("wrote", args.output)
    print("full-table maximum reduction = %.3f%% at T=%.6g K, n_H=%.6g m^-3" %
          (max_suppression, t_on[i_max[1]], n_on[i_max[0]]))
    print("Model C7 maximum absolute relative effect = %.3e%%" % c7_max)


if __name__ == "__main__":
    main()
