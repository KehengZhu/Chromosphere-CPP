#!/usr/bin/env python3
"""Three figures comparing the release scheme before and after the MC3 branch fix
and with the experimental hydrostatic reconstruction (ISO_HSE_RECON=1).

Reads the A/B/C runs of docs/studies/numerics/hydrostatic_reconstruction_and_mc3_fix.md
(A = pre-fix binary, B = MC3 fixed, C = B + ISO_HSE_RECON=1) from --runs and writes
to --out:

  hse_recon_static_residual.pdf  column at rest, conduction off, 20 s: t = 0
                                 max|F^rho| [kg m^-2 s^-1] and 20 s max|V| [m/s],
                                 500 and 1000 cells (static_tt0_{A,B,C}_N*_20s);
  hse_recon_evaporation.pdf      conduction-driven 1000 s runs (prod_{A,B,C}_N*_1000s):
                                 mean V in 2130-2150 km [m/s] and relative column
                                 mass change [1e-4] against time [s];
  hse_recon_evaporation_bars.pdf the same 500-cell runs at t = 1000 s, as bars:
                                 mean V in 2130-2150 km [m/s] and the dissipative
                                 share |f_diff/f_total| of the base-face mass flux [%].

Usage:
  util/release/plot_hse_recon_comparison.py [--runs DIR] [--out DIR]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LABEL = {"A": "original", "B": "limiter fixed", "C": "limiter fixed + hydrostatic"}
COLOR = {"A": "#4c72b0", "B": "#dd8452", "C": "#55a868"}
WINDOW_KM = (2130.0, 2150.0)
BASE_KM = 1600.0


def frames(path: Path) -> list[tuple[float, np.ndarray]]:
    """All `# t = ...` frames of a sidecar as (time [s], rows) pairs."""
    out, rows, t = [], None, None
    with open(path) as handle:
        for line in handle:
            if line.startswith("# t = "):
                if rows is not None:
                    out.append((t, np.array(rows, float)))
                t, rows = float(line.split()[3]), []
            elif rows is not None and not line.startswith("#") and line.strip():
                rows.append(line.split())
    if rows:
        out.append((t, np.array(rows, float)))
    return out


def static_residual(runs: Path, out: Path) -> None:
    stat = {}
    for n in (500, 1000):
        for k in "ABC":
            run = runs / f"static_tt0_{k}_N{n}_20s.txt"
            face = frames(Path(str(run) + ".faceflux"))[0][1]
            cells = frames(Path(str(run) + ".gamma_diag"))[-1][1]
            stat[k, n] = (np.max(np.abs(face[:, 16])), np.max(np.abs(cells[:, 1])))
    fig, ax = plt.subplots(1, 2, figsize=(8, 3.5))
    x, w = np.arange(2), 0.26
    for j, k in enumerate("ABC"):
        for m in (0, 1):
            ax[m].bar(x + (j - 1) * w, [stat[k, n][m] for n in (500, 1000)], w,
                      color=COLOR[k], label=LABEL[k])
    titles = ["spurious mass flux at t = 0\n[kg m$^{-2}$ s$^{-1}$]",
              "largest spurious speed after 20 s\n[m/s]"]
    for m, title in enumerate(titles):
        ax[m].set_yscale("log")
        ax[m].set_xticks(x)
        ax[m].set_xticklabels(["500 cells", "1000 cells"])
        ax[m].set_title(title, fontsize=10)
        ax[m].set_xlim(-0.5, 1.65)
        ax[m].text(1 + w, stat["C", 1000][m] * 6, "round-off", ha="center", va="bottom",
                   fontsize=8, color=COLOR["C"], zorder=5,
                   bbox=dict(fc="white", ec="none", pad=0.3))
    handles, labels = ax[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=3, fontsize=9, frameon=False, loc="upper center")
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(out / "hse_recon_static_residual.pdf")
    plt.close(fig)


def evaporation(runs: Path, out: Path) -> None:
    fig, ax = plt.subplots(1, 2, figsize=(8, 3.2))
    for n, style in ((500, "-"), (1000, "--")):
        for k in "ABC":
            path = runs / f"prod_{k}_N{n}_1000s.txt.gamma_diag"
            with open(path) as handle:
                handle.readline()
                tops = np.array(handle.readline().split(), float)   # cell upper faces [km]
            edges = np.concatenate([[BASE_KM], tops])
            dz = np.diff(edges)
            centre = 0.5 * (edges[1:] + edges[:-1])
            window = (centre >= WINDOW_KM[0]) & (centre <= WINDOW_KM[1])
            series = frames(path)
            t = np.array([f[0] for f in series])
            v = np.array([f[1][window, 1].mean() for f in series])
            mass = np.array([np.sum(f[1][:, 0] * dz) for f in series])
            ax[0].plot(t, v, style, color=COLOR[k], label=f"{LABEL[k]}, {n}")
            ax[1].plot(t, (mass / mass[0] - 1) * 1e4, style, color=COLOR[k])
    ax[0].set_xlabel("time [s]")
    ax[0].set_ylabel("upflow speed, 2130-2150 km [m/s]")
    ax[1].set_xlabel("time [s]")
    ax[1].set_ylabel("column mass change [10$^{-4}$]")
    ax[0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(out / "hse_recon_evaporation.pdf")
    plt.close(fig)


def evaporation_bars(runs: Path, out: Path) -> None:
    speed, share = {}, {}
    for k in "ABC":
        path = runs / f"prod_{k}_N500_1000s.txt.gamma_diag"
        with open(path) as handle:
            handle.readline()
            tops = np.array(handle.readline().split(), float)
        # window on the listed upper-face heights, as pressure_reconstruction_diag.py
        window = (tops >= WINDOW_KM[0]) & (tops <= WINDOW_KM[1])
        speed[k] = frames(path)[-1][1][window, 1].mean()
        base = frames(runs / f"prod_{k}_N500_1000s.txt.faceflux")[-1][1][0]
        share[k] = abs(base[15] / base[16]) * 100.0          # f_diff / f_total
    fig, ax = plt.subplots(1, 2, figsize=(8, 3.2))
    x = np.arange(3)
    for m, (vals, fmt, title) in enumerate((
            (speed, "{:.2f}", "evaporation speed at 1000 s\n[m/s]"),
            (share, "{:.1f} %", "numerical error at the bottom face\n[% of the mass flux]"))):
        heights = [vals[k] for k in "ABC"]
        ax[m].bar(x, heights, 0.6, color=[COLOR[k] for k in "ABC"])
        for xi, h in zip(x, heights):
            ax[m].text(xi, h, fmt.format(h), ha="center", va="bottom", fontsize=10)
        ax[m].set_xticks(x)
        ax[m].set_xticklabels(["original", "limiter\nfixed", "+ hydrostatic"])
        ax[m].set_title(title, fontsize=10)
        ax[m].set_ylim(0, 1.2 * max(heights))
        ax[m].spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / "hse_recon_evaporation_bars.pdf")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--runs", type=Path, default=Path("outputs/model_column/hse_recon"))
    parser.add_argument("--out", type=Path, default=Path("visualization/model_column"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 11})
    static_residual(args.runs, args.out)
    evaporation(args.runs, args.out)
    evaporation_bars(args.runs, args.out)


if __name__ == "__main__":
    main()
