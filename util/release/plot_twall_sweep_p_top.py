#!/usr/bin/env python3
"""Absolute pressure in the top of the column for every wall of the T_wall sweep.

One snapshot, one line per imposed conduction-wall temperature, absolute p in Pa
on a linear axis, restricted to h >= --hmin (default 2100 km). Over the full
column p falls 6.6x and the 13 walls sit within a line width of one another; over
the top ~50 km the whole wall-to-wall spread is resolved, so no normalization is
needed. Each run's transition-region base (first cell above --t-tr, default
7.5 kK) is marked with a dot: that is where its p(h) turns from the steep
chromospheric slope to the nearly flat hot-TR slope.

Run discovery, acceptance and colours are those of util/release/twall_sweep_diag.py and
util/release/animate_twall_sweep.py, so this frame matches the overlay movies.

Usage:
    python util/release/plot_twall_sweep_p_top.py [--time 2000] [--hmin 2100]
        [--run T_WALL_K=path] [--out FILE.png]
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from animate_twall_sweep import collect_entries, T_WALL_RELEASE
from twall_sweep_diag import IP, IT, console_ok, load_gamma_diag, wall_colors


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="outputs/model_column/twall_sweep")
    ap.add_argument("--cfl", default="0.50")
    ap.add_argument("--run", action="append", default=[],
                    help="explicit T_wall_K=path entry, added to the glob")
    ap.add_argument("--time", type=float, default=2000.0,
                    help="snapshot time [s]; the nearest snapshot is used")
    ap.add_argument("--hmin", type=float, default=2100.0,
                    help="lower edge of the plotted height range [km]")
    ap.add_argument("--t-tr", type=float, default=7.5e3,
                    help="temperature marking the TR base [K]")
    ap.add_argument("--out", default="visualization/model_column/"
                                     "twall_sweep_N500_p_abs_top_t2000.png")
    args = ap.parse_args()

    runs = []
    for t_wall, path in collect_entries(args.dir, args.cfl, args.run):
        if console_ok(path)[0] != "end_time":
            print(f"{t_wall/1e3:>5g} kK  EXCLUDED: did not reach end time")
            continue
        h, times, data = load_gamma_diag(path)
        k = int(np.argmin(np.abs(times - args.time)))
        runs.append((t_wall, h, times[k], data[k, :, IP], data[k, :, IT]))
    if not runs:
        raise SystemExit("no accepted runs")
    t_snap = runs[0][2]
    worst = max(abs(r[2] - t_snap) for r in runs)

    col, _, _ = wall_colors([r[0] for r in runs])
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for t_wall, h, _, p, T in runs:
        sel = h >= args.hmin
        lw = 2.4 if t_wall == T_WALL_RELEASE else 1.2
        ax.plot(h[sel], p[sel], color=col[t_wall], lw=lw, label=f"{t_wall/1e3:g}")
        hot = np.nonzero(T > args.t_tr)[0]
        if hot.size and h[hot[0]] >= args.hmin:
            ax.plot(h[hot[0]], p[hot[0]], "o", color=col[t_wall], ms=5, mec="k",
                    mew=0.5)
    ax.set_xlim(args.hmin, max(r[1].max() for r in runs))
    ax.set_xlabel("height [km]")
    ax.set_ylabel("p [Pa]")
    ax.set_title(f"Absolute pressure above {args.hmin:g} km, t = {t_snap:.0f} s  "
                 f"(dot: TR base, first cell with T > {args.t_tr/1e3:g} kK)",
                 fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(title=r"$T_{\rm wall}$ [kK]", ncol=2, fontsize=8, title_fontsize=9)
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"{len(runs)} runs, snapshot t = {t_snap:.1f} s (worst mismatch "
          f"{worst:.3f} s); wrote {args.out}")


if __name__ == "__main__":
    main()
