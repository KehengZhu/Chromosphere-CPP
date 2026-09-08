#!/usr/bin/env python3
"""Wall-temperature response of `model_column` under the unanchored outer pressure.

Four runs, the 2x2 of {release outer closure, `ISO_OUTER_P_NEUMANN=2` unanchored
`dp/ds = -rho g`} x {22 kK, 100 kK conduction wall}, reduced to the four quantities
that decide whether the unanchored closure can be used for a wall-temperature study:

    (a) V_top(t)        the evaporation signal itself
    (b) p_top(t)/p_top(0)   the pressure level the release pins and this closure does not
    (c) M(t)/M(0)       column mass budget
    (d) V(h) quasi-steady   mean over the last `--avg-window` seconds

Panel (b) is the point of the figure: with no external pressure datum the level is
free, so any evaporation number measured under this closure is quoted against a
moving reference. Panels (a) and (d) say how much of the `T_wall` response survives.

Writes the figure and a JSON of the tabulated metrics. Loaders and the
`.gamma_diag` column layout are imported from `util/twall_sweep_diag.py`, so this
figure and the sweep figures always read the runs the same way.

Usage:
  .venv/bin/python util/outer_bc_twall_diag.py [--avg-window 500]
      [--out visualization/model_column/outer_bc_twall_response.png]
      [--json outputs/model_column/outer_bc_swap/outer_bc_twall_metrics.json]
"""
import argparse
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from twall_sweep_diag import (IRHO, IP, IT, IV, cell_weights, console_ok,
                              load_gamma_diag)

# (label, path, colour, linestyle) — colour carries T_wall, dash carries the
# closure, the same encoding as `util/animate_outer_bc_compare.py --preset twall`.
RUNS = [
    ("release BC, 22 kK", "outputs/model_column/twall_sweep/twall_22kK_cfl0.50.txt",
     "#1f6fb2", "-"),
    ("release BC, 100 kK", "outputs/model_column/twall_sweep/twall_100kK_cfl0.25.txt",
     "#c0392b", "-"),
    ("unanchored $dp/ds=-\\rho g$, 22 kK",
     "outputs/model_column/outer_bc_swap/hse_22kK_cfl0.25.txt", "#1f6fb2", "--"),
    ("unanchored $dp/ds=-\\rho g$, 100 kK",
     "outputs/model_column/outer_bc_swap/hse_100kK_cfl0.25.txt", "#c0392b", "--"),
]


def analyse(label, path, avg_window):
    h, times, data = load_gamma_diag(path)
    term, fallbacks = console_ok(path)
    ds = cell_weights(h)
    mass = (data[:, :, IRHO] * ds).sum(axis=1)
    v_top = data[:, -1, IV]
    p_top = data[:, -1, IP]
    window = times >= times[-1] - avg_window
    # Drift of the quasi-steady velocity ACROSS the averaging window, the same
    # non-stationarity measure the wall-temperature sweep recap reports. Without an
    # anchored pressure the column has no steady state at all, so the window mean is
    # a snapshot of a moving solution and must be quoted with this number.
    first, second = np.array_split(np.flatnonzero(window), 2)
    drift = float(v_top[second].mean() / v_top[first].mean() - 1.0)
    return dict(label=label, path=path, term=term, fallbacks=fallbacks,
                h=h, times=times, v_top=v_top, p_top=p_top, mass=mass,
                v_of_h=data[window, :, IV].mean(axis=0),
                v_top_qs=float(v_top[window].mean()),
                v_top_end=float(v_top[-1]),
                v_top_drift_window=drift,
                p_drift=float(p_top[-1] / p_top[0] - 1.0),
                mass_drift=float(mass[-1] / mass[0] - 1.0),
                t_end=float(times[-1]),
                T_top_end=float(data[-1, -1, IT]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--avg-window", type=float, default=500.0,
                    help="quasi-steady average window [s] at the end of the run")
    ap.add_argument("--out", default="visualization/model_column/"
                                     "outer_bc_twall_response.png")
    ap.add_argument("--json", default="outputs/model_column/outer_bc_swap/"
                                      "outer_bc_twall_metrics.json")
    args = ap.parse_args()

    runs = [analyse(lab, p, args.avg_window) for lab, p, _, _ in RUNS]
    for r in runs:
        print(f"{r['label']:36s} term={r['term']} t_end={r['t_end']:7.1f} "
              f"V_top(qs)={r['v_top_qs']:+10.3f} m/s  p drift={r['p_drift']*100:+7.2f} %"
              f"  mass drift={r['mass_drift']*100:+8.4f} %"
              f"  V drift/window={r['v_top_drift_window']*100:+7.2f} %")

    # Effective two-point exponent of the wall-temperature response, V ~ T_wall^n,
    # per closure. Two points cannot establish a power law; this is the exponent the
    # sweep's fitted one (2.91 over 15-100 kK) is compared against, nothing more.
    pairs = {"release": (runs[0], runs[1]), "unanchored": (runs[2], runs[3])}
    exponents = {}
    for name, (cool, hot) in pairs.items():
        n = float(np.log(hot["v_top_qs"] / cool["v_top_qs"]) / np.log(100.0 / 22.0))
        exponents[name] = n
        print(f"{name:>12s}: V(100 kK)/V(22 kK) = "
              f"{hot['v_top_qs'] / cool['v_top_qs']:7.2f}  ->  n = {n:.3f}")
    for name, (_, hot) in pairs.items():
        if name == "release":
            continue
        print(f"suppression vs release: 22 kK {runs[0]['v_top_qs']/runs[2]['v_top_qs']:.1f}x"
              f"   100 kK {runs[1]['v_top_qs']/runs[3]['v_top_qs']:.1f}x")

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 9,
                         "axes.linewidth": 0.6, "figure.dpi": 150,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.4))

    for r, (_, _, colour, ls) in zip(runs, RUNS):
        style = dict(color=colour, ls=ls, lw=1.6, label=r["label"])
        ax[0, 0].plot(r["times"], r["v_top"], **style)
        ax[0, 1].plot(r["times"], r["p_top"] / r["p_top"][0], **style)
        ax[1, 0].plot(r["times"], r["mass"] / r["mass"][0], **style)
        ax[1, 1].plot(r["h"], r["v_of_h"], **style)

    ax[0, 0].set(xlabel="t [s]", ylabel="$V_{\\rm top}$ [m/s]",
                 title="(a) top-cell velocity")
    ax[0, 0].set_yscale("symlog", linthresh=1.0)
    ax[0, 0].axhline(0.0, c="k", lw=0.5)
    ax[0, 1].set(xlabel="t [s]", ylabel="$p_{\\rm top}(t)/p_{\\rm top}(0)$",
                 title="(b) outer pressure level — pinned by the release BC, free here")
    ax[0, 1].axhline(1.0, c="k", lw=0.5)
    ax[1, 0].set(xlabel="t [s]", ylabel="$M(t)/M(0)$", title="(c) column mass")
    ax[1, 0].axhline(1.0, c="k", lw=0.5)
    ax[1, 1].set(xlabel="height [km]", ylabel="$V$ [m/s]",
                 title=f"(d) quasi-steady $V(h)$ (mean over last {args.avg_window:.0f} s)")
    ax[1, 1].set_yscale("symlog", linthresh=1.0)
    ax[1, 1].axhline(0.0, c="k", lw=0.5)
    for a in ax.ravel():
        a.grid(True, alpha=0.3, which="both")

    # One figure-level legend under all four panels: inside panel (a) it sits on top
    # of the 100 kK curves, which are the ones the panel exists to show.
    handles, legend_labels = ax[0, 0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="lower center", ncol=2, frameon=False,
               fontsize=9)

    fig.suptitle("model_column wall-temperature response with and without an "
                 "externally anchored outer pressure — N=500/R4, double state",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0.075, 1, 0.955])
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out)
    print("wrote", args.out)

    if args.json:
        os.makedirs(os.path.dirname(args.json) or ".", exist_ok=True)
        drop = ("h", "times", "v_top", "p_top", "mass", "v_of_h")
        with open(args.json, "w") as fh:
            json.dump({"avg_window_s": args.avg_window,
                       "runs": [{k: v for k, v in r.items() if k not in drop}
                                for r in runs]}, fh, indent=2)
        print("wrote", args.json)


if __name__ == "__main__":
    main()
