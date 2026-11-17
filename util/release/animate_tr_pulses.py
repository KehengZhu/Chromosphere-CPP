#!/usr/bin/env python3
"""Side-by-side movie of the transition-region (TR) pulses in two model_column runs.

Made for the hot-wall TR pulses: at T_wall = 100 kK the 7.0-7.5 kK base of the
TR is narrower than one N500 cell, and every time it enters a new cell the hot
layer above it pulses. One column per run, from the face-flux sidecar (cell
state every 1780 steps), zoomed on the top of the domain:
  row 1  T(h), dots at the cell centres, the 7.0-7.5 kK band shaded;
  row 2  V(h);
  row 3  the pulse, V minus its 60 s running median at each cell;
  row 4  the pulse 3 km above the 7 kK isotherm against time, with a cursor.
Frames are rendered in parallel (util/lib/_anim_parallel.py).

Usage:
    python util/release/animate_tr_pulses.py --run "LABEL=path.txt" --run "LABEL=path.txt"
        [--t0 300] [--t1 500] [--out FILE.mp4] [--fps 8]
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
from scipy.ndimage import median_filter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), "util", "lib"))
from _anim_parallel import save_frames_parallel
from plot_outer_bc_and_tr_pulses import HP_WINDOW_S, T_FOOT, foot_height, load_cells


def prepare(label, path, t0, t1):
    c = load_cells(path, t0, t1)
    dt = float(np.median(np.diff(c["t"])))
    win = max(3, int(round(HP_WINDOW_S / dt)) | 1)
    c["hp"] = c["v"] - median_filter(c["v"], size=(win, 1), mode="nearest")
    c["foot"] = np.array([foot_height(c["h"], row, T_FOOT) for row in c["T"]])
    k = [int(np.searchsorted(c["h"], f + 3.0)) for f in c["foot"]]
    c["trace"] = c["hp"][np.arange(len(k)), k]
    c["label"] = label
    return c


def _render_frame(k, tmpdir, panels, t, lims, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(4, len(panels), figsize=(7.2 * len(panels), 10.5),
                           sharex="row", sharey="row",
                           gridspec_kw=dict(height_ratios=[1.2, 1, 1, 0.9]))
    for j, p in enumerate(panels):
        h, sel = p["h"], p["sel"]
        ax[0, j].plot(h[sel], p["T"][sel], "o-", ms=3, lw=1.0, color="k")
        ax[0, j].axhspan(7000, 7500, color="C1", alpha=0.3, lw=0)
        ax[0, j].set(yscale="log", ylim=lims["T"], title=f"{p['label']}   (t = {t:.1f} s)")
        ax[1, j].plot(h[sel], p["v"][sel], lw=1.2, color="C3")
        ax[1, j].set(ylim=lims["v"])
        ax[2, j].plot(h[sel], p["hp"][sel], lw=1.2, color="C0")
        ax[2, j].axhline(0, color="0.6", lw=0.6)
        ax[2, j].set(ylim=lims["hp"], xlabel="height [km]")
        for r in range(3):
            ax[r, j].axvline(p["foot"], color="C1", lw=0.8, ls=":")
            ax[r, j].set_xlim(lims["h"])
        ax[3, j].plot(p["t_all"], p["trace_all"], lw=0.8, color="C0")
        ax[3, j].axvline(t, color="k", lw=1.0)
        ax[3, j].set(ylim=lims["hp"], xlabel="t [s]")
    ax[0, 0].set_ylabel("T [K]\n(dots = cell centres)")
    ax[1, 0].set_ylabel("velocity [m/s]")
    ax[2, 0].set_ylabel("pulse = velocity minus\n60 s running median [m/s]")
    ax[3, 0].set_ylabel("pulse 3 km above\nthe 7 kK isotherm [m/s]")
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=100)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", action="append", required=True,
                    help="LABEL=path.txt (repeat; split on the last '=')")
    ap.add_argument("--t0", type=float, default=300.0)
    ap.add_argument("--t1", type=float, default=500.0)
    ap.add_argument("--out", default="visualization/model_column/tr_pulses_N500_vs_N1000_100kK.mp4")
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--title", default="Transition-region pulses at a 100 kK wall")
    a = ap.parse_args()
    runs = [prepare(*spec.rsplit("=", 1), a.t0, a.t1) for spec in a.run]
    times = runs[0]["t"]                                 # the coarsest capture grid
    lo = min(r["foot"].min() for r in runs) - 2.0
    hi = max(r["h"][-1] for r in runs) + 0.5
    def span(key, q):
        vals = np.concatenate([r[key][:, (r["h"] >= lo)].ravel() for r in runs])
        return np.percentile(vals, [q, 100 - q])
    v0, v1 = span("v", 0.2)
    hp = np.abs(span("hp", 0.2)).max() * 1.1
    lims = dict(h=(lo, hi), T=(5.0e3, 1.2e5), v=(v0 - 50, v1 + 50), hp=(-hp, hp))
    frame_args = []
    for k, t in enumerate(times):
        panels = []
        for r in runs:
            i = int(np.abs(r["t"] - t).argmin())
            sel = r["h"] >= lo
            panels.append(dict(label=r["label"], h=r["h"], sel=sel, T=r["T"][i],
                               v=r["v"][i], hp=r["hp"][i], foot=r["foot"][i],
                               t_all=r["t"], trace_all=r["trace"]))
        frame_args.append((k, panels, float(t), lims, a.title))
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    save_frames_parallel(_render_frame, frame_args, a.out, a.fps)


if __name__ == "__main__":
    main()
