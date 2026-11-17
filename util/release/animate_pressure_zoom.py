#!/usr/bin/env python3
"""Single-panel movie of p(h) near the top of the column, several runs overlaid.

Reads each run's .gamma_diag snapshots and draws pressure against height above
--h-min, one line per run, on a fixed linear axis. Frames are rendered in
parallel (util/lib/_anim_parallel.py).

Usage:
    python util/release/animate_pressure_zoom.py --run "LABEL=path.txt" [--run ...]
        [--h-min 2100] [--out FILE.mp4] [--fps 20] [--max-frames 200] [--title TEXT]
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), "util", "lib"))
from _anim_parallel import save_frames_parallel
from twall_sweep_diag import IP, load_gamma_diag

STYLES = [dict(color="#2f3640", lw=2.6, ls="-"), dict(color="#c0392b", lw=1.8, ls="-"),
          dict(color="#1f6fb2", lw=1.2, ls="--"), dict(color="#27ae60", lw=1.2, ls=":")]


def _render_frame(k, tmpdir, h, curves, labels, t, ylim, xlim, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5))
    for y, lab, st in zip(curves, labels, STYLES):
        ax.plot(h, y, label=lab, **st)
    ax.set(xlim=xlim, ylim=ylim, xlabel="height [km]", ylabel="pressure [Pa]",
           title=f"{title}   t = {t:6.0f} s")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=120)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", action="append", required=True,
                    help="LABEL=path.txt (repeat; split on the last '=')")
    ap.add_argument("--h-min", type=float, default=2100.0)
    ap.add_argument("--out", default="visualization/model_column/pressure_zoom.mp4")
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--max-frames", type=int, default=200)
    ap.add_argument("--title", default="Pressure near the top")
    a = ap.parse_args()
    runs = []
    for spec in a.run:
        label, path = spec.rsplit("=", 1)
        h, t, d = load_gamma_diag(path)
        runs.append((label, np.asarray(h), t, d[:, :, IP]))
    h = runs[0][1]
    sel = h >= a.h_min
    times = min((r[2] for r in runs), key=lambda x: x[-1])          # overlap window
    if len(times) > a.max_frames:
        times = times[np.unique(np.linspace(0, len(times) - 1, a.max_frames).round().astype(int))]
    picks = [np.abs(r[2][None, :] - times[:, None]).argmin(axis=1) for r in runs]
    allp = np.concatenate([r[3][p][:, sel].ravel() for r, p in zip(runs, picks)])
    pad = 0.03 * (allp.max() - allp.min())
    ylim = (allp.min() - pad, allp.max() + pad)
    xlim = (a.h_min, h[sel][-1] + 0.5)
    frame_args = [(k, h[sel], [r[3][p[k], sel] for r, p in zip(runs, picks)],
                   [r[0] for r in runs], float(times[k]), ylim, xlim, a.title)
                  for k in range(len(times))]
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    save_frames_parallel(_render_frame, frame_args, a.out, a.fps)


if __name__ == "__main__":
    main()
