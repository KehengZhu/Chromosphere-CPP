#!/usr/bin/env python3
"""Render an MP4 of the time evolution from chromo_main output.txt."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plot_output import read_frames, primitives
from _anim_parallel import save_frames_parallel


# (row, col, prims_index, panel_title, log_yscale)
_PANELS = [
    (0, 0, 1, r"$n_n$ (m$^{-3}$)", True),
    (0, 1, 3, r"$u$ (m/s)",        False),
    (0, 2, 5, r"$P_n$ (Pa)",       True),
    (0, 3, 7, r"$T_n$ (K)",        False),
    (1, 0, 0, r"$n_i$ (m$^{-3}$)", True),
    (1, 1, 2, r"$v$ (m/s)",        False),
    (1, 2, 4, r"$P_i$ (Pa)",       True),
    (1, 3, 6, r"$T_i$ (K)",        False),
]

# Number-density-weighted bulk velocity (n_n u + n_i v)/(n_n + n_i).
# primitives() order is (ni, nn, v, u, ...): ni=0, nn=1, v=2, u=3.
_MEAN_V_TITLE = r"$(n_n u + n_i v)/(n_n + n_i)$  (m/s)"


def _mean_velocity(prims_k):
    ni, nn, v, u = prims_k[0], prims_k[1], prims_k[2], prims_k[3]
    return (nn * u + ni * v) / (nn + ni)


def _render_frame(k, tmpdir, xx, prims_k, t, step, nF, ylims):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(15, 9.5))
    gs = fig.add_gridspec(3, 4)
    fig.suptitle(
        f"Chromosphere column  |  t = {t:6.2f} s, step {step:5d}  ({k+1}/{nF})"
    )
    for row, col, idx, title, log in _PANELS:
        ax = fig.add_subplot(gs[row, col])
        ax.plot(xx, prims_k[idx], lw=1.2)
        ax.set_title(title)
        ax.set_xlim(xx.min(), xx.max())
        ax.set_ylim(*ylims[title])
        if log:
            ax.set_yscale("log")
        ax.grid(True, alpha=0.3)
        if row == 1:
            ax.set_xlabel("height (km)")

    # Bottom row: number-density-weighted bulk velocity, spanning all columns.
    axm = fig.add_subplot(gs[2, :])
    axm.plot(xx, _mean_velocity(prims_k), lw=1.2, color="C3")
    axm.set_title(_MEAN_V_TITLE)
    axm.set_xlim(xx.min(), xx.max())
    axm.set_ylim(*ylims[_MEAN_V_TITLE])
    axm.grid(True, alpha=0.3)
    axm.set_xlabel("height (km)")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=120)
    plt.close(fig)


def main():
    in_path  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "build", "output.txt")
    out_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "..", "visualization", "chromo_evolution.mp4")
    fps      = int(sys.argv[3]) if len(sys.argv) > 3 else 15
    os.makedirs(os.path.join(HERE, "..", "visualization"), exist_ok=True)

    xx, frames = read_frames(in_path)
    nF = len(frames)
    print(f"read {nF} frames, ns = {xx.size}, t in [{frames[0][0]:.1f}, {frames[-1][0]:.1f}] s")

    prims = [primitives(fr[2]) for fr in frames]

    def lim(arr_list, log=False):
        a = np.concatenate(arr_list)
        if log:
            a = a[a > 0]
        lo, hi = float(np.min(a)), float(np.max(a))
        if log:
            return lo / 1.2, hi * 1.2
        pad = 0.05 * (hi - lo if hi > lo else max(abs(hi), 1.0))
        return lo - pad, hi + pad

    ylims = {title: lim([p[idx] for p in prims], log=log)
             for _, _, idx, title, log in _PANELS}
    ylims[_MEAN_V_TITLE] = lim([_mean_velocity(p) for p in prims], log=False)

    save_frames_parallel(
        _render_frame,
        [(k, xx, prims[k], frames[k][0], frames[k][1], nF, ylims) for k in range(nF)],
        out_path, fps,
    )


if __name__ == "__main__":
    main()
