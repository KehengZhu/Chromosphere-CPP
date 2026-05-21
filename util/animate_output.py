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


def _render_frame(k, tmpdir, xx, prims_k, t, step, nF, ylims):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 4, figsize=(15, 7))
    fig.suptitle(
        f"Chromosphere — Model C7  |  t = {t:6.2f} s, step {step:5d}  ({k+1}/{nF})"
    )
    for row, col, idx, title, log in _PANELS:
        ax = axes[row, col]
        ax.plot(xx, prims_k[idx], lw=1.2)
        ax.set_title(title)
        ax.set_xlim(xx.min(), xx.max())
        ax.set_ylim(*ylims[title])
        if log:
            ax.set_yscale("log")
        ax.grid(True, alpha=0.3)
    for ax in axes[1, :]:
        ax.set_xlabel("height (km)")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=120)
    plt.close(fig)


def main():
    in_path  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "build", "output.txt")
    out_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "visualization", "model_c7_evolution.mp4")
    fps      = int(sys.argv[3]) if len(sys.argv) > 3 else 15
    os.makedirs(os.path.join(HERE, "visualization"), exist_ok=True)

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

    save_frames_parallel(
        _render_frame,
        [(k, xx, prims[k], frames[k][0], frames[k][1], nF, ylims) for k in range(nF)],
        out_path, fps,
    )


if __name__ == "__main__":
    main()
