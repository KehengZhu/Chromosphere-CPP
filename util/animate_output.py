#!/usr/bin/env python3
"""Render an MP4 of the time evolution from build/output.txt."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from plot_output import read_frames, primitives


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    in_path  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "build", "output.txt")
    out_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "model_c7_evolution.mp4")
    fps      = int(sys.argv[3]) if len(sys.argv) > 3 else 15

    xx, frames = read_frames(in_path)
    nF = len(frames)
    print(f"read {nF} frames, ns = {xx.size}, t in [{frames[0][0]:.1f}, {frames[-1][0]:.1f}] s")

    # Pre-compute primitives for every frame so axis limits are stable across the run.
    prims = [primitives(fr[2]) for fr in frames]
    # Stack each field across frames, then pick global min/max with a 5% margin.
    def lim(arr_list, log=False):
        a = np.concatenate(arr_list)
        if log:
            a = a[a > 0]
        lo, hi = float(np.min(a)), float(np.max(a))
        if log:
            return lo / 1.2, hi * 1.2
        pad = 0.05 * (hi - lo if hi > lo else max(abs(hi), 1.0))
        return lo - pad, hi + pad

    nn_all = [p[1] for p in prims]
    u_all  = [p[3] for p in prims]
    pn_all = [p[5] for p in prims]
    Tn_all = [p[7] for p in prims]
    ni_all = [p[0] for p in prims]
    v_all  = [p[2] for p in prims]
    pi_all = [p[4] for p in prims]
    Ti_all = [p[6] for p in prims]

    fig, axes = plt.subplots(2, 4, figsize=(15, 7))
    title = fig.suptitle("")

    panel_specs = [
        (axes[0, 0], nn_all, r"$n_n$ (m$^{-3}$)", "log"),
        (axes[0, 1], u_all,  r"$u$ (m/s)",        None),
        (axes[0, 2], pn_all, r"$P_n$ (Pa)",       "log"),
        (axes[0, 3], Tn_all, r"$T_n$ (K)",        None),
        (axes[1, 0], ni_all, r"$n_i$ (m$^{-3}$)", "log"),
        (axes[1, 1], v_all,  r"$v$ (m/s)",        None),
        (axes[1, 2], pi_all, r"$P_i$ (Pa)",       "log"),
        (axes[1, 3], Ti_all, r"$T_i$ (K)",        None),
    ]
    lines = []
    for ax, series, ttl, ys in panel_specs:
        (ln,) = ax.plot(xx, series[0], lw=1.2)
        ax.set_title(ttl)
        ax.set_xlim(xx.min(), xx.max())
        ax.set_ylim(*lim(series, log=(ys == "log")))
        if ys:
            ax.set_yscale(ys)
        ax.grid(True, alpha=0.3)
        lines.append(ln)
    for ax in axes[1, :]:
        ax.set_xlabel("height (km)")

    fig.tight_layout(rect=[0, 0, 1, 0.94])

    def update(k):
        t, step, _ = frames[k]
        ni, nn, v, u, p_i, p_n, Ti, Tn = prims[k]
        for ln, y in zip(lines, [nn, u, p_n, Tn, ni, v, p_i, Ti]):
            ln.set_ydata(y)
        title.set_text(f"Chromosphere — Model C7  |  t = {t:6.2f} s, step {step:5d}  ({k+1}/{nF})")
        return [*lines, title]

    anim = animation.FuncAnimation(fig, update, frames=nF, interval=1000.0/fps, blit=False)

    writer = animation.FFMpegWriter(fps=fps, codec="libx264",
                                    extra_args=["-pix_fmt", "yuv420p"])
    anim.save(out_path, writer=writer, dpi=120)
    print(f"wrote {out_path}  ({nF} frames @ {fps} fps = {nF/fps:.1f} s)")


if __name__ == "__main__":
    main()
