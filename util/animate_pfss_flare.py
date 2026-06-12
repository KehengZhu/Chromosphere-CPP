#!/usr/bin/env python3
"""Evolution movie for a flaring PFSS field line (PFSS_FLARE run).

Reads a chromo_main output (default outputs/pfss_lineA_flare.txt) and renders a
4-panel vertical-profile animation vs height — T_i, V, densities (n_i, n_n), and
ionization fraction — with the beam deposition window shaded and a beam-status
clock. Frames are rendered in parallel and stitched to MP4 (util/_anim_parallel).

Temperature uses the full cons2prim (state.cpp): p_i = 2/3 E_i - 1/3 rho_i V^2
- 2/3 rho_i phi_g, phi_g = g (h - h_base) from the cell-center height — so T is
accurate (plot_output.py / visualize_flare_frames.py drop the phi_g term).

Usage:
  python util/animate_pfss_flare.py [output.txt] [out.mp4]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _anim_parallel import save_frames_parallel

m_i = 1.6726219e-27
m_n = m_i
k_b = 1.380649e-23
g_si = 274.0
H_BASE_KM = 1003.0

# Beam window (FLARE_* defaults used by the PFSS_FLARE overlay).
H_LO, H_HI = 1600.0, 2000.0
T_ON, T_OFF = 2.0, 12.0

# Fixed axis/value ranges so the animation does not jump frame-to-frame.
T_MIN, T_MAX = 3.0e3, 6.0e7
V_MIN, V_MAX = -650.0, 300.0       # km/s
N_MIN, N_MAX = 1.0e13, 1.0e20      # m^-3
MAX_FRAMES = 500                    # subsample longer runs to keep movie ~length


def read_frames(path):
    with open(path) as f:
        raw = [ln for ln in f.read().splitlines() if ln.strip()]
    ns, neq = (int(x) for x in raw[0].split())
    heights = np.fromstring(raw[1], sep=" ")
    frames, i = [], 2
    while i < len(raw):
        toks = raw[i].split()
        t = float(toks[3])
        i += 1
        data = np.zeros((ns, neq))
        ok = True
        for j in range(ns):
            tk = raw[i + j].split()
            if len(tk) != neq:
                ok = False
                break
            data[j] = [float(x) for x in tk]
        i += ns
        if ok:
            frames.append((t, data))
    return heights, frames


def cell_center_km(heights):
    prev = np.concatenate(([H_BASE_KM], heights[:-1]))
    return 0.5 * (heights + prev)


def primitives(xn, phi_g):
    ni = xn[:, 0] / m_i
    nn = xn[:, 1] / m_n
    v = xn[:, 2] / xn[:, 0]
    p_i = 2.0 / 3.0 * xn[:, 4] - 1.0 / 3.0 * xn[:, 0] * v * v - 2.0 / 3.0 * xn[:, 0] * phi_g
    Ti = p_i / (2.0 * ni * k_b)
    return ni, nn, v, Ti


def _render_frame(k, tmpdir, H, t, Ti, V, ni, nn, tmax):
    """Render one frame to tmpdir/frame_{k:06d}.png. Top-level for pickling."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    phase = "PRE-BEAM" if t < T_ON else ("BEAM ON" if t <= T_OFF else "POST-BEAM")
    pcol = {"PRE-BEAM": "0.4", "BEAM ON": "orangered", "POST-BEAM": "navy"}[phase]
    fig.suptitle(f"Flaring PFSS line (pfss_lineA)   t = {t:6.2f} s    [{phase}]",
                 fontsize=13, color=pcol)

    a_T, a_V, a_n, a_f = ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1]

    a_T.semilogy(H, np.clip(Ti, T_MIN, T_MAX), color="crimson", lw=1.6)
    a_T.set_ylabel(r"$T_i$ [K]"); a_T.set_ylim(T_MIN, T_MAX)

    a_V.plot(H, V / 1e3, color="C0", lw=1.6)
    a_V.axhline(0, color="k", lw=0.6)
    a_V.set_ylabel(r"$V$ [km s$^{-1}$]  (+up / $-$down)"); a_V.set_ylim(V_MIN, V_MAX)

    a_n.semilogy(H, np.clip(ni, N_MIN, None), color="C3", lw=1.5, label=r"$n_i$")
    a_n.semilogy(H, np.clip(nn, N_MIN, None), color="C0", lw=1.5, label=r"$n_n$")
    a_n.set_ylabel(r"number density [m$^{-3}$]"); a_n.set_ylim(N_MIN, N_MAX)
    a_n.set_xlabel("height [km]"); a_n.legend(loc="upper right", fontsize=9)

    f_ion = ni / np.maximum(ni + nn, 1.0)
    a_f.plot(H, f_ion, color="purple", lw=1.6)
    a_f.set_ylabel(r"ionization fraction $n_i/(n_i+n_n)$"); a_f.set_ylim(-0.03, 1.05)
    a_f.set_xlabel("height [km]")

    for a in (a_T, a_V, a_n, a_f):
        a.axvspan(H_LO, H_HI, color="orange", alpha=0.12)
        a.grid(True, alpha=0.25)
        a.set_xlim(H[0], H[-1])

    # progress bar in the suptitle region
    fig.text(0.5, 0.945, "", ha="center")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=110)
    plt.close(fig)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    in_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "outputs", "pfss_lineA_flare.txt")
    out_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "visualization", "pfss_lineA_flare_evolution.mp4")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    H, frames = read_frames(in_path)
    phi_g = g_si * (cell_center_km(H) - H_BASE_KM) * 1.0e3

    # Subsample to keep the movie a reasonable length.
    if len(frames) > MAX_FRAMES:
        idx = np.linspace(0, len(frames) - 1, MAX_FRAMES).round().astype(int)
        frames = [frames[i] for i in idx]
    tmax = frames[-1][0]
    print(f"{in_path}: {len(H)} cells, {len(frames)} frames, t<= {tmax:.2f} s")

    frame_args = []
    for k, (t, xn) in enumerate(frames):
        ni, nn, v, Ti = primitives(xn, phi_g)
        frame_args.append((k, H, t, Ti, v, ni, nn, tmax))

    save_frames_parallel(_render_frame, frame_args, out_path, fps=30)


if __name__ == "__main__":
    main()
