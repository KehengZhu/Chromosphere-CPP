#!/usr/bin/env python3
"""Combined Stage 1 (relaxation) + Stage 2 (heat conduction) movie for model_isentropic.

Concatenates the single-fluid isentropic relaxation run and the heat-conduction run
into one publication-styled MP4 with a stage banner. Nature-figure styling: Arial
sans, clean spines, restrained palette (red = heating hero T, blue = flow V, teal =
density), hero T panel + two subordinate panels, common y-axes across both stages so
the Stage 1 -> Stage 2 contrast is honest.

Usage:
    python util/animate_iso_combined.py            # default inputs/output
    python util/animate_iso_combined.py <stage1.txt> <stage2b.txt> <out.mp4> [fps]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _anim_parallel import save_frames_parallel

# ── Nature-figure mandatory font + SVG rules ────────────────────────────────
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['svg.fonttype'] = 'none'

MP = 1.6726219e-27
KB = 1.380649e-23
G = 273.95
F_ION = 1e-3

# Restrained palette (one signal family + neutral + accent).
C_T = "#B64342"      # heating hero (red)
C_V = "#0F4D92"      # flow (blue)
C_RHO = "#42949E"    # density (teal)
C_P = "#7C6CCF"      # pressure (violet)
C_DOWN = "#E9A6A1"   # downflow fill (pale red, directional cue)
C_AXIS = "#4D4D4D"
C_GRID = "#CFCECE"
BANNER_BG = {1: "#DCE6F2", 2: "#F6CFCB"}   # cool (relax) / warm (heating)
BANNER_FG = {1: "#234E78", 2: "#8A2E2C"}


def load(fn):
    f = open(fn)
    ns, neq = map(int, f.readline().split())
    h = np.array(list(map(float, f.readline().split())))
    frames = []
    cur = None
    t = None

    def flush():
        if cur is None:
            return
        rows = [r for r in cur if len(r) == neq]
        if len(rows) == ns:
            frames.append((t, np.array(rows)))

    for line in f:
        if line.startswith('#'):
            flush()
            t = float(line.split('t =')[1].split('step')[0])
            cur = []
        else:
            try:
                cur.append(list(map(float, line.split())))
            except ValueError:
                cur.append([])
    flush()
    return h, frames


def primitives(U, phi):
    rho_i, rho_n = U[:, 0], U[:, 1]
    V = U[:, 2] / rho_i
    Un = U[:, 3] / rho_n
    p_i = 2 / 3 * U[:, 4] - 1 / 3 * rho_i * V * V - 2 / 3 * rho_i * phi
    p_n = 2 / 3 * U[:, 5] - 1 / 3 * rho_n * Un * Un - 2 / 3 * rho_n * phi
    rho = rho_i + rho_n
    T = (p_i + p_n) / (rho * ((1 + F_ION) * KB / MP))
    return V, T, rho, (p_i + p_n)


def pick(frames, t_targets):
    """Nearest-frame lookup for a list of target times."""
    times = np.array([f[0] for f in frames])
    out = []
    for tt in t_targets:
        k = int(np.argmin(np.abs(times - tt)))
        out.append(frames[k])
    return out


def _style_ax(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_color(C_AXIS)
    ax.spines['bottom'].set_color(C_AXIS)
    ax.tick_params(colors=C_AXIS, labelsize=11)
    ax.grid(True, color=C_GRID, lw=0.6, alpha=0.6)
    ax.set_axisbelow(True)


def _render_frame(k, tmpdir, h, V, T, rho, p, t, stage, nF, ylims):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

    fig = plt.figure(figsize=(12, 7.4), dpi=140)
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1],
                          left=0.075, right=0.975, top=0.785, bottom=0.085,
                          hspace=0.33, wspace=0.235)

    # ── persistent header: states single fluid ──────────────────────────────
    fig.text(0.5, 0.962, "Single-fluid isentropic chromosphere",
             ha='center', va='center', fontsize=19, fontweight='bold',
             color="#272727")
    fig.text(0.5, 0.917,
             "single fluid (ion+neutral slaved)  ·  isentropic initial condition  ·  imposed transition-region boundary",
             ha='center', va='center', fontsize=10.5, color=C_AXIS)

    # ── stage banner (changes between stages) ────────────────────────────────
    if stage == 1:
        label = "Stage 1 — Relaxation"
        sub = "heat conduction off"
    else:
        label = "Stage 2 — Heat conduction enabled"
        sub = "conductive flux q(T) imposed at the top boundary"
    fig.text(0.5, 0.862, f"{label}      t = {t:.0f} s",
             ha='center', va='center', fontsize=15, fontweight='bold',
             color=BANNER_FG[stage],
             bbox=dict(boxstyle="round,pad=0.5", facecolor=BANNER_BG[stage],
                       edgecolor=BANNER_FG[stage], linewidth=1.2))
    fig.text(0.5, 0.818, sub, ha='center', va='center', fontsize=9.5,
             style='italic', color=BANNER_FG[stage])

    # ── Temperature (top-left) ───────────────────────────────────────────────
    axT = fig.add_subplot(gs[0, 0])
    axT.plot(h, T / 1e3, color=C_T, lw=2.6)
    axT.set_ylim(*ylims["T"])
    axT.set_xlim(h.min(), h.max())
    axT.set_ylabel("temperature  [kK]", fontsize=12.5, color=C_AXIS)
    _style_ax(axT)
    axT.text(-0.135, 1.03, "a", transform=axT.transAxes, fontsize=16,
             fontweight='bold', color="#272727")
    # show where heat enters during Stage 2
    if stage == 2:
        axT.annotate("conductive\nflux  q(T)", xy=(h.max(), ylims["T"][1] * 0.80),
                     xytext=(h.max() - 0.34 * (h.max() - h.min()), ylims["T"][1] * 0.86),
                     fontsize=10.5, color=BANNER_FG[2], ha='center', va='center',
                     arrowprops=dict(arrowstyle="-|>", color=BANNER_FG[2], lw=1.8))

    # ── Velocity (top-right) ─────────────────────────────────────────────────
    axV = fig.add_subplot(gs[0, 1])
    axV.axhline(0, color=C_AXIS, lw=0.8)
    axV.fill_between(h, V / 1e3, 0, where=(V < 0), color=C_DOWN, alpha=0.55,
                     linewidth=0)
    axV.plot(h, V / 1e3, color=C_V, lw=2.0)
    axV.set_ylim(*ylims["V"])
    axV.set_xlim(h.min(), h.max())
    axV.set_ylabel("velocity  [km s$^{-1}$]", fontsize=12, color=C_AXIS)
    _style_ax(axV)
    # minor gridlines + a tick every 0.1 km s^-1 on y (so +0.1 is marked) and on x.
    axV.yaxis.set_minor_locator(MultipleLocator(0.1))
    axV.xaxis.set_minor_locator(AutoMinorLocator(2))
    axV.grid(True, which='minor', color=C_GRID, lw=0.45, alpha=0.32)
    axV.tick_params(which='minor', length=3, colors=C_AXIS)
    axV.text(-0.135, 1.03, "b", transform=axV.transAxes, fontsize=16,
             fontweight='bold', color="#272727")
    # direction cues (static, low-key)
    axV.text(0.97, 0.08, "downflow (condensation)", transform=axV.transAxes,
             ha='right', va='bottom', fontsize=8.5, color="#8A2E2C")
    axV.text(0.97, 0.90, "upflow (evaporation)", transform=axV.transAxes,
             ha='right', va='top', fontsize=8.5, color=C_V)

    # ── Density (bottom-left) ────────────────────────────────────────────────
    axR = fig.add_subplot(gs[1, 0])
    axR.semilogy(h, rho, color=C_RHO, lw=2.0)
    axR.set_ylim(*ylims["rho"])
    axR.set_xlim(h.min(), h.max())
    axR.set_xlabel("height  [km]", fontsize=12.5, color=C_AXIS)
    axR.set_ylabel(r"density  [kg m$^{-3}$]", fontsize=12, color=C_AXIS)
    _style_ax(axR)
    axR.text(-0.135, 1.03, "c", transform=axR.transAxes, fontsize=16,
             fontweight='bold', color="#272727")

    # ── Pressure (bottom-right) ──────────────────────────────────────────────
    axP = fig.add_subplot(gs[1, 1])
    axP.semilogy(h, p, color=C_P, lw=2.0)
    axP.set_ylim(*ylims["p"])
    axP.set_xlim(h.min(), h.max())
    axP.set_xlabel("height  [km]", fontsize=12.5, color=C_AXIS)
    axP.set_ylabel("pressure  [Pa]", fontsize=12, color=C_AXIS)
    _style_ax(axP)
    axP.text(-0.135, 1.03, "d", transform=axP.transAxes, fontsize=16,
             fontweight='bold', color="#272727")

    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"))
    plt.close(fig)


def main():
    s1 = sys.argv[1] if len(sys.argv) > 1 else "outputs/_archive/iso_stage1.txt"
    s2 = sys.argv[2] if len(sys.argv) > 2 else "outputs/_archive/iso_stage2b.txt"
    out_path = sys.argv[3] if len(sys.argv) > 3 else "visualization/_archive/iso_combined_stage1_stage2.mp4"
    fps = int(sys.argv[4]) if len(sys.argv) > 4 else 24
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    h1, f1 = load(s1)
    h2, f2 = load(s2)
    phi = G * (h1 - h1[0]) * 1000.0  # same grid for both (ns=400)

    # Time sampling: cluster early (power>1) to capture the transient / TR formation.
    N1, T1_disp = 70, 450.0     # Stage 1: relaxation then steady hold
    N2, T2_disp = 150, 350.0    # Stage 2: TR formation + downflow front
    t1 = (np.linspace(0, 1, N1) ** 1.4) * T1_disp
    t2 = (np.linspace(0, 1, N2) ** 1.5) * T2_disp
    fr1 = pick(f1, t1)
    fr2 = pick(f2, t2)

    # Assemble frame records: (stage, t, U). Add short holds at the seam + end.
    records = [(1, t, U) for (t, U) in fr1]
    records += [records[-1]] * 12                 # hold relaxed Stage 1
    records += [(2, t, U) for (t, U) in fr2]
    records += [records[-1]] * 14                 # hold final Stage 2

    prims = [primitives(U, phi) for (_, _, U) in records]  # (V, T, rho, p)

    def lim_lin(vals, pad=0.06):
        a = np.concatenate(vals)
        a = a[np.isfinite(a)]
        lo, hi = float(a.min()), float(a.max())
        d = pad * (hi - lo)
        return lo - d, hi + d

    def lim_log(vals):
        a = np.concatenate(vals)
        a = a[np.isfinite(a) & (a > 0)]
        return float(a.min()) / 1.4, float(a.max()) * 1.4

    ylims = {
        "T": lim_lin([pr[1] / 1e3 for pr in prims]),
        "V": lim_lin([pr[0] / 1e3 for pr in prims], pad=0.08),
        "rho": lim_log([pr[2] for pr in prims]),
        "p": lim_log([pr[3] for pr in prims]),
    }
    nF = len(records)
    print(f"Stage 1: {len(fr1)} frames (t<= {fr1[-1][0]:.0f}s); "
          f"Stage 2: {len(fr2)} frames (t<= {fr2[-1][0]:.0f}s); total {nF} @ {fps} fps")
    print(f"y-lims  T[kK]={ylims['T']}  V[km/s]={ylims['V']}  rho={ylims['rho']}  p={ylims['p']}")

    save_frames_parallel(
        _render_frame,
        [(k, h1, prims[k][0], prims[k][1], prims[k][2], prims[k][3],
          records[k][1], records[k][0], nF, ylims) for k in range(nF)],
        out_path, fps,
    )


if __name__ == "__main__":
    main()
