#!/usr/bin/env python3
"""Evolution movie for the corona-as-boundary gentle-evaporation run (model_gentle;
docs/gentle_evaporation_plan.md, option A). A C7 column extended into a *resolved*
corona (h ≈ 1.0 → ~11.5 Mm), with the coronal conductive flux q(T) imposed at the
top and ramped to drive gentle conduction-driven evaporation.

Seven profile panels vs height — temperature, total gas pressure, charged-fluid (ion)
velocity V and bulk (mass-weighted) velocity (n_i V + n_n U)/(n_i + n_n) side by side
(both upflow/downflow shaded), n_e/n_n, ionization fraction, field-aligned conductive
flux q∥ — on a SPLIT x-axis that zooms the sub-Mm chromosphere+TR (where the TR burns down)
and compresses the extended corona (where the evaporated mass accumulates). A live
phase chip (RELAXING / q-RAMP / EVAPORATION) states the driver state; a bottom
timeline marks the slow-motion window.

Styling follows the project's Nature-figure convention: restrained palette (one
neutral, one signal-blue, one accent-red family), thin despined axes, sans-serif
type, direct labels over legends.

Usage:
  python util/animate_gentle_column.py <run.txt> <out.mp4> [t_on] [ramp] [enhance]
Env:
  SLOW_T0, SLOW_T1   slow-motion window [s] (default [t_on, t_on+1000])
  SLOW_FACTOR        relative slowdown inside the window (default 4)
  ANIM_MAX_FRAMES    output frame budget (default 500); ANIM_FPS (default 25)
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _anim_parallel import save_frames_parallel

M_I = 1.6726219e-27
K_B = 1.380649e-23
KAPPA0 = 1.0e-11  # SI Spitzer parallel conductivity: kappa = KAPPA0 * T^2.5 [W m^-1 K^-1]
CNI, CNN, CNV, CNU, CEI, CEN = 0, 1, 2, 3, 4, 5
MAX_FRAMES = int(os.environ.get("ANIM_MAX_FRAMES", "500"))

# --- restrained palette: one neutral (ink/grid), one signal (blue), one accent (red) ---
INK       = "#1b1b1b"   # text, spines, profile lines where hue is not informative
PRIMARY   = "#1f4e79"   # signal family: main profiles (T, P, n_e)
PRIMARY_L = "#86a8cc"   # lighter signal: secondary line (n_n)
ACCENT    = "#b4341f"   # accent family: the driver (q∥) and upflow direction
GRID      = "#c9ced6"
ZOOM_END_MM, CHROMO_FRAC = 2.6, 0.5   # split: chromo+TR gets CHROMO_FRAC of the width
REGION_TINT = {"chromo": "#f1efea", "TR": "#fbe4cf", "corona": "#e9eff6"}
PHASE_COL   = {"RELAXING": "#5b6675", "q(T) RAMP": "#cc7a1f", "EVAPORATION": PRIMARY}


def _apply_style():
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "axes.edgecolor": INK,
        "axes.linewidth": 0.8,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "xtick.color": INK, "ytick.color": INK,
        "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
        "legend.frameon": False,
    })


_apply_style()


def read_frames(path):
    with open(path) as f:
        raw = [ln for ln in f.read().splitlines() if ln.strip()]
    ns, neq = (int(x) for x in raw[0].split())
    h_km = np.fromstring(raw[1], sep=" ")
    frames, i = [], 2
    while i < len(raw):
        toks = raw[i].split(); t = float(toks[3]); i += 1
        if i + ns > len(raw):
            break
        d = np.array([np.fromstring(raw[i + j], sep=" ") for j in range(ns)])
        i += ns
        frames.append((t, d))
    return h_km, frames


def primitives(xn):
    ni = xn[:, CNI] / M_I
    nn = xn[:, CNN] / M_I
    v = xn[:, CNV] / xn[:, CNI]
    vn = xn[:, CNU] / np.maximum(xn[:, CNN], 1e-30)
    p_i = 2.0 / 3.0 * xn[:, CEI] - 1.0 / 3.0 * M_I * ni * v * v    # charged (ion+e) pressure
    p_n = 2.0 / 3.0 * xn[:, CEN] - 1.0 / 3.0 * M_I * nn * vn * vn  # neutral pressure
    Ti = p_i / (2.0 * ni * K_B)
    v_tot = (ni * v + nn * vn) / np.maximum(ni + nn, 1e-30)   # bulk (mass-weighted) velocity
    return ni, nn, v, Ti, p_i + p_n, v_tot


def phase_of(t, t_on, ramp, enhance):
    """Return (label, q-enhancement amp) for the q(T)-ramp timeline."""
    if enhance == 1.0 or t < t_on:
        return "RELAXING", 1.0
    if t <= t_on + ramp:
        w = 0.5 * (1.0 - np.cos(np.pi * (t - t_on) / ramp))
        return "q(T) RAMP", 1.0 + (enhance - 1.0) * w
    return "EVAPORATION", enhance


def region_seps(h_Mm, Ti0):
    """chromo|TR (T=2e4) and TR|corona (T=5e5) boundaries from the initial profile."""
    def cross(level):
        j = int(np.argmax(Ti0 >= level))
        return float(h_Mm[j]) if Ti0[j] >= level else float(h_Mm[-1])
    return cross(2.0e4), cross(5.0e5)


def _render(k, tmpdir, h_km, prim, t, t_on, ramp, enhance,
            ylims, kd, kp, xticks, seps, t_max, slow):
    import matplotlib.pyplot as plt
    _apply_style()
    ni, nn, v, Ti, p_tot, v_tot = prim
    hMm = h_km / 1e3
    q_cond = -KAPPA0 * np.power(np.clip(Ti, 1.0, None), 2.5) * np.gradient(Ti, h_km * 1e3)
    fwd = lambda x: np.interp(x, kd, kp)        # piecewise-linear split transform
    inv = lambda x: np.interp(x, kp, kd)
    ct, cb = seps

    fig = plt.figure(figsize=(13.0, 14.4))
    gs = fig.add_gridspec(4, 2, left=0.075, right=0.985, top=0.912, bottom=0.072,
                          hspace=0.46, wspace=0.205)
    a_T  = fig.add_subplot(gs[0, 0]); a_P  = fig.add_subplot(gs[0, 1])
    a_V  = fig.add_subplot(gs[1, 0]); a_Vt = fig.add_subplot(gs[1, 1])
    a_n  = fig.add_subplot(gs[2, 0]); a_f  = fig.add_subplot(gs[2, 1])
    a_q  = fig.add_subplot(gs[3, 0])
    fig.add_subplot(gs[3, 1]).axis("off")   # 8th slot unused (7 panels)

    # ---- (a) temperature ----
    a_T.semilogy(hMm, np.clip(Ti, 3e3, None), color=PRIMARY, lw=1.6)
    a_T.set_ylabel(r"$T_i$  [K]"); a_T.set_ylim(*ylims["T"])
    a_T.set_title("temperature — TR burns down, corona heats")

    # ---- (b) total gas pressure ----
    a_P.semilogy(hMm, np.clip(p_tot, 1e-4, None), color=PRIMARY, lw=1.6)
    a_P.set_ylabel(r"$P$  [Pa]"); a_P.set_ylim(*ylims["P"])
    a_P.set_title("total gas pressure")

    # ---- (c) velocity: upflow (accent) vs downflow (signal) ----
    vk = v / 1e3
    a_V.fill_between(hMm, 0, vk, where=vk > 0, interpolate=True, color=ACCENT, alpha=0.30, lw=0)
    a_V.fill_between(hMm, 0, vk, where=vk < 0, interpolate=True, color=PRIMARY, alpha=0.18, lw=0)
    a_V.plot(hMm, vk, color=INK, lw=1.3)
    a_V.axhline(0, color=INK, lw=0.7)
    a_V.set_ylabel(r"$V$  [km s$^{-1}$]"); a_V.set_ylim(*ylims["V"])
    a_V.set_title(r"charged-fluid (ion) velocity $V$")
    a_V.text(0.985, 0.93, "upflow", color=ACCENT, fontsize=9, ha="right", va="top",
             transform=a_V.transAxes)
    a_V.text(0.985, 0.07, "downflow", color=PRIMARY, fontsize=9, ha="right", va="bottom",
             transform=a_V.transAxes)

    # ---- (d) bulk (mass-weighted) velocity  (n_i V + n_n U)/(n_i + n_n) ----
    vtk = v_tot / 1e3
    a_Vt.fill_between(hMm, 0, vtk, where=vtk > 0, interpolate=True, color=ACCENT, alpha=0.30, lw=0)
    a_Vt.fill_between(hMm, 0, vtk, where=vtk < 0, interpolate=True, color=PRIMARY, alpha=0.18, lw=0)
    a_Vt.plot(hMm, vtk, color=INK, lw=1.3)
    a_Vt.axhline(0, color=INK, lw=0.7)
    a_Vt.set_ylabel(r"$V_\mathrm{tot}$  [km s$^{-1}$]"); a_Vt.set_ylim(*ylims["V"])
    a_Vt.set_title(r"bulk velocity  $(n_iV+n_nU)/(n_i+n_n)$")

    # ---- (e) density: direct-labelled, no legend ----
    a_n.semilogy(hMm, np.clip(ni, 1e11, None), color=PRIMARY, lw=1.6)
    a_n.semilogy(hMm, np.clip(nn, 1e11, None), color=PRIMARY_L, lw=1.4, ls="--")
    a_n.set_ylabel(r"density  [m$^{-3}$]"); a_n.set_ylim(*ylims["n"])
    a_n.set_title("density — corona fills")
    a_n.text(0.985, 0.93, r"$n_e$", color=PRIMARY, fontsize=11, fontweight="bold",
             ha="right", va="top", transform=a_n.transAxes)
    a_n.text(0.985, 0.80, r"$n_n$", color=PRIMARY_L, fontsize=11, fontweight="bold",
             ha="right", va="top", transform=a_n.transAxes)

    # ---- (f) ionization fraction ----
    a_f.plot(hMm, ni / np.maximum(ni + nn, 1.0), color=PRIMARY, lw=1.6)
    a_f.set_ylabel("ionization fraction"); a_f.set_ylim(-0.03, 1.05)
    a_f.set_xlabel("height  [Mm]"); a_f.set_title("ionization")

    # ---- (g) conductive flux — the driver ----
    a_q.fill_between(hMm, 0, q_cond, color=ACCENT, alpha=0.16, lw=0)
    a_q.plot(hMm, q_cond, color=ACCENT, lw=1.6)
    a_q.axhline(0, color=INK, lw=0.7)
    a_q.set_ylabel(r"$q_\parallel$  [W m$^{-2}$]   (+ up)"); a_q.set_ylim(*ylims["q"])
    a_q.set_xlabel("height  [Mm]"); a_q.set_title("conductive flux — the driver (ramped at top)")

    # ---- shared split x-axis, region bands, panel tags ----
    tags = "abcdefg"
    for ax, tag in zip((a_T, a_P, a_V, a_Vt, a_n, a_f, a_q), tags):
        ax.set_xscale("function", functions=(fwd, inv))
        ax.set_xlim(kd[0], kd[-1])
        ax.set_xticks([tk for tk in xticks if kd[0] <= tk <= kd[-1]])
        for x0, x1, rt in ((kd[0], ct, "chromo"), (ct, cb, "TR"), (cb, kd[-1], "corona")):
            ax.axvspan(x0, x1, color=REGION_TINT[rt], alpha=0.85, lw=0, zorder=0)
        for sp in (ct, cb):
            ax.axvline(sp, color="dimgray", ls="--", lw=0.8, alpha=0.55, zorder=1)
        ax.text(-0.085, 1.04, tag, transform=ax.transAxes, fontsize=14,
                fontweight="bold", va="top", ha="left", color=INK)

    # ---- header: title, phase chip, time readout ----
    phase, amp = phase_of(t, t_on, ramp, enhance)
    fig.text(0.075, 0.965, "Gentle conduction-driven evaporation — corona resolved as a volume",
             fontsize=15, fontweight="bold", color=INK, va="center")
    fig.text(0.075, 0.945, f"  {phase}    q(T) ×{amp:.2f}  ", fontsize=12, color="white",
             fontweight="bold", va="center", ha="left",
             bbox=dict(boxstyle="round,pad=0.32", fc=PHASE_COL[phase], ec="none"))
    fig.text(0.985, 0.945, f"t = {t:6.0f} s", fontsize=13, color=INK, va="center", ha="right")

    # ---- bottom timeline with the slow-motion window ----
    a_time = fig.add_axes([0.075, 0.026, 0.91, 0.012])
    a_time.set_xlim(0, t_max); a_time.set_ylim(0, 1)
    s0, s1, sf = slow
    a_time.axvspan(s0, s1, color=ACCENT, alpha=0.16, lw=0)
    a_time.hlines(0.5, 0, t_max, color=GRID, lw=2.2)
    a_time.hlines(0.5, 0, t, color=PRIMARY, lw=2.2)
    a_time.plot([t], [0.5], "o", color=PRIMARY, ms=6)
    a_time.text((s0 + s1) / 2, 1.15, f"slow-motion ×{sf:.0f}", color=ACCENT, fontsize=9,
                ha="center", va="bottom")
    a_time.set_yticks([]); a_time.set_xlabel("simulation time  [s]", fontsize=9.5, labelpad=2)
    for sp in ("top", "right", "left"):
        a_time.spines[sp].set_visible(False)
    a_time.tick_params(labelsize=8.5, length=3)

    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=130)
    plt.close(fig)


def main():
    run = sys.argv[1] if len(sys.argv) > 1 else "outputs/gentle_v3_evap.txt"
    out = sys.argv[2] if len(sys.argv) > 2 else "visualization/gentle_v3_evolution.mp4"
    t_on    = float(sys.argv[3]) if len(sys.argv) > 3 else 3000.0
    ramp    = float(sys.argv[4]) if len(sys.argv) > 4 else 300.0
    enhance = float(sys.argv[5]) if len(sys.argv) > 5 else 3.0
    os.makedirs(os.path.dirname(out), exist_ok=True)

    h_km, frames = read_frames(run)
    times = np.array([fr[0] for fr in frames])
    t_max = float(times[-1])

    # --- slow-motion time-remap: weight frames inside [SLOW_T0, SLOW_T1] by SLOW_FACTOR,
    #     then sample MAX_FRAMES uniformly in cumulative-weight space → more output frames
    #     (slower playback) there. Each output frame keeps its true sim-time, so phase
    #     labels/annotations stay correct. ---
    s0 = float(os.environ.get("SLOW_T0", t_on))
    s1 = float(os.environ.get("SLOW_T1", t_on + 1000.0))
    sf = float(os.environ.get("SLOW_FACTOR", 4.0))
    w = np.where((times >= s0) & (times <= s1), sf, 1.0)
    cumw = np.cumsum(w); cumw = cumw / cumw[-1]
    sel = np.clip(np.searchsorted(cumw, np.linspace(0, 1, MAX_FRAMES)), 0, len(frames) - 1)

    prims = [primitives(frames[i][1]) for i in sel]
    sel_t = [float(times[i]) for i in sel]
    print(f"{run}: {len(h_km)} cells, {len(frames)}→{len(sel)} frames "
          f"(slow ×{sf:.0f} in [{s0:.0f},{s1:.0f}] s), t<= {t_max:.0f} s")

    # --- axis limits over the shown frames ---
    def lim(idx, log, lo_floor=0.0):
        a = np.concatenate([p[idx] for p in prims])
        if log:
            a = a[a > lo_floor]
            return a.min() / 1.3, a.max() * 1.3
        lo, hi = a.min(), a.max(); pad = 0.08 * (hi - lo if hi > lo else 1.0)
        return lo - pad, hi + pad
    ylims = {"T": lim(3, True, 3e3), "P": lim(4, True, 1e-4), "n": lim(0, True, 1e11)}
    nn_lo = np.concatenate([p[1] for p in prims]); nn_lo = nn_lo[nn_lo > 1e11]
    ylims["n"] = (min(ylims["n"][0], nn_lo.min() / 1.3), ylims["n"][1])
    vmax = max(abs(np.concatenate([np.concatenate([p[2], p[5]]) for p in prims])).max() / 1e3, 5.0)
    ylims["V"] = (-1.12 * vmax, 1.12 * vmax)
    qall = np.concatenate([
        -KAPPA0 * np.power(np.clip(p[3], 1.0, None), 2.5) * np.gradient(p[3], h_km * 1e3)
        for p in prims])
    qmax = max(abs(qall).max(), 1.0)
    ylims["q"] = (-1.15 * qmax, 1.15 * qmax)

    # --- split x-axis transform + region boundaries from the first frame ---
    h0, htop = h_km[0] / 1e3, h_km[-1] / 1e3
    kd = np.array([h0, ZOOM_END_MM, htop])
    kp = np.array([0.0, CHROMO_FRAC, 1.0])
    xticks = np.array([1.0, 1.5, 2.0, 2.5, 4.0, 6.0, 8.0, 10.0])
    seps = region_seps(h_km / 1e3, primitives(frames[0][1])[3])

    fps = int(os.environ.get("ANIM_FPS", "25"))
    save_frames_parallel(
        _render,
        [(k, h_km, prims[k], sel_t[k], t_on, ramp, enhance,
          ylims, kd, kp, xticks, seps, t_max, (s0, s1, sf)) for k in range(len(sel))],
        out, fps)


if __name__ == "__main__":
    main()
