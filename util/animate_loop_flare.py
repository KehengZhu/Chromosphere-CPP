#!/usr/bin/env python3
"""Evolution movie for a flaring loop (full chromosphere→corona domain).

Topology-aware (read from the .dat [META] topology= key):

  * full  — the WHOLE loop footpoint A→apex→footpoint B. Profiles are plotted vs
    ARC LENGTH s (NOT height: on a closed loop height folds the two legs onto each
    other). The split x-axis zooms BOTH footpoints (thin sub-Mm chromosphere+TR)
    and compresses the extended coronal middle, laying the loop out flat
    [footpoint A | apex | footpoint B]. The beam lands on footpoint A only.

  * closed/open (half loop) — profiles vs TRUE HEIGHT (from the .dat phi_g), single
    split: zoom the sub-Mm chromosphere+TR, compress the corona.

Region boundaries (chromosphere|TR|corona, both legs) are marked with translucent
dashed vlines; the apex (full loop) with a dotted line. (Split-axis convention noted
in CLAUDE.md → Visualization.) Frames rendered in parallel, stitched to MP4.

Usage:
  python util/animate_loop_flare.py [loop.dat] [flare.txt] [out.mp4]

Examples:
# Closed ensemble loop
FLARE_DELTA=4 FLARE_E_CUT=20 FLARE_T_ON=2 FLARE_DUR=10 python3 util/animate_loop_flare.py \
  scenarios/data/event_ensemble/c05_101Mm.dat outputs/event_ensemble/c05_101Mm.txt \
  util/visualization/event_ensemble_c05_101Mm_closed.mp4

# Open ensemble line
FLARE_DELTA=4 FLARE_E_CUT=20 FLARE_T_ON=2 FLARE_DUR=10 python3 util/animate_loop_flare.py \
  scenarios/data/event_ensemble/o00_1040Mm.dat outputs/event_ensemble/o00_1040Mm.txt \
  util/visualization/event_ensemble_o00_1040Mm_open.mp4

# Event starter (full loop, AR 13768)
FLARE_DELTA=4 FLARE_E_CUT=20 FLARE_T_ON=2 FLARE_DUR=10 python3 util/animate_loop_flare.py \
  scenarios/data/event_20240801_AR13768.dat outputs/event_20240801_AR13768_flare.txt \
  util/visualization/event_20240801_AR13768_flare.mp4

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

# Beam timing window [s] and thick-target parameters. Read from the same FLARE_*
# env vars the scenario uses so the overlaid deposition profile Q(s) matches the
# run it is animating (fall back to the scenario defaults if unset).
# The self-consistent beam is injected at the apex and attenuates down each leg;
# its deposition profile Q(s) is recomputed per frame and overlaid (orange) so you
# can watch it migrate up the loop as the corona fills. (See physics.hpp.)
T_ON       = float(os.environ.get("FLARE_T_ON", 2.0))
T_OFF      = T_ON + float(os.environ.get("FLARE_DUR", 10.0))
BEAM_E_CUT = float(os.environ.get("FLARE_E_CUT", 20.0))   # low-energy cutoff E_c [keV]
BEAM_DELTA = float(os.environ.get("FLARE_DELTA", 5.0))    # spectral index δ

T_MIN, T_MAX = 3.0e3, 5.0e7
V_MIN, V_MAX = -700.0, 1300.0
N_MIN, N_MAX = 1.0e11, 1.0e20
MAX_FRAMES = 500

# Split x-axis: keep this much arc/height near EACH footpoint at full (zoomed) scale,
# then allocate the plot width so the chromosphere(+TR) gets CHROMO_PLOT_FRAC of it and
# the (extended) corona the rest — the corona is compressed by whatever factor that
# implies, independent of loop length. A continuous piecewise-linear transform
# (np.interp over knots) maps data-x → plot-x so both the thin chromosphere+TR and the
# extended corona are legible. CHROMO_PLOT_FRAC=0.6 → chromosphere:corona = 6:4.
ZOOM_MM = 2.0
CHROMO_PLOT_FRAC = 0.90


def read_meta(path):
    meta, in_meta = {}, False
    for ln in open(path):
        s = ln.strip()
        if s.startswith("[META]"):
            in_meta = True; continue
        if s.startswith("[") and not s.startswith("[META]"):
            in_meta = False
        if in_meta and "=" in s and not s.startswith("#"):
            k, v = s.split("=", 1)
            meta[k.strip()] = v.strip()
    return meta


def parse_dat(path):
    """Return per-cell ds[m], arc-length center[m], phi_g[J/kg], height[m]."""
    rows, in_cells = [], False
    for ln in open(path):
        s = ln.strip()
        if s.startswith("[CELLS]"):
            in_cells = True; continue
        if s.startswith("[GHOSTS]"):
            break
        if not in_cells or not s or s.startswith("#"):
            continue
        rows.append([float(x) for x in s.split()])
    a = np.array(rows)
    ds = a[:, 1]
    phi_c = 0.5 * (a[:, 4] + a[:, 5])
    s_arc = np.cumsum(ds) - 0.5 * ds          # cell-center arc length [m]
    return ds, s_arc, phi_c, phi_c / g_si


def read_frames(path):
    raw = [l for l in open(path).read().splitlines() if l.strip()]
    ns, neq = map(int, raw[0].split())
    i, frames = 2, []
    while i < len(raw):
        t = float(raw[i].split()[3]); i += 1
        d = np.zeros((ns, neq)); ok = True
        for j in range(ns):
            tk = raw[i + j].split()
            if len(tk) != neq:
                ok = False; break
            d[j] = [float(x) for x in tk]
        i += ns
        if ok:
            frames.append((t, d))
    return frames


def primitives(xn, phi_g):
    ni = xn[:, 0] / m_i
    nn = xn[:, 1] / m_n
    v = xn[:, 2] / xn[:, 0]
    p_i = 2/3*xn[:, 4] - 1/3*xn[:, 0]*v*v - 2/3*xn[:, 0]*phi_g
    Ti = p_i / (2.0 * ni * k_b)
    return ni, nn, v, Ti


def thick_target_Q(ni, nn, ds, phi_g):
    """Thick-target deposition profile Q(s), peak-normalized (matches physics.hpp):
    inject at the apex (max φ_g), column N(s)=∫n ds from the apex down each leg,
    Q ∝ n_tot·(1+N/N_c)^{-δ/2}, N_c = 2e21·E_c²."""
    ntot = ni + nn
    a = int(np.argmax(phi_g))
    N = np.zeros_like(ntot)
    for i in range(a - 1, -1, -1):
        N[i] = N[i + 1] + 0.5 * (ntot[i + 1] * ds[i + 1] + ntot[i] * ds[i])
    for i in range(a + 1, len(ntot)):
        N[i] = N[i - 1] + 0.5 * (ntot[i] * ds[i] + ntot[i - 1] * ds[i - 1])
    N_c = 2.0e21 * BEAM_E_CUT ** 2
    w = ntot * (1.0 + N / N_c) ** (-0.5 * BEAM_DELTA)
    return w / w.max() if w.max() > 0 else w


def _round_ticks(lo, hi, target=4):
    """A few 'nice' round tick values in (lo, hi]."""
    if hi <= lo:
        return []
    raw = (hi - lo) / target
    mag = 10.0 ** np.floor(np.log10(raw))
    step = min([1, 2, 2.5, 5, 10], key=lambda m: abs(m * mag - raw)) * mag
    start = np.ceil((lo + 1e-9) / step) * step
    return list(np.arange(start, hi + step * 0.5, step))


def _spread_ticks(cand, kd, kp, min_frac=0.06):
    """Keep ticks whose PLOT positions are ≥ min_frac of the plot width apart, so the
    compressed corona never crams its labels. Endpoints are always kept."""
    cand = sorted({float(c) for c in cand if kd[0] - 1e-9 <= c <= kd[-1] + 1e-9})
    if not cand:
        return cand
    span = kp[-1] - kp[0]
    keep, last_p = [], -1e18
    for c in cand:
        p = float(np.interp(c, kd, kp))
        if p - last_p >= min_frac * span:
            keep.append(c); last_p = p
    if cand[-1] not in keep:                                # always label the far end
        if float(np.interp(cand[-1], kd, kp)) - last_p < min_frac * span and len(keep) > 1:
            keep[-1] = cand[-1]
        else:
            keep.append(cand[-1])
    return keep


def build_transform(xmax, topology):
    """Return (knots_data, knots_plot, xticks) for the piecewise-linear x transform.

    Plot width is split by CHROMO_PLOT_FRAC: the zoomed chromosphere(+TR) gets that
    fraction, the corona the remainder (corona:chromo = (1-f)/f), so e.g. f=0.6 lays the
    plot out 60% chromosphere / 40% corona regardless of how long the corona is. Ticks
    are spread in PLOT space so the compressed corona stays legible."""
    z = ZOOM_MM
    corona_over_chromo = (1.0 - CHROMO_PLOT_FRAC) / CHROMO_PLOT_FRAC
    if topology == "full":
        Stot = xmax
        chromo_w = 2.0 * z                                 # two footpoint zooms
        corona_w = chromo_w * corona_over_chromo           # → 6:4 plot split
        mid = z + corona_w                                 # plot-x at the far-leg knot
        kd = np.array([0.0, z, Stot - z, Stot])            # data knots
        kp = np.array([0.0, z, mid, mid + z])              # plot knots
        cand = ([0, z] + _round_ticks(z, Stot - z, 4) + [Stot / 2.0]
                + [Stot - z, Stot])
        return kd, kp, _spread_ticks(cand, kd, kp)
    # half loop / open: zoom [0, z], corona [z, xmax] compressed to the (1-f) fraction
    corona_w = z * corona_over_chromo
    kd = np.array([0.0, z, xmax])
    kp = np.array([0.0, z, z + corona_w])
    cand = [0, 0.5, 1, 1.5, z] + _round_ticks(z, xmax, 4) + [xmax]
    return kd, kp, _spread_ticks(cand, kd, kp)


def region_layout(Ti0, x, topology, Stot):
    """Separator x-positions, region labels, and apex-x from the initial T profile.

    Boundaries: chromosphere|TR where T first exceeds 2e4 K (going up from the
    footpoint), TR|corona where it reaches 5e5 K. On a full loop the far leg mirrors
    leg A about the apex."""
    def first_cross(level):
        j = int(np.argmax(Ti0 >= level))
        return float(x[j]) if Ti0[j] >= level else float(x[-1])

    ct, cb = first_cross(2.0e4), first_cross(5.0e5)
    if topology == "full":
        apex = Stot / 2.0
        seps = [ct, cb, Stot - cb, Stot - ct]
        labels = [(ct / 2, "chr"), ((ct + cb) / 2, "TR"),
                  ((cb + apex) / 2, "corona"),
                  (Stot - (ct + cb) / 2, "TR"), (Stot - ct / 2, "chr")]
        return seps, labels, apex
    seps = [ct, cb]
    labels = [(ct / 2, "chromo"), ((ct + cb) / 2, "TR"), ((cb + x[-1]) / 2, "corona")]
    return seps, labels, None


def _render_frame(k, tmpdir, xMm, t, Ti, V, ni, nn, ds_arr, phi_arr,
                  kd, kp, xticks, seps, labels, apex_x, xlabel, vlabel, title_prefix):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    fwd = lambda x: np.interp(x, kd, kp)        # piecewise-linear, monotonic
    inv = lambda x: np.interp(x, kp, kd)
    beam_on = (T_ON <= t <= T_OFF)
    Qbeam = thick_target_Q(ni, nn, ds_arr, phi_arr) if beam_on else None

    fig, ax = plt.subplots(2, 2, figsize=(12.5, 8), sharex=True)
    phase = "PRE-BEAM" if t < T_ON else ("BEAM ON" if t <= T_OFF else "POST-BEAM")
    pcol = {"PRE-BEAM": "0.4", "BEAM ON": "orangered", "POST-BEAM": "navy"}[phase]
    fig.suptitle(f"{title_prefix}   t = {t:6.2f} s    [{phase}]", fontsize=13, color=pcol)
    a_T, a_V, a_n, a_f = ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1]

    a_T.semilogy(xMm, np.clip(Ti, T_MIN, T_MAX), color="crimson", lw=1.4)
    a_T.set_ylabel(r"$T_i$ [K]"); a_T.set_ylim(T_MIN, T_MAX)

    a_V.plot(xMm, V / 1e3, color="C0", lw=1.4); a_V.axhline(0, color="k", lw=0.6)
    a_V.set_ylabel(vlabel); a_V.set_ylim(V_MIN, V_MAX)

    a_n.semilogy(xMm, np.clip(ni, N_MIN, None), color="C3", lw=1.3, label=r"$n_e$")
    a_n.semilogy(xMm, np.clip(nn, N_MIN, None), color="C0", lw=1.3, label=r"$n_n$")
    a_n.set_ylabel(r"density [m$^{-3}$]"); a_n.set_ylim(N_MIN, N_MAX)
    a_n.set_xlabel(xlabel); a_n.legend(loc="upper right", fontsize=9)

    a_f.plot(xMm, ni / np.maximum(ni + nn, 1.0), color="purple", lw=1.4)
    a_f.set_ylabel("ionization fraction"); a_f.set_ylim(-0.03, 1.05)
    a_f.set_xlabel(xlabel)

    # Region bands (chromosphere | TR | corona, mirrored on a full loop) — light
    # background tints make the regions legible even where the TR is razor-thin, so
    # the boundaries don't read as one smudged pair of dashed lines.
    edges = [kd[0]] + list(seps) + [kd[-1]]
    rtypes = (["chromo", "TR", "corona", "TR", "chromo"] if len(seps) == 4
              else ["chromo", "TR", "corona"])
    rcol = {"chromo": "#cfe0f2", "TR": "#d6efd6", "corona": "#fae1c6"}

    for a in (a_T, a_V, a_n, a_f):
        a.set_xscale("function", functions=(fwd, inv))    # zoom footpoints, compress corona
        a.set_xlim(kd[0], kd[-1])
        a.set_xticks([tk for tk in xticks if kd[0] <= tk <= kd[-1]])
        for x0, x1, rt in zip(edges[:-1], edges[1:], rtypes):
            a.axvspan(x0, x1, color=rcol[rt], alpha=0.55, lw=0, zorder=0)
        for sp in seps:
            a.axvline(sp, color="dimgray", ls="--", lw=0.9, alpha=0.6)     # region boundaries
        if apex_x is not None:
            a.axvline(apex_x, color="seagreen", ls=":", lw=1.3, alpha=0.65)  # apex = beam injection
        a.grid(True, alpha=0.18)

    # Live thick-target beam deposition Q(s) (peak-normalized) shaded on the T panel,
    # only while the beam is on — it injects at the apex and migrates as the loop fills.
    if Qbeam is not None:
        aQ = a_T.twinx()
        aQ.fill_between(xMm, 0.0, Qbeam, color="orange", alpha=0.30, lw=0)
        aQ.set_ylim(0.0, 4.0); aQ.set_yticks([])     # fill occupies bottom ~1/4
        aQ.set_zorder(a_T.get_zorder() - 1); a_T.patch.set_visible(False)

    trans = a_T.get_xaxis_transform()
    for xpos, lab in labels:
        a_T.text(xpos, 1.02, lab, transform=trans, ha="center", va="bottom",
                 fontsize=8, color="dimgray", alpha=0.85)
    if apex_x is not None:
        a_T.text(apex_x, 0.98, "apex", transform=trans, ha="center", va="top",
                 fontsize=8, color="seagreen", alpha=0.9)
    if beam_on:
        a_T.text(0.015, 0.93, "beam: thick-target (orange = deposition Q)",
                 transform=a_T.transAxes, fontsize=7.5, color="darkorange")

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=110)
    plt.close(fig)


def main():
    dat = sys.argv[1] if len(sys.argv) > 1 else "scenarios/data/loop_full_fine.dat"
    flare = sys.argv[2] if len(sys.argv) > 2 else "outputs/loop_full_flare.txt"
    out = sys.argv[3] if len(sys.argv) > 3 else "util/visualization/loop_full_flare_evolution.mp4"
    os.makedirs(os.path.dirname(out), exist_ok=True)

    meta = read_meta(dat)
    topology = meta.get("topology", "closed")
    ds, s_arc, phi_g, height = parse_dat(dat)
    Stot = float(ds.sum())                       # total arc length [m]

    if topology == "full":
        xMm = s_arc / 1e6
        xmax = Stot / 1e6
        xlabel = "arc length s [Mm]   (footpoint A | apex | footpoint B)"
        vlabel = r"$V$ [km/s]  (+ along $s$: A$\to$apex$\to$B)"
        title_prefix = f"Full closed loop ({xmax:.0f} Mm): apex-injected thick-target beam"
    else:
        xMm = height / 1e6
        xmax = float(height.max() / 1e6)
        xlabel = "height above footpoint [Mm]  (zoom <2 Mm | corona compressed)"
        vlabel = r"$V$ [km/s] (+up the loop)"
        title_prefix = f"Half loop ({xmax:.0f} Mm): flare"

    frames = read_frames(flare)
    # Optional 4th arg: cap the movie at t_max seconds (the run may extend well past
    # the window of interest, e.g. a 100-min run when only the 25-min cooldown matters).
    if len(sys.argv) > 4:
        t_cap = float(sys.argv[4])
        frames = [f for f in frames if f[0] <= t_cap]
    if len(frames) > MAX_FRAMES:
        idx = np.linspace(0, len(frames) - 1, MAX_FRAMES).round().astype(int)
        frames = [frames[i] for i in idx]
    print(f"{flare}: topology={topology}, {len(xMm)} cells, {len(frames)} frames, "
          f"t<= {frames[-1][0]:.1f} s")

    kd, kp, xticks = build_transform(xmax, topology)
    _ni0, _nn0, _v0, Ti0 = primitives(frames[0][1], phi_g)
    seps, labels, apex_x = region_layout(Ti0, xMm, topology, xmax)
    sep_str = ", ".join(f"{s:.2f}" for s in seps)
    print(f"  region separators [Mm]: {sep_str}" + (f"; apex @ {apex_x:.1f}" if apex_x else ""))

    # Phase-dependent playback speed (the movie stitches at a constant fps, so retime
    # by repeating / dropping frames): pre-beam + beam-on play 2x SLOWER (each frame
    # held twice); post-beam plays 1.5x FASTER (keep 2 of every 3 frames).
    playback, post_i = [], 0
    for fr in frames:
        if fr[0] <= T_OFF:                     # pre-beam + beam-on
            playback += [fr, fr]
        else:                                  # post-beam
            if post_i % 3 != 2:                # drop every 3rd → 1.5x faster
                playback.append(fr)
            post_i += 1
    print(f"  retimed: {len(frames)} -> {len(playback)} playback frames "
          f"(pre/beam 2x slower, post-beam 1.5x faster)")

    args = []
    for k, (t, xn) in enumerate(playback):
        ni, nn, v, Ti = primitives(xn, phi_g)
        args.append((k, xMm, t, Ti, v, ni, nn, ds, phi_g, kd, kp, xticks,
                     seps, labels, apex_x, xlabel, vlabel, title_prefix))
    save_frames_parallel(_render_frame, args, out, fps=30)


if __name__ == "__main__":
    main()
