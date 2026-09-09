#!/usr/bin/env python3
"""Twelve-panel movie / publication snapshot of an EXPERIMENTAL two-temperature run.

Reads the `<run>.twotemp` sidecar written by the `model_column_2t` scenario
(columns: rho v T_e T_i x_eq n_e n_HI p_total p_e p_i kappa_e kappa_i Q_ei tau_eq,
all SI, `Q_ei > 0` heats the electrons, `tau_eq` in seconds) and, optionally, the
release single-temperature `<run>.gamma_diag` sidecar, whose `T`, `v` and `rho`
are overlaid for comparison.  The `.gamma_diag` velocity column is named `v_cm`
for historical reasons but holds SI m/s (see docs/doxygen/pages/io_formats.md),
so no unit conversion is applied to either file.

The panels answer the study's question directly:
  a,b  T_e and T_i (plus the release T), full column and zoomed on the refined top
  c    the decoupling T_e - T_i in the zoom, against the sidecar's own
       ten-significant-digit print resolution, which is the honest noise floor
  d    |T_e - T_i|/T_i over the whole column, log, with that same floor drawn
  e,f  the collisional exchange Q_ei: |Q_ei| over the column against the exchange
       a print-floor temperature difference alone would produce, then signed in the zoom
  g,h  the equilibration time tau_eq and the two conduction channels kappa_e, kappa_i
  i-l  hydro context: rho, v, x, and the p_total/p_e/p_i split

Panels c, d, e, f and k also carry a dashed marker at the single physical
T_e = T_i sign change, which sits at the ionization front: below it the electron
pool is paying for ionization and runs cold, above it the electrons are
conductively heated from the outer wall and hand energy to the heavies.

The model_column_2t domain is only 1600-2153 km, so a linear height axis is
honest and no split/compressed axis is used (that convention is for runs that
span into the corona).  Because the mesh is R4-refined near the top and the
decoupling lives there, the top of the domain gets its own zoom panels and is
shaded on every full-column panel.

Usage:
  .venv/bin/python util/animate_two_temp.py <run.txt.twotemp> [options]

    --ref PATH        release .gamma_diag sidecar to overlay (matched by nearest time)
    --out PATH        movie path, or PNG path with --snapshot
    --fps N           movie frame rate (default 12)
    --snapshot        write a publication PNG + PDF of one frame instead of a movie
    --time T          frame to snapshot, in physical seconds (default: last frame)
    --stride N        use every Nth frame of the movie (default 1)
    --max-frames N    subsample the movie down to N frames
    --zoom LO HI      zoom-panel height window in km (default 2080 2153)
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/chromosphere2026-mpl")

import numpy as np
import matplotlib

matplotlib.use("Agg")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _anim_parallel import save_frames_parallel  # noqa: E402

# `.twotemp` column layout, from chromo_main.cpp's two_temp_diag block.
TT = dict(rho=0, v=1, Te=2, Ti=3, x=4, ne=5, nHI=6, p=7, pe=8, pi=9,
          kappa_e=10, kappa_i=11, Q_ei=12, tau_eq=13)
# Release `.gamma_diag` column layout, same file shape, ten columns.
GD = dict(rho=0, v=1, T=2, x=3, ne=4, nHI=5, p=6, gamma1=7, kappa=8)

DEFAULT_ZOOM = (2080.0, 2153.0)
OUT_DIR = Path("visualization/model_column_2t")

C_E = "#c0392b"      # electrons
C_I = "#1f6fb2"      # heavy particles (protons + neutral H)
C_REF = "#2f3640"    # release single-temperature run
C_DELTA = "#8e44ad"
C_Q = "#d35400"
C_TAU = "#16a085"
C_NEUTRAL = "#2f3640"
C_ZOOM = "#f0c419"
C_FLOOR = "#95a5a6"

# chromo_main writes the sidecar with ostream precision(10), i.e. ten significant
# digits, so a temperature DIFFERENCE below one unit in the tenth digit of the
# temperature itself is print round-off and carries no physics.
PRINT_SIG_DIGITS = 10
K_B = 1.380649e-23  # [J K^-1], for g_ei = 1.5 n_e k_B / tau_eq


def _apply_style(base_font=9.0):
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": base_font,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_diag(path, ncol_expected=None):
    """Stream a `.twotemp` / `.gamma_diag` sidecar.

    Returns (height_km, meta, frames) with frames a list of (t, step, ns x ncol).
    A trailing partially written frame is dropped, so this is safe to run against
    a file a live simulation is still appending to.
    """
    with open(path) as stream:
        raw = [line for line in stream.read().splitlines() if line.strip()]
    if len(raw) < 2:
        raise ValueError(f"{path}: no header")
    ns, ncol = (int(value) for value in raw[0].split())
    if ncol_expected is not None and ncol != ncol_expected:
        raise ValueError(f"{path}: expected {ncol_expected} columns, found {ncol}")
    height_km = np.fromstring(raw[1], sep=" ")
    if height_km.size != ns:
        raise ValueError(f"{path}: header says ns={ns}, heights row has {height_km.size}")

    meta, frames, i = {}, [], 2
    while i < len(raw):
        line = raw[i]
        if line.startswith("# t ="):
            parts = line.split()
            t, step = float(parts[3]), int(parts[6])
            i += 1
            if i + ns > len(raw):
                break                       # truncated trailing frame
            rows, complete = [], True
            for j in range(ns):
                row = np.fromstring(raw[i + j], sep=" ")
                if row.size != ncol:
                    complete = False
                    break
                rows.append(row)
            if not complete:
                break
            frames.append((t, step, np.asarray(rows)))
            i += ns
            continue
        if line.startswith("#") and "=" in line:
            key, _, value = line[1:].strip().partition("=")
            meta[key.strip()] = value.strip()
        i += 1

    if not frames:
        raise ValueError(f"{path}: no complete frames")
    return height_km, meta, frames


def match_reference(times, ref_times):
    """Nearest-time reference frame index for each 2T frame, plus the worst offset."""
    ref_times = np.asarray(ref_times)
    index = np.array([int(np.abs(ref_times - t).argmin()) for t in times])
    return index, float(np.max(np.abs(ref_times[index] - np.asarray(times))))


# ---------------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------------

def print_floor(temperature):
    """One unit in the last printed significant digit of `temperature` [K]."""
    safe = np.maximum(np.abs(temperature), 1.0)
    return 10.0 ** (np.floor(np.log10(safe)) - (PRINT_SIG_DIGITS - 1))


def cell_widths_m(height_km):
    """Cell widths [m] from the sidecar's cumulative upper-face heights [km]."""
    first = height_km[0] - (height_km[1] - height_km[0])
    return np.diff(height_km, prepend=first) * 1.0e3


def derived(block):
    """Everything the panels need beyond the raw columns, for one frame."""
    Te, Ti = block[:, TT["Te"]], block[:, TT["Ti"]]
    delta = Te - Ti
    floor = print_floor(np.maximum(Te, Ti))
    # Q_ei = g_ei (T_i - T_e) and tau_eq = 1.5 n_e k_B / g_ei, both printed, so the
    # exchange a print-floor temperature difference alone would produce is exact.
    g_ei = 1.5 * block[:, TT["ne"]] * K_B / block[:, TT["tau_eq"]]
    return {
        "Te": Te, "Ti": Ti, "delta": delta, "floor": floor,
        "relative": delta / Ti, "rel_floor": floor / Ti,
        "Q": block[:, TT["Q_ei"]], "Q_floor": g_ei * floor,
    }


def sound_crossing_s(block, ds_m):
    """Cell acoustic crossing time ds / sqrt(5/3 p/rho) [s].

    The release/2T CFL step is a fraction of this, so it is the timescale tau_eq
    has to be compared against to judge whether the exchange is stiff.
    """
    c_s = np.sqrt((5.0 / 3.0) * block[:, TT["p"]] / block[:, TT["rho"]])
    return ds_m / c_s


def crossing_km(height_km, delta, floor):
    """Height [km] of the single physical T_e = T_i sign change, or None.

    Taken as the zero crossing BETWEEN the most negative and the most positive
    cell, which in this column is the ionization front: below it the electron
    pool is paying for ionization and T_e < T_i, above it the electrons are
    conductively heated from the outer wall and hand energy to the heavies.
    Sign flips elsewhere are print-floor flicker, so the whole marker is
    suppressed unless the profile is an order of magnitude above that floor.
    """
    if np.nanmax(np.abs(delta)) < 10.0 * float(np.median(floor)):
        return None
    lo, hi = int(np.argmin(delta)), int(np.argmax(delta))
    if delta[lo] >= 0.0 or delta[hi] <= 0.0:
        return None
    a, b = (lo, hi) if lo < hi else (hi, lo)
    flips = np.flatnonzero(np.diff(np.sign(delta[a:b + 1])) != 0)
    if flips.size == 0:
        return None
    i = a + int(flips[0])
    left, right = delta[i], delta[i + 1]
    weight = left / (left - right) if left != right else 0.5
    return float(height_km[i] + weight * (height_km[i + 1] - height_km[i]))


def _mark_crossing(ax, height, label=False):
    """Vertical marker at the T_e = T_i sign change."""
    if height is None:
        return
    ax.axvline(height, color="#34495e", lw=0.9, ls=(0, (4, 2)), zorder=3,
               label=(rf"$T_e = T_i$ at {height:.1f} km" if label else None))


def fmt_time(seconds):
    """Human time label that does not round a sub-second frame away to 0.0."""
    return f"{seconds:.1f}" if abs(seconds) >= 1.0 else f"{seconds:.4g}"


def _time_tag(seconds):
    """Filename-safe time label; keeps sub-second frames distinguishable."""
    if seconds >= 1.0:
        return f"{seconds:.0f}s"
    return ("%.3gs" % seconds).replace(".", "p")


def _positive(values):
    """Copy with non-positive / non-finite entries blanked, for log axes."""
    out = np.asarray(values, dtype=float).copy()
    out[~np.isfinite(out) | (out <= 0.0)] = np.nan
    return out


def _sym_limits(arrays, pad=1.35, floor=1e-30):
    top = max(float(np.max(np.abs(a[np.isfinite(a)]))) for a in arrays)
    top = max(top, floor)
    return -top * pad, top * pad


def _log_limits(arrays, lo_pad=0.6, hi_pad=1.7):
    lo = min(float(np.min(a[np.isfinite(a) & (a > 0)])) for a in arrays)
    hi = max(float(np.max(a[np.isfinite(a) & (a > 0)])) for a in arrays)
    return lo * lo_pad, hi * hi_pad


def _lin_limits(arrays, pad=0.06):
    lo = min(float(np.min(a[np.isfinite(a)])) for a in arrays)
    hi = max(float(np.max(a[np.isfinite(a)])) for a in arrays)
    span = hi - lo if hi > lo else max(abs(hi), 1.0)
    return lo - pad * span, hi + pad * span


def build_context(height_km, meta, blocks, ref_blocks, zoom, label):
    """Global y-limits and symlog thresholds, so nothing jumps between frames."""
    zoom_mask = (height_km >= zoom[0]) & (height_km <= zoom[1])
    ds_m = cell_widths_m(height_km)
    der = [derived(b) for b in blocks]
    delta = [d["delta"] for d in der]
    t_sound = [sound_crossing_s(b, ds_m) for b in blocks]

    temps = [b[:, TT["Te"]] for b in blocks] + [b[:, TT["Ti"]] for b in blocks]
    if ref_blocks:
        temps += [r[:, GD["T"]] for r in ref_blocks]
    vels = [b[:, TT["v"]] for b in blocks]
    if ref_blocks:
        vels += [r[:, GD["v"]] for r in ref_blocks]

    def floored_log(values, floors, lo_pad=0.2, hi_pad=3.0):
        """Log limits that always keep the print-resolution floor in view."""
        hi = max(float(np.nanmax(np.abs(a))) for a in values)
        lo = min(float(np.nanmin(f[f > 0])) for f in floors)
        return lo * lo_pad, max(hi, lo) * hi_pad

    q = [b[:, TT["Q_ei"]] for b in blocks]
    lims = {
        "T": _lin_limits(temps),
        "T_zoom": _lin_limits([a[zoom_mask] for a in temps]),
        "delta_zoom": _lin_limits([d[zoom_mask] for d in delta], pad=0.15),
        "rel": floored_log([d["relative"] for d in der],
                           [d["rel_floor"] for d in der]),
        "Q": floored_log([d["Q"] for d in der], [d["Q_floor"] for d in der]),
        "Q_zoom": _sym_limits([a[zoom_mask] for a in q]),
        "tau": _log_limits([b[:, TT["tau_eq"]] for b in blocks] + t_sound),
        "kappa": _log_limits([b[:, TT["kappa_e"]] for b in blocks]
                             + [b[:, TT["kappa_i"]] for b in blocks]),
        "rho": _log_limits([b[:, TT["rho"]] for b in blocks]
                           + ([r[:, GD["rho"]] for r in ref_blocks] if ref_blocks else [])),
        "v": _lin_limits(vels),
        "x": _lin_limits([b[:, TT["x"]] for b in blocks], pad=0.08),
        "p": _log_limits([b[:, TT["p"]] for b in blocks]
                         + [b[:, TT["pe"]] for b in blocks]
                         + [b[:, TT["pi"]] for b in blocks]),
    }
    # Symlog linear band of the signed zoom panel: the print-resolution exchange,
    # so everything that collapses into the band is round-off by construction.
    lin_q = max(float(np.median(np.concatenate(
        [d["Q_floor"][zoom_mask] for d in der]))), 1e-30)
    return {
        "h": height_km,
        "ds_m": ds_m,
        "zoom": zoom,
        "lims": lims,
        "linthresh": {"Q_zoom": lin_q},
        "label": label,
        "bc": meta.get("outer_Ti_neumann", "?"),
        "has_ref": bool(ref_blocks),
    }


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------

def _decorate(ax, ctx, xlim, ylabel, title, letter, zoom_span=True):
    if zoom_span and xlim[0] < ctx["zoom"][0]:
        ax.axvspan(*ctx["zoom"], color=C_ZOOM, alpha=0.12, lw=0, zorder=0)
    ax.set_xlim(*xlim)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize="small", loc="left", color="#57606f")
    ax.grid(axis="both", color="0.9", lw=0.5, zorder=0)
    ax.text(0.012, 0.955, letter, transform=ax.transAxes, va="top",
            fontweight="bold", fontsize="medium")


def draw_figure(ctx, payload, figsize=(17.0, 10.0), font=9.0):
    """Render one frame; returns the matplotlib figure."""
    import matplotlib.pyplot as plt

    _apply_style(font)
    h = ctx["h"]
    zoom = ctx["zoom"]
    full = (h[0], h[-1])
    lims, lt = ctx["lims"], ctx["linthresh"]
    d = payload["block"]
    ref = payload.get("ref")

    dv = derived(d)
    Te, Ti, delta, floor = dv["Te"], dv["Ti"], dv["delta"], dv["floor"]
    relative = dv["relative"]
    cross = crossing_km(h, delta, floor)

    fig, axes = plt.subplots(3, 4, figsize=figsize)

    # -- a, b: the two temperatures, full column and refined top -------------
    for col, (ax, xlim, key) in enumerate((
            (axes[0, 0], full, "T"), (axes[0, 1], zoom, "T_zoom"))):
        if ref is not None:
            ax.plot(h, ref[:, GD["T"]] * 1e-3, color=C_REF, lw=1.6, ls="--",
                    alpha=0.75, label=r"release 1T $T$")
        ax.plot(h, Ti * 1e-3, color=C_I, lw=1.6, label=r"$T_i$ (heavy)")
        ax.plot(h, Te * 1e-3, color=C_E, lw=1.0, label=r"$T_e$")
        _decorate(ax, ctx, xlim, r"$T$  [kK]",
                  "full column" if col == 0 else f"zoom {zoom[0]:.0f}-{zoom[1]:.0f} km",
                  "ab"[col])
        ax.set_ylim(lims[key][0] * 1e-3, lims[key][1] * 1e-3)
        if col == 0:
            ax.legend(loc="upper left", fontsize="small", bbox_to_anchor=(0.06, 0.98))

    # -- c: the decoupling itself, against the print resolution --------------
    ax = axes[0, 2]
    ax.fill_between(h, -floor, floor, color=C_FLOOR, alpha=0.30, lw=0,
                    label="sidecar print resolution")
    ax.axhline(0.0, color="#95a5a6", lw=0.6)
    ax.plot(h, delta, color=C_DELTA, lw=1.4)
    _mark_crossing(ax, cross, label=True)
    _decorate(ax, ctx, zoom, r"$T_e - T_i$  [K]",
              f"decoupling, zoom {zoom[0]:.0f}-{zoom[1]:.0f} km", "c")
    ax.set_ylim(*lims["delta_zoom"])
    ax.legend(loc="lower left", fontsize="small")

    # -- d: relative decoupling over the whole column ------------------------
    ax = axes[0, 3]
    ax.fill_between(h, lims["rel"][0], dv["rel_floor"], color=C_FLOOR,
                    alpha=0.25, lw=0, label="below print resolution")
    ax.semilogy(h, _positive(np.abs(relative)), color=C_DELTA, lw=1.2)
    _mark_crossing(ax, cross)
    _decorate(ax, ctx, full, r"$|T_e - T_i| / T_i$",
              "relative decoupling, full column", "d")
    ax.set_ylim(*lims["rel"])
    ax.legend(loc="upper left", fontsize="small", bbox_to_anchor=(0.06, 0.99))

    # -- e: exchange magnitude against the exchange the print floor alone gives
    ax = axes[1, 0]
    ax.fill_between(h, lims["Q"][0], _positive(dv["Q_floor"]), color=C_FLOOR,
                    alpha=0.25, lw=0, label="below print resolution")
    ax.semilogy(h, _positive(np.abs(dv["Q"])), color=C_Q, lw=1.2)
    _mark_crossing(ax, cross)
    _decorate(ax, ctx, full, r"$|Q_{ei}|$  [W m$^{-3}$]",
              "collisional exchange magnitude, full column", "e")
    ax.set_ylim(*lims["Q"])
    ax.legend(loc="upper left", fontsize="small", bbox_to_anchor=(0.06, 0.99))

    # -- f: signed exchange where it is real --------------------------------
    ax = axes[1, 1]
    ax.set_yscale("symlog", linthresh=lt["Q_zoom"])
    ax.fill_between(h, -dv["Q_floor"], dv["Q_floor"], color=C_FLOOR, alpha=0.25,
                    lw=0, label="print resolution")
    ax.axhline(0.0, color="#95a5a6", lw=0.6)
    ax.plot(h, dv["Q"], color=C_Q, lw=1.2)
    _mark_crossing(ax, cross)
    _decorate(ax, ctx, zoom, r"$Q_{ei}$  [W m$^{-3}$]",
              f"signed exchange ($>0$ heats $e^-$), zoom "
              f"{zoom[0]:.0f}-{zoom[1]:.0f} km", "f")
    ax.set_ylim(*lims["Q_zoom"])
    ax.legend(loc="upper left", fontsize="small", bbox_to_anchor=(0.06, 0.99))

    # -- g: equilibration time against the timescale it must beat ------------
    ax = axes[1, 2]
    ax.semilogy(h, sound_crossing_s(d, ctx["ds_m"]), color=C_REF, lw=1.3, ls="--",
                label=r"cell crossing $\Delta s / c_s$")
    ax.semilogy(h, d[:, TT["tau_eq"]], color=C_TAU, lw=1.6,
                label=r"$\tau_{\rm eq}$")
    _decorate(ax, ctx, full, "time  [s]",
              "electron-heavy equilibration vs. the hydro timescale", "g")
    ax.set_ylim(*lims["tau"])
    ax.legend(loc="center left", fontsize="small", bbox_to_anchor=(0.06, 0.58))

    # -- h: the two conduction channels -------------------------------------
    ax = axes[1, 3]
    ax.semilogy(h, d[:, TT["kappa_i"]], color=C_I, lw=1.4,
                label=r"$\kappa_i$ (neutral H, on $T_i$)")
    ax.semilogy(h, d[:, TT["kappa_e"]], color=C_E, lw=1.4,
                label=r"$\kappa_e$ (Spitzer, on $T_e$)")
    _decorate(ax, ctx, full, r"$\kappa$  [W m$^{-1}$ K$^{-1}$]",
              "conduction channels", "h")
    ax.set_ylim(*lims["kappa"])
    ax.legend(loc="lower left", fontsize="small")

    # -- i, j, k, l: hydro context ------------------------------------------
    ax = axes[2, 0]
    if ref is not None:
        ax.semilogy(h, ref[:, GD["rho"]], color=C_REF, lw=1.6, ls="--", alpha=0.75)
    ax.semilogy(h, d[:, TT["rho"]], color=C_NEUTRAL, lw=1.2)
    _decorate(ax, ctx, full, r"$\rho$  [kg m$^{-3}$]", "mass density", "i")
    ax.set_ylim(*lims["rho"])

    ax = axes[2, 1]
    v = d[:, TT["v"]]
    ax.axhline(0.0, color="#95a5a6", lw=0.6)
    ax.fill_between(h, 0.0, v, where=(v >= 0), color=C_I, alpha=0.20, lw=0)
    ax.fill_between(h, 0.0, v, where=(v < 0), color=C_E, alpha=0.20, lw=0)
    if ref is not None:
        ax.plot(h, ref[:, GD["v"]], color=C_REF, lw=1.4, ls="--", alpha=0.75,
                label="release 1T")
        ax.legend(loc="upper left", fontsize="small")
    ax.plot(h, v, color=C_I, lw=1.2)
    _decorate(ax, ctx, full, r"$V$  [m s$^{-1}$]", "field-aligned velocity", "j")
    ax.set_ylim(*lims["v"])

    ax = axes[2, 2]
    ax.plot(h, d[:, TT["x"]], color=C_NEUTRAL, lw=1.4, label=r"$x_{\rm Saha}(T_e)$")
    ax.plot(h, d[:, TT["pe"]] / d[:, TT["p"]], color=C_E, lw=1.2, ls=":",
            label=r"$p_e/p$")
    _mark_crossing(ax, cross)
    _decorate(ax, ctx, full, "fraction", "ionization and electron pressure share", "k")
    ax.set_ylim(*lims["x"])
    ax.legend(loc="upper left", fontsize="small", bbox_to_anchor=(0.06, 0.99))

    ax = axes[2, 3]
    ax.semilogy(h, d[:, TT["p"]], color=C_NEUTRAL, lw=1.4, label=r"$p$")
    ax.semilogy(h, d[:, TT["pi"]], color=C_I, lw=1.2, label=r"$p_i$")
    ax.semilogy(h, d[:, TT["pe"]], color=C_E, lw=1.2, label=r"$p_e$")
    _decorate(ax, ctx, full, r"$p$  [Pa]", "pressure split", "l")
    ax.set_ylim(*lims["p"])
    ax.legend(loc="lower left", fontsize="small", ncol=3)

    for ax in axes[2, :]:
        ax.set_xlabel("height  [km]")

    resolved = np.abs(delta) > floor
    peak = int(np.abs(delta).argmax())
    headline = (
        rf"$\max|T_e-T_i| = {np.abs(delta[peak]):.3g}$ K at $h = {h[peak]:.1f}$ km"
        rf"   |   $\max|T_e-T_i|/T_i = {np.abs(relative).max():.2e}$"
        rf"   |   {int(resolved.sum())}/{delta.size} cells above print resolution"
    )
    if cross is not None:
        headline += rf"   |   $T_e = T_i$ crossing at $h = {cross:.1f}$ km"
    ref_note = ""
    if ref is not None:
        ref_note = (rf"   |   release $T$ overlaid at "
                    rf"$t = {fmt_time(payload['ref_t'])}$ s")
    fig.suptitle(
        f"EXPERIMENTAL two-temperature column  |  {ctx['label']}  |  "
        f"outer $T_i$ BC: {'Neumann' if ctx['bc'] == '1' else 'Dirichlet'}  |  "
        f"t = {fmt_time(payload['t'])} s, step {payload['step']}"
        f"{payload.get('counter','')}"
        "\n" + headline + ref_note,
        fontsize="medium", y=0.995, va="top",
    )
    fig.tight_layout(h_pad=1.2, w_pad=1.8)
    return fig


def _render_frame(k, tmpdir, ctx, payload):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = draw_figure(ctx, payload)
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("twotemp", type=Path, help="<run>.twotemp sidecar")
    parser.add_argument("--ref", type=Path, default=None,
                        help="release <run>.gamma_diag sidecar to overlay")
    parser.add_argument("--out", type=Path, default=None,
                        help="movie path, or PNG path with --snapshot")
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--snapshot", action="store_true",
                        help="write one publication-quality PNG + PDF, no movie")
    parser.add_argument("--time", type=float, default=None,
                        help="snapshot time in seconds (default: last frame)")
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--max-frames", type=int, default=None)
    parser.add_argument("--zoom", type=float, nargs=2, default=list(DEFAULT_ZOOM),
                        metavar=("LO", "HI"))
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    height_km, meta, frames = read_diag(args.twotemp, ncol_expected=14)
    label = os.path.basename(str(args.twotemp)).replace(".txt.twotemp", "")
    print(f"{args.twotemp}: {len(frames)} complete frames, {height_km.size} cells, "
          f"t = {frames[0][0]:.1f} -> {frames[-1][0]:.1f} s, "
          f"outer_Ti_neumann={meta.get('outer_Ti_neumann', '?')}")

    ref_frames = None
    if args.ref is not None:
        ref_h, _ref_meta, ref_frames = read_diag(args.ref, ncol_expected=10)
        if ref_h.size != height_km.size or not np.allclose(ref_h, height_km, atol=1e-3):
            sys.exit(f"{args.ref}: mesh differs from {args.twotemp}; cannot overlay")
        print(f"{args.ref}: {len(ref_frames)} frames, "
              f"t = {ref_frames[0][0]:.1f} -> {ref_frames[-1][0]:.1f} s")

    if args.snapshot:
        if args.time is None:
            pick = len(frames) - 1
        else:
            pick = int(np.abs(np.array([f[0] for f in frames]) - args.time).argmin())
        selected = [frames[pick]]
    else:
        selected = frames[::max(1, args.stride)]
        if args.max_frames and len(selected) > args.max_frames:
            take = np.linspace(0, len(selected) - 1, args.max_frames).astype(int)
            selected = [selected[i] for i in take]

    ref_blocks = None
    if ref_frames is not None:
        ref_times = [f[0] for f in ref_frames]
        index, worst = match_reference([f[0] for f in selected], ref_times)
        ref_blocks = [ref_frames[i][2] for i in index]
        ref_times_used = [ref_times[i] for i in index]
        if worst > 1.0:
            print(f"WARNING: reference frames matched by nearest time; worst "
                  f"offset {worst:.1f} s (the two runs do not share a cadence)")

    blocks = [f[2] for f in selected]
    ctx = build_context(height_km, meta, blocks, ref_blocks, tuple(args.zoom), label)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if args.snapshot:
        out = args.out or (OUT_DIR / f"{label}_t{_time_tag(selected[0][0])}.png")
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        payload = {"block": blocks[0], "t": selected[0][0], "step": selected[0][1]}
        if ref_blocks is not None:
            payload["ref"] = ref_blocks[0]
            payload["ref_t"] = ref_times_used[0]
        fig = draw_figure(ctx, payload, figsize=(17.0, 10.0), font=9.5)
        stem = Path(out).with_suffix("")
        fig.savefig(stem.with_suffix(".png"), dpi=300)
        fig.savefig(stem.with_suffix(".pdf"))
        import matplotlib.pyplot as plt
        plt.close(fig)
        block = blocks[0]
        delta = block[:, TT["Te"]] - block[:, TT["Ti"]]
        print(f"wrote {stem}.png and {stem}.pdf  "
              f"(t = {fmt_time(selected[0][0])} s, "
              f"max|T_e-T_i| = {np.abs(delta).max():.4g} K, "
              f"max relative = {np.abs(delta / block[:, TT['Ti']]).max():.4e})")
        return

    out = args.out or (OUT_DIR / f"{label}.mp4")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    frame_args = []
    for k, (t, step, block) in enumerate(selected):
        payload = {"block": block, "t": t, "step": step,
                   "counter": f"  ({k + 1}/{len(selected)})"}
        if ref_blocks is not None:
            payload["ref"] = ref_blocks[k]
            payload["ref_t"] = ref_times_used[k]
        frame_args.append((k, ctx, payload))
    save_frames_parallel(_render_frame, frame_args, str(out), args.fps)


if __name__ == "__main__":
    main()
