#!/usr/bin/env python3
"""Upper-boundary conduction-wall temperature sensitivity for model_column.

Reads the ``<run>.gamma_diag`` sidecar (columns
``rho v T x_eq n_e n_HI p Gamma1 kappa_physical kappa_solver``) and the
``<run>.outercond`` sidecar (``t step T_top T_wall kappa_face ds_face
area_ratio q_face``) of a sweep produced by ``util/run_twall_sweep.sh``, and
answers one question: how does the field-aligned velocity respond to the imposed
conduction-wall temperature T_wall, with the C7 initial condition and every other
piece of physics and numerics held fixed?

Products (into visualization/model_column/ by default):

  twall_V_of_t.png        V(t) at representative upper-chromosphere/TR heights
  twall_V_of_h.png        quasi-steady V(h), full column and top-of-domain zoom
  twall_summary.png       quasi-steady V, mass flux, q_wall and T_top vs T_wall
  twall_support.png       mass flux rho*V(h), T(h), column mass drift

and a JSON of every scalar metric.

"Quasi-steady" means the time average over the last ``--avg-window`` seconds of
the run, which averages over the residual acoustic ringing without assuming the
column has reached a true steady state (it has not at 4000 s; see the column-mass
drift panel).

Usage:
    python util/twall_sweep_diag.py --dir outputs/model_column/twall_sweep \
        [--cfl 0.50] [--heights 1800 2000 2100] [--avg-window 500] \
        [--out-dir visualization/model_column] [--json metrics.json]
"""
import argparse
import glob
import json
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

COLUMNS = ["rho", "v", "T", "x_eq", "n_e", "n_HI", "p", "Gamma1",
           "kappa_physical", "kappa_solver"]
IV = COLUMNS.index("v")
IRHO = COLUMNS.index("rho")
IT = COLUMNS.index("T")
IP = COLUMNS.index("p")
IX = COLUMNS.index("x_eq")

CHI_H_J = 2.179872361e-18     # hydrogen ionization energy [J]
M_H_KG = 1.6726219e-27        # hydrogen mass [kg]
GAMMA_FROZEN = 5.0 / 3.0      # the release Riemann solver's frozen index


# --------------------------------------------------------------------------
# loaders (the .gamma_diag reader follows util/upper_bc_decouple_diag.py)
# --------------------------------------------------------------------------
def load_gamma_diag(path):
    """-> (heights_km, times[nt], data[nt, ns, ncol]); stops at the first bad frame."""
    diag = path if path.endswith(".gamma_diag") else path + ".gamma_diag"
    with open(diag) as handle:
        ns, ncol = map(int, handle.readline().split())
        h = np.array([float(x) for x in handle.readline().split()])
        frames, rows, t = [], None, None

        def flush():
            if rows is None or len(rows) != ns:
                return False
            arr = np.array(rows)
            if not np.isfinite(arr).all():
                return False
            frames.append((t, arr))
            return True

        for line in handle:
            if line.startswith("# t ="):
                if rows is not None and not flush():
                    break
                t = float(line.split("t =")[1].split("step")[0])
                rows = []
            elif line.startswith("#"):
                continue
            elif rows is not None:
                parts = line.split()
                if len(parts) == ncol:
                    rows.append([float(x) for x in parts])
        else:
            flush()
    if not frames:
        raise SystemExit(f"no finite frames in {diag}")
    times = np.array([f[0] for f in frames])
    data = np.stack([f[1] for f in frames])
    return h, times, data


def load_outercond(path):
    """-> dict of the outer-face conduction sidecar columns."""
    oc = path if path.endswith(".outercond") else path + ".outercond"
    if not os.path.exists(oc):
        return None
    keys = ["t", "step", "T_top", "T_wall", "kappa_face", "ds_face",
            "area_ratio", "q_face"]
    rows = []
    with open(oc) as handle:                    # tolerate a truncated last line
        for line in handle:                     # (the sidecar of a live run)
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != len(keys):
                continue
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    if not rows:
        return None
    raw = np.array(rows)
    return {k: raw[:, i] for i, k in enumerate(keys)}


def console_ok(path):
    """-> (terminated_at_end_time, godunov_fallbacks) from the run console log."""
    log = path.replace(".txt", ".console.log")
    if not os.path.exists(log):
        return None, None
    text = open(log, errors="replace").read()
    term = re.search(r"^termination=(\S+)", text, re.M)
    fb = re.search(r"^godunov\.fallbacks=(\d+)", text, re.M)
    return (term.group(1) if term else None, int(fb.group(1)) if fb else None)


# --------------------------------------------------------------------------
def cell_weights(h_km):
    """Cell widths [m] from the cell-centre heights (works on the refined mesh)."""
    edges = np.empty(len(h_km) + 1)
    edges[1:-1] = 0.5 * (h_km[1:] + h_km[:-1])
    edges[0] = h_km[0] - (edges[1] - h_km[0])
    edges[-1] = h_km[-1] + (h_km[-1] - edges[-2])
    return np.diff(edges) * 1.0e3


def analyse(label, t_wall, path, heights, avg_window):
    h, times, data = load_gamma_diag(path)
    ds = cell_weights(h)
    oc = load_outercond(path)
    term, fallbacks = console_ok(path)

    v = data[:, :, IV]                      # [nt, ns] m/s
    rho = data[:, :, IRHO]
    mass_flux = rho * v                     # kg m^-2 s^-1
    mass = (rho * ds[None, :]).sum(axis=1)  # kg m^-2

    window = times >= (times[-1] - avg_window)
    # The preceding equal-length window: the fractional change between the two
    # is how far from a true steady state the run still is at t_end. The hot-wall
    # cases are still relaxing downward at 4000 s, so this number has to be
    # reported alongside every "quasi-steady" value.
    prev = (times >= (times[-1] - 2.0 * avg_window)) & ~window
    v_qs = v[window].mean(axis=0)
    rho_qs_top = float(rho[window, -1].mean())
    f_qs = mass_flux[window].mean(axis=0)
    T_qs = data[window, :, IT].mean(axis=0)
    p_qs = data[window, :, IP].mean(axis=0)

    # Outer-BC health: ISO_VCAP caps the outer GHOST velocity at 0.1 c_s. If the
    # live top-cell |V| approaches that cap the outflow face is constraining the
    # solution, and any apparent saturation would be the boundary, not physics.
    g1 = data[:, -1, COLUMNS.index("Gamma1")]
    c_s = np.sqrt(g1 * data[:, -1, IP] / data[:, -1, IRHO])
    vcap_frac = np.abs(v[:, -1]) / (0.1 * c_s)

    idx = {hh: int(np.argmin(np.abs(h - hh))) for hh in heights}
    out = {
        "label": label, "t_wall_K": t_wall, "path": path,
        "termination": term, "godunov_fallbacks": fallbacks,
        "ns": len(h), "t_end": float(times[-1]), "n_frames": len(times),
        "avg_window_s": avg_window,
        "V_top_qs": float(v_qs[-1]),
        "V_top_prev_window": float(v[prev, -1].mean()) if prev.any() else None,
        "V_top_drift_frac": (float((v_qs[-1] - v[prev, -1].mean()) / v_qs[-1])
                             if prev.any() and v_qs[-1] != 0.0 else None),
        "V_top_final": float(v[-1, -1]),
        "F_top_qs": float(f_qs[-1]),
        "T_top_qs": float(T_qs[-1]),
        "p_top_qs": float(p_qs[-1]),
        "V_max_qs": float(v_qs[np.argmax(np.abs(v_qs))]),
        "h_V_max_qs": float(h[np.argmax(np.abs(v_qs))]),
        "cs_top_qs": float(c_s[window].mean()),
        "mach_top_qs": float((v[window, -1] / c_s[window]).mean()),
        "vcap_frac_qs": float(vcap_frac[window].mean()),
        "vcap_frac_max": float(vcap_frac.max()),
        "mass_0": float(mass[0]), "mass_end": float(mass[-1]),
        "mass_rel_drift": float(mass[-1] / mass[0] - 1.0),
        "at_height": {str(hh): {"h_actual_km": float(h[i]),
                                "V_qs": float(v_qs[i]),
                                "F_qs": float(f_qs[i]),
                                "T_qs": float(T_qs[i])}
                      for hh, i in idx.items()},
    }
    if oc is not None:
        w = oc["t"] >= (oc["t"][-1] - avg_window)
        out["q_wall_qs"] = float(oc["q_face"][w].mean())
        out["q_wall_t0"] = float(oc["q_face"][0])
        out["kappa_wall_qs"] = float(oc["kappa_face"][w].mean())
        out["dT_wall_top_qs"] = float((oc["T_wall"] - oc["T_top"])[w].mean())
        # Steady-state enthalpy budget: whatever the outer face conducts in has
        # nowhere to go in the release (no radiative sink, no coronal heating),
        # so a quasi-steady column must export it as the enthalpy flux of the
        # evaporating material. The ionization term x*chi_H/m_H = 1.3e9 J/kg at
        # x = 1 dominates the translational term 2.5 p/rho ~ 9e5 J/kg by three
        # orders of magnitude, so the prediction is essentially
        #     V ~ q_wall / (rho_top * x_top * chi_H/m_H).
        x_top = float(data[window, -1, IX].mean())
        h_therm = GAMMA_FROZEN / (GAMMA_FROZEN - 1.0) * p_qs[-1] / rho_qs_top
        h_ion = x_top * CHI_H_J / M_H_KG
        out["x_top_qs"] = x_top
        out["h_thermal_top"] = float(h_therm)
        out["h_ion_top"] = float(h_ion)
        out["rho_top_qs"] = float(rho_qs_top)
        out["V_pred_enthalpy"] = float(out["q_wall_qs"]
                                       / (rho_qs_top * (h_therm + h_ion)))
        out["V_ratio_pred"] = float(out["V_top_qs"] / out["V_pred_enthalpy"])
    return out, dict(h=h, times=times, v=v, rho=rho, mass_flux=mass_flux,
                     mass=mass, v_qs=v_qs, f_qs=f_qs, T_qs=T_qs, oc=oc, idx=idx)


# --------------------------------------------------------------------------
def wall_colors(t_walls):
    """Perceptually monotone colour per T_wall on a log scale.

    Not a diverging map: the mid-range of coolwarm is near-white and the 18-30 kK
    curves — the ones nearest the 22 kK release baseline — became invisible.
    """
    norm = mcolors.LogNorm(vmin=min(t_walls), vmax=max(t_walls))
    cmap = plt.get_cmap("plasma")
    return {t: cmap(0.9 * norm(t)) for t in t_walls}, norm, cmap


def plot_v_of_t(runs, heights, out_png, avg_window):
    t_walls = [r["meta"]["t_wall_K"] for r in runs]
    col, norm, cmap = wall_colors(t_walls)
    panels = list(heights) + ["top"]
    fig, axes = plt.subplots(len(panels), 1, figsize=(9.0, 2.5 * len(panels)),
                             sharex=True)
    for ax, panel in zip(np.atleast_1d(axes), panels):
        for r in runs:
            d, m = r["data"], r["meta"]
            if panel == "top":
                series, hh = d["v"][:, -1], d["h"][-1]
            else:
                i = d["idx"][panel]
                series, hh = d["v"][:, i], d["h"][i]
            ax.plot(d["times"], series, color=col[m["t_wall_K"]], lw=1.2,
                    label=f"{m['t_wall_K']/1e3:g} kK")
        ax.axhline(0.0, color="0.4", lw=0.7, ls=":")
        ax.axvspan(runs[0]["data"]["times"][-1] - avg_window,
                   runs[0]["data"]["times"][-1], color="0.85", zorder=0)
        # symlog: the 8-15 kK responses are sub-m/s while the 100 kK top cell
        # starts at 3 km/s, so a linear axis shows only the hot walls.
        ax.set_yscale("symlog", linthresh=1.0)
        ax.set_ylabel("V [m/s]")
        ax.set_title(f"h = {hh:.1f} km" + (" (top cell)" if panel == "top" else ""),
                     fontsize=9, loc="left")
        ax.grid(alpha=0.3, which="both")
    np.atleast_1d(axes)[-1].set_xlabel("t [s]")
    handles, labels = np.atleast_1d(axes)[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=min(7, len(labels)), fontsize=8,
               loc="lower center", title="$T_{\\rm wall}$", title_fontsize=8)
    fig.suptitle("model_column: V(t) vs imposed conduction-wall temperature "
                 "(same C7 IC, decoupled hydro ghost T)", fontsize=10)
    fig.tight_layout(rect=(0, 0.07, 1, 0.98))
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def plot_v_of_h(runs, out_png, avg_window):
    t_walls = [r["meta"]["t_wall_K"] for r in runs]
    col, _, _ = wall_colors(t_walls)
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.6))
    for r in runs:
        d, m = r["data"], r["meta"]
        c = col[m["t_wall_K"]]
        lbl = f"{m['t_wall_K']/1e3:g} kK"
        axes[0].plot(d["h"], d["v_qs"], color=c, lw=1.3, label=lbl)
        axes[1].plot(d["h"], d["v_qs"], color=c, lw=1.3)
        axes[2].plot(d["h"], np.sign(d["v_qs"]) * np.abs(d["v_qs"]),
                     color=c, lw=1.3)
    axes[0].set_title("full column", fontsize=9, loc="left")
    axes[1].set_title("top 40 km (resolved TR)", fontsize=9, loc="left")
    axes[1].set_xlim(2114, 2154)
    axes[2].set_title("full column, symlog V", fontsize=9, loc="left")
    axes[2].set_yscale("symlog", linthresh=1.0)
    for ax in axes:
        ax.axhline(0.0, color="0.4", lw=0.7, ls=":")
        ax.set_xlabel("h [km]")
        ax.set_ylabel("V [m/s]")
        ax.grid(alpha=0.3)
    axes[0].legend(ncol=2, fontsize=8, title="$T_{\\rm wall}$", title_fontsize=8)
    fig.suptitle(f"model_column: quasi-steady V(h), mean over the last "
                 f"{avg_window:g} s of a 4000 s run", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def plot_summary(metrics, heights, out_png, avg_window):
    tw = np.array([m["t_wall_K"] for m in metrics]) / 1.0e3
    order = np.argsort(tw)
    tw = tw[order]
    g = lambda key: np.array([metrics[i].get(key, np.nan) for i in order])
    gh = lambda hh, key: np.array(
        [metrics[i]["at_height"][str(hh)][key] for i in order])

    fig, axes = plt.subplots(2, 3, figsize=(16.0, 8.0))
    ax = axes[0, 0]
    ax.plot(tw, g("V_top_qs"), "o-", color="crimson", label="top cell")
    for hh, c in zip(heights, ["tab:blue", "tab:green", "tab:orange"]):
        ax.plot(tw, gh(hh, "V_qs"), "s--", color=c, ms=4, label=f"{hh:g} km")
    ax.axhline(0.0, color="0.4", lw=0.8, ls=":")
    ax.axvline(22.0, color="0.6", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=0.1)
    ax.set_xlabel("$T_{\\rm wall}$ [kK]")
    ax.set_ylabel("quasi-steady V [m/s]")
    ax.set_title("velocity response (symlog)", fontsize=9, loc="left")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    # log-log power-law view over the upflow branch, plus the LOCAL slope: the
    # single question the sweep exists to answer is whether the response is
    # linear, a power law, or saturating, and only the local slope shows that.
    ax = axes[0, 2]
    vt = g("V_top_qs")
    up = vt > 0
    ax.loglog(tw[up], vt[up], "o", color="crimson")
    fit = up & (tw >= 15.0)
    if fit.sum() >= 2:
        sl, ic = np.polyfit(np.log(tw[fit]), np.log(vt[fit]), 1)
        ax.loglog(tw[fit], np.exp(ic) * tw[fit] ** sl, "-", color="0.4",
                  label=f"$T_{{\\rm wall}}\\geq 15$ kK: "
                        f"$V \\propto T_{{\\rm wall}}^{{{sl:.2f}}}$")
        ax.legend(fontsize=8, loc="upper left")
    ax.set_xlabel("$T_{\\rm wall}$ [kK]")
    ax.set_ylabel("quasi-steady $V_{\\rm top}$ [m/s]")
    ax.set_title("upflow branch, log-log + local slope", fontsize=9, loc="left")
    ax.grid(alpha=0.3, which="both")
    if up.sum() >= 3:
        lt, lv = np.log(tw[up]), np.log(vt[up])
        n_loc = np.diff(lv) / np.diff(lt)
        t_mid = np.exp(0.5 * (lt[1:] + lt[:-1]))
        ax2 = ax.twinx()
        ax2.plot(t_mid, n_loc, "^--", color="tab:blue", ms=4, lw=1.0)
        ax2.set_ylabel(r"local slope $d\ln V/d\ln T_{\rm wall}$",
                       color="tab:blue", fontsize=9)
        ax2.tick_params(axis="y", labelcolor="tab:blue")
        ax2.axhline(1.0, color="tab:blue", lw=0.6, ls=":")
        ax2.set_ylim(0, max(6.0, float(np.nanmax(n_loc)) * 1.1))

    # enthalpy-budget closure
    ax = axes[1, 2]
    if np.isfinite(g("V_pred_enthalpy")).any():
        ax.plot(tw, g("V_top_qs"), "o-", color="crimson", label="measured")
        ax.plot(tw, g("V_pred_enthalpy"), "k--s", ms=4,
                label=r"$q_{\rm wall}/[\rho_{\rm top}(h+x\chi_H/m_H)]$")
        ax.axhline(0.0, color="0.4", lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=0.1)
        ax.set_xlabel("$T_{\\rm wall}$ [kK]")
        ax.set_ylabel("$V_{\\rm top}$ [m/s]")
        ax.set_title("steady enthalpy-flux budget", fontsize=9, loc="left")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3, which="both")

    ax = axes[0, 1]
    ax.plot(tw, g("F_top_qs"), "o-", color="crimson", label="top cell")
    for hh, c in zip(heights, ["tab:blue", "tab:green", "tab:orange"]):
        ax.plot(tw, gh(hh, "F_qs"), "s--", color=c, ms=4, label=f"{hh:g} km")
    ax.axhline(0.0, color="0.4", lw=0.8, ls=":")
    ax.axvline(22.0, color="0.6", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1.0e-13)
    ax.set_xlabel("$T_{\\rm wall}$ [kK]")
    ax.set_ylabel(r"quasi-steady $\rho V$ [kg m$^{-2}$ s$^{-1}$]")
    ax.set_title("mass flux", fontsize=9, loc="left")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1, 0]
    if np.isfinite(g("q_wall_qs")).any():
        ax.plot(tw, g("q_wall_qs"), "o-", color="k")
        ax.axhline(0.0, color="0.4", lw=0.8, ls=":")
        ax.axvline(22.0, color="0.6", lw=0.8, ls="--")
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=1.0e-2)
        ax.set_xlabel("$T_{\\rm wall}$ [kK]")
        ax.set_ylabel("$q_{\\rm wall}$ [W m$^{-2}$]  (+ = into column)")
        ax.set_title("outer-face conductive flux", fontsize=9, loc="left")
        ax.grid(alpha=0.3, which="both")

    ax = axes[1, 1]
    ax.plot(tw, g("T_top_qs") / 1.0e3, "o-", color="tab:purple", label="$T_{\\rm top}$")
    ax.plot(tw, tw, ":", color="0.5", label="$T_{\\rm wall}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("$T_{\\rm wall}$ [kK]")
    ax.set_ylabel("quasi-steady $T$ [kK]")
    ax.set_title("top-cell temperature", fontsize=9, loc="left")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    fig.suptitle(f"model_column upper-boundary sensitivity: response vs "
                 f"$T_{{\\rm wall}}$ (last {avg_window:g} s of 4000 s)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def plot_support(runs, out_png, avg_window):
    t_walls = [r["meta"]["t_wall_K"] for r in runs]
    col, _, _ = wall_colors(t_walls)
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.6))
    for r in runs:
        d, m = r["data"], r["meta"]
        c = col[m["t_wall_K"]]
        lbl = f"{m['t_wall_K']/1e3:g} kK"
        axes[0].plot(d["h"], d["f_qs"], color=c, lw=1.3, label=lbl)
        axes[1].plot(d["h"], d["T_qs"] / 1.0e3, color=c, lw=1.3)
        axes[2].plot(d["times"], d["mass"] / d["mass"][0] - 1.0, color=c, lw=1.3)
    axes[0].set_yscale("symlog", linthresh=1.0e-12)
    axes[0].set_xlabel("h [km]")
    axes[0].set_ylabel(r"$\rho V$ [kg m$^{-2}$ s$^{-1}$]")
    axes[0].set_title("quasi-steady mass flux", fontsize=9, loc="left")
    axes[0].legend(ncol=2, fontsize=7)
    axes[1].set_xlabel("h [km]")
    axes[1].set_ylabel("T [kK]")
    axes[1].set_yscale("log")
    axes[1].set_title("quasi-steady temperature", fontsize=9, loc="left")
    axes[2].set_xlabel("t [s]")
    axes[2].set_ylabel("$M(t)/M(0) - 1$")
    axes[2].set_title("column mass drift", fontsize=9, loc="left")
    for ax in axes:
        ax.axhline(0.0, color="0.4", lw=0.7, ls=":")
        ax.grid(alpha=0.3)
    fig.suptitle("model_column upper-boundary sensitivity: supporting diagnostics",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="outputs/model_column/twall_sweep")
    ap.add_argument("--cfl", default="0.50",
                    help="only analyse runs with this CFL tag (keeps the "
                         "numerics uniform); use 'all' to take every run")
    ap.add_argument("--run", action="append", default=[],
                    help="explicit T_wall_K=path entry, added to the glob")
    ap.add_argument("--heights", type=float, nargs="*",
                    default=[1800.0, 2000.0, 2100.0])
    ap.add_argument("--avg-window", type=float, default=500.0)
    ap.add_argument("--out-dir", default="visualization/model_column")
    ap.add_argument("--prefix", default="twall")
    ap.add_argument("--json", default=None)
    ap.add_argument("--keep-incomplete", action="store_true",
                    help="also plot runs that did not reach CHROMO_T_END")
    args = ap.parse_args()

    entries = []
    for path in sorted(glob.glob(os.path.join(args.dir, "twall_*kK_cfl*.txt"))):
        mm = re.search(r"twall_(\d+)kK_cfl([0-9.]+)\.txt$", os.path.basename(path))
        if not mm:
            continue
        if args.cfl != "all" and mm.group(2) != args.cfl:
            continue
        entries.append((float(mm.group(1)) * 1.0e3, path, f"cfl{mm.group(2)}"))
    for spec in args.run:
        t_wall, path = spec.split("=", 1)
        entries.append((float(t_wall), path, "explicit"))
    if not entries:
        raise SystemExit(f"no runs matched in {args.dir} (cfl={args.cfl})")
    entries.sort()

    runs, metrics, rejected = [], [], []
    for t_wall, path, tag in entries:
        label = f"{t_wall/1e3:g}kK"
        meta, data = analyse(label, t_wall, path, args.heights, args.avg_window)
        meta["cfl_tag"] = tag
        # A run that did not reach CHROMO_T_END is not a 4000 s result and must
        # not sit on the same axes as one that did. Report it, then drop it.
        if not args.keep_incomplete and meta["termination"] != "end_time":
            rejected.append(meta)
            print(f"{label:>7s} {tag:>8s}  EXCLUDED: termination="
                  f"{meta['termination']} at t={meta['t_end']:.1f} s")
            continue
        runs.append({"meta": meta, "data": data})
        metrics.append(meta)
        print(f"{label:>7s} {tag:>8s}  term={meta['termination']} "
              f"fallbacks={meta['godunov_fallbacks']} "
              f"V_top={meta['V_top_qs']:+9.3f} m/s  "
              f"V_max={meta['V_max_qs']:+9.3f} @ {meta['h_V_max_qs']:.0f} km  "
              f"T_top={meta['T_top_qs']:8.1f} K  "
              f"q_wall={meta.get('q_wall_qs', float('nan')):+.3e} W/m2  "
              f"dM/M={meta['mass_rel_drift']:+.3e}  "
              f"|V|/Vcap={meta['vcap_frac_max']:.3f}  "
              f"V/V_pred={meta.get('V_ratio_pred', float('nan')):+6.2f}  "
              f"drift={(meta['V_top_drift_frac'] or float('nan')):+7.3f}")

    os.makedirs(args.out_dir, exist_ok=True)
    p = lambda name: os.path.join(args.out_dir, f"{args.prefix}_{name}")
    plot_v_of_t(runs, args.heights, p("V_of_t.png"), args.avg_window)
    plot_v_of_h(runs, p("V_of_h.png"), args.avg_window)
    plot_summary(metrics, args.heights, p("summary.png"), args.avg_window)
    plot_support(runs, p("support.png"), args.avg_window)
    print("wrote", p("V_of_t.png"), p("V_of_h.png"), p("summary.png"),
          p("support.png"), sep="\n  ")

    if args.json:
        os.makedirs(os.path.dirname(args.json) or ".", exist_ok=True)
        with open(args.json, "w") as fh:
            json.dump({"dir": args.dir, "cfl": args.cfl,
                       "avg_window_s": args.avg_window,
                       "heights_km": args.heights, "runs": metrics,
                       "excluded": rejected}, fh, indent=2)
        print("  " + args.json)


if __name__ == "__main__":
    main()
