#!/usr/bin/env python3
"""Reconstruction-comparison diagnostic for the model_column Gamma/Saha path.

Answers one question: does reconstructing (ln rho, V, ln p) instead of
(ln rho, V, ln T) reduce the artificial face-pressure mismatch at the Saha
ionization transition, and does that carry through to the late-time velocity /
rho*V grid mode and the bulk evaporation/conduction observables?

Metric definitions (Q_V, Q_M, F_eff, Q_eff, boundary layer, mass drift) are the
same ones util/focused_long_resolution_diag.py uses, so numbers are directly
comparable with the existing N500/N1000 release evidence. The face-pressure
mismatch is computed from the SAME estimator for every run --- the authoritative
Saha pressure p(rho,T) evaluated at the captured one-sided reconstructed face
states --- so runs written before the .faceflux gained explicit p_L/p_R columns
are handled identically to new ones.

Usage:
  util/pressure_reconstruction_diag.py \
      --run label=outputs/model_column/base.txt [...] \
      [--targets 100,500,1000] [--json out.json]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

# Roughness/window conventions shared with the existing release diagnostics.
WINDOW = (2130.0, 2150.0)          # velocity / rho*V roughness window [km]
P_WINDOW = (2139.5, 2143.5)        # Saha-transition pressure-mismatch window [km]

# Legacy 26-column layout; runs that also carry the reconstructed pressures
# append p_cell p_L p_R.
FF_BASE = ("cell_km face_km rho_cell v_cell T_cell rho_L rho_R v_L v_R T_L T_R "
           "cs_L cs_R a_face f_central f_diff f_total eq_residual_mass "
           "r_rho phi_plus_rho r_ip1_rho phi_minus_rho r_v phi_plus_v "
           "r_T phi_plus_T").split()
FF_EXTRA = ["p_cell", "p_L", "p_R"]

K_B = 1.380649e-23
M_E = 9.1093837015e-31
M_H = 1.6726219e-27
H_PLANCK = 6.62607015e-34
CHI_H = 2.179872361e-18
SAHA_LOG_C = 1.5 * np.log(2.0 * np.pi * M_E * K_B / (H_PLANCK * H_PLANCK))


def saha_x(n_h: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Pure-H Saha fraction, same branch split as src/eos.cpp."""
    log_a = SAHA_LOG_C + 1.5 * np.log(T) - CHI_H / (K_B * T) - np.log(n_h)
    out = np.empty_like(log_a)
    hi = log_a >= 0.0
    out[hi] = 2.0 / (1.0 + np.sqrt(1.0 + 4.0 * np.exp(-log_a[hi])))
    u = np.exp(0.5 * log_a[~hi])
    out[~hi] = 2.0 * u / (u + np.sqrt(u * u + 4.0))
    return out


def saha_pressure(rho: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Authoritative total pressure p_i + p_n = (1+x) n_H k_B T."""
    n_h = np.asarray(rho, dtype=float) / M_H
    T = np.asarray(T, dtype=float)
    return (1.0 + saha_x(n_h, T)) * n_h * K_B * T


def roughness(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.mean(np.abs(np.diff(x, 2))) / max(np.mean(np.abs(x)), 1e-300))


def relstd(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.std(x) / max(abs(np.mean(x)), 1e-300))


def outer_boundary_layer(h_km: np.ndarray, T: np.ndarray, frac: float = 0.1) -> dict:
    s = np.asarray(h_km) * 1e3
    G = np.abs(np.gradient(np.asarray(T), s))
    gmax = float(G.max())
    i = len(G) - 1
    while i > 0 and G[i - 1] >= frac * gmax:
        i -= 1
    return {"bl_thickness_km": float(h_km[-1] - h_km[i]),
            "bl_n_cells": int(len(G) - i),
            "h_Gmax_km": float(h_km[int(np.argmax(G))])}


class Nearest:
    """Keep, per target time, the record whose time is closest to it."""

    def __init__(self, targets):
        self.targets = list(targets)
        self.best = {}

    def offer(self, t, value):
        for target in self.targets:
            d = abs(t - target)
            if target not in self.best or d < self.best[target][0]:
                self.best[target] = (d, t, value)

    def result(self):
        return {t: self.best[t][2] for t in self.targets if t in self.best}


def stream_gamma(base: Path, targets) -> dict:
    path = Path(str(base) + ".gamma_diag")
    if not path.exists():
        return {}
    near = Nearest(targets)
    with path.open() as fh:
        ns, ncol = map(int, fh.readline().split())
        h = np.fromstring(fh.readline(), sep=" ")
        cur_t, rows = None, []

        def flush():
            nonlocal rows
            if cur_t is None or len(rows) != ns:
                rows = []
                return
            arr = np.asarray(rows, dtype=float)
            m = (h >= WINDOW[0]) & (h <= WINDOW[1])
            v = arr[m, 1]
            metrics = {"t_gamma": float(cur_t),
                       "mean_V_ms": float(v.mean()),
                       "min_V_ms": float(v.min()),
                       "max_V_ms": float(v.max()),
                       "Q_V": roughness(v),
                       "n_negative_V": int(np.count_nonzero(v < 0)),
                       "T_top_K": float(arr[-1, 2]),
                       "p_top_Pa": float(arr[-1, 6]),
                       "max_abs_V_global_ms": float(np.max(np.abs(arr[:, 1])))}
            metrics.update(outer_boundary_layer(h, arr[:, 2]))
            near.offer(cur_t, metrics)
            rows = []

        for line in fh:
            if line.startswith("# t ="):
                flush()
                cur_t = float(line.split("t =")[1].split("step")[0])
            elif line.startswith("#"):
                continue
            elif cur_t is not None:
                vals = np.fromstring(line, sep=" ")
                if vals.size == ncol:
                    rows.append(vals)
        flush()
    return near.result()


def stream_faceflux(base: Path, targets) -> tuple[dict, dict]:
    path = Path(str(base) + ".faceflux")
    if not path.exists():
        return {}, {}
    near = Nearest(targets)
    meta: dict = {}
    ref = None
    mass0 = None
    cur_t = cur_step = None
    rows: list = []
    ncol = None

    def flush():
        nonlocal rows, ref, mass0
        if cur_t is None or not rows:
            rows = []
            return
        arr = np.asarray(rows, dtype=float)
        names = FF_BASE + (FF_EXTRA if arr.shape[1] == len(FF_BASE) + 3 else [])
        if arr.shape[1] != len(names):
            rows = []
            return
        c = {n: arr[:, k] for k, n in enumerate(names)}
        widths_km = np.concatenate(([2.0 * (c["face_km"][0] - c["cell_km"][0])],
                                    np.diff(c["face_km"])))
        mass = float(np.sum(c["rho_cell"] * widths_km * 1e3))
        if ref is None:
            ref = c["f_total"].copy()
            mass0 = mass
        if len(ref) != len(c["f_total"]):
            rows = []
            return

        mc = (c["cell_km"] >= WINDOW[0]) & (c["cell_km"] <= WINDOW[1])
        mf = (c["face_km"] >= WINDOW[0]) & (c["face_km"] <= WINDOW[1])
        M = c["rho_cell"] * c["v_cell"]
        f_eff = c["f_total"] - ref

        # Reconstructed face-pressure mismatch through the Saha transition. One
        # estimator for every run: the authoritative p(rho,T) at the captured
        # one-sided corrector states.
        pw = (c["face_km"] >= P_WINDOW[0]) & (c["face_km"] <= P_WINDOW[1])
        pl = saha_pressure(c["rho_L"][pw], c["T_L"][pw])
        pr = saha_pressure(c["rho_R"][pw], c["T_R"][pw])
        eps = np.abs(pr - pl) / (0.5 * (pr + pl))
        imax = int(np.argmax(eps)) if eps.size else 0
        metrics = {
            "t_faceflux": float(cur_t), "step_faceflux": int(cur_step),
            "M_mean": float(M[mc].mean()), "Q_M": roughness(M[mc]),
            "F_eff_mean": float(f_eff[mf].mean()), "Q_eff": roughness(f_eff[mf]),
            "column_mass_kgm2": mass,
            "column_mass_rel_change": float(mass / mass0 - 1.0),
            "eps_p_median": float(np.median(eps)) if eps.size else float("nan"),
            "eps_p_max": float(eps.max()) if eps.size else float("nan"),
            "eps_p_max_km": float(c["face_km"][pw][imax]) if eps.size else float("nan"),
            "eps_p_n_faces": int(pw.sum()),
        }
        if eps.size:
            rho_at = c["rho_L"][pw][imax]
            t_at = c["T_L"][pw][imax]
            metrics.update({
                "eps_p_max_T_K": float(t_at),
                "eps_p_max_x_eq": float(saha_x(np.array([rho_at / M_H]),
                                               np.array([t_at]))[0]),
                "eps_p_max_v_L_ms": float(c["v_L"][pw][imax]),
            })
        if arr.shape[1] == len(FF_BASE) + 3:   # direct check of the capture
            pw_l, pw_r = c["p_L"][pw], c["p_R"][pw]
            eps_direct = np.abs(pw_r - pw_l) / (0.5 * (pw_r + pw_l))
            metrics["eps_p_max_direct"] = float(eps_direct.max()) if eps.size else 0.0
        near.offer(cur_t, metrics)
        rows = []

    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("# t "):
                    flush()
                    parts = line.split()
                    cur_t, cur_step = float(parts[3]), int(parts[6])
                else:
                    for tok in line[1:].split():
                        if "=" in tok and not tok.startswith("columns"):
                            k, v = tok.split("=", 1)
                            meta.setdefault(k, v)
                continue
            vals = np.fromstring(line, sep=" ")
            if ncol is None:
                ncol = vals.size
            if vals.size == ncol:
                rows.append(vals)
        flush()
    return meta, near.result()


def stream_outercond(base: Path, targets) -> dict:
    """Interpolate the outer conductive flux and integrate it to each target."""
    path = Path(str(base) + ".outercond")
    if not path.exists():
        return {}
    out: dict = {}
    prev = None
    last_row = None
    integral = 0.0
    remaining = sorted(targets)
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            v = np.fromstring(line, sep=" ")
            if v.size < 12:
                continue
            t, q_phys, q_total = float(v[0]), float(v[9]), float(v[11])
            if prev is not None:
                t0, q0 = prev[0], prev[1]
                while remaining and remaining[0] <= t:
                    target = remaining.pop(0)
                    frac = 0.0 if t == t0 else (target - t0) / (t - t0)
                    q_at = q0 + frac * (q_phys - q0)
                    out[target] = {
                        "q_phys_W_m2": q_at,
                        "cumulative_phys_J_m2": integral
                        + 0.5 * (q0 + q_at) * (target - t0),
                        "T_top_bc_K": float(v[2]), "q_total_W_m2": q_total}
                integral += 0.5 * (q0 + q_phys) * (t - t0)
            prev = (t, q_phys)
            last_row = v
    # A target that falls between the last captured record and the run's end time
    # (the capture stride rarely lands exactly on it) is reported at that record.
    if remaining and prev is not None:
        for target in list(remaining):
            if target - prev[0] <= 0.02 * max(target, 1.0):
                out[target] = {"q_phys_W_m2": prev[1],
                               "cumulative_phys_J_m2": integral,
                               "T_top_bc_K": float(last_row[2]),
                               "q_total_W_m2": float(last_row[11]),
                               "cumulative_at_t_s": prev[0]}
    return out


def analyze(label: str, base: str, targets) -> dict:
    p = Path(base)
    gamma = stream_gamma(p, targets)
    meta, ff = stream_faceflux(p, targets)
    cond = stream_outercond(p, targets)
    per_time = {}
    for t in targets:
        row = {}
        row.update(gamma.get(t, {}))
        row.update(ff.get(t, {}))
        row.update(cond.get(t, {}))
        if row:
            per_time[t] = row
    return {"label": label, "base": base, "faceflux_meta": meta, "times": per_time}


ROWS = [("t_gamma", "t [s]", "{:.1f}"),
        ("eps_p_max", "eps_p max", "{:.4%}"),
        ("eps_p_median", "eps_p med", "{:.4%}"),
        ("eps_p_max_km", "eps_p max @ [km]", "{:.2f}"),
        ("eps_p_max_x_eq", "  x_eq there", "{:.4f}"),
        ("Q_V", "Q_V", "{:.5f}"),
        ("Q_M", "Q_M", "{:.5f}"),
        ("Q_eff", "Q_eff", "{:.2e}"),
        ("mean_V_ms", "mean V [m/s]", "{:.3f}"),
        ("min_V_ms", "min V [m/s]", "{:.3f}"),
        ("max_V_ms", "max V [m/s]", "{:.3f}"),
        ("n_negative_V", "n(V<0)", "{:d}"),
        ("F_eff_mean", "F_eff [kg/m2/s]", "{:.4e}"),
        ("T_top_K", "T_top [K]", "{:.1f}"),
        ("p_top_Pa", "p_top [Pa]", "{:.6g}"),
        ("bl_thickness_km", "BL thick [km]", "{:.3f}"),
        ("q_phys_W_m2", "q_phys [W/m2]", "{:.4f}"),
        ("cumulative_phys_J_m2", "E_phys [J/m2]", "{:.2f}"),
        ("column_mass_rel_change", "mass drift", "{:.3e}")]


def report(results, targets) -> str:
    lines = []
    for t in targets:
        present = [r for r in results if t in r["times"]]
        if not present:
            continue
        lines.append(f"\n=== target t = {t:g} s ===")
        width = max(len(r["label"]) for r in present) + 2
        header = " " * 20 + "".join(r["label"].ljust(width) for r in present)
        lines.append(header)
        for key, name, fmt in ROWS:
            cells = []
            for r in present:
                v = r["times"][t].get(key)
                cells.append(("n/a" if v is None else fmt.format(v)).ljust(width))
            if all(c.strip() == "n/a" for c in cells):
                continue
            lines.append(name.ljust(20) + "".join(cells))
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True,
                    metavar="LABEL=PATH", help="run to analyze (repeatable)")
    ap.add_argument("--targets", default="100,500,1000,2000,4000")
    ap.add_argument("--json", default=None)
    args = ap.parse_args()
    targets = [float(x) for x in args.targets.split(",")]

    results = []
    for spec in args.run:
        label, _, path = spec.partition("=")
        results.append(analyze(label, path, targets))
    print(report(results, targets))
    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=2))
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
