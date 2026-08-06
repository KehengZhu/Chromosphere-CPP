#!/usr/bin/env python3
"""Compare CFL-sweep runs of model_column against a finer-time reference.

Parses the production outputs directly -- the conserved-variable snapshot file
written by chromo_main and its `.gamma_diag` sidecar (rho, v, T, x_eq, n_e,
n_HI, p_total, Gamma1, kappa) -- plus the stdout log for solver counters. None
of the solver equations are reimplemented here; every physical field either
comes straight out of a file column or is a product/ratio of file columns.

Usage
-----
    .venv/bin/python util/analyze_cfl_sweep.py \
        --reference outputs/model_column/cfl_sweep/cfl_reference_0125_20s.txt \
        --candidate 0.25=outputs/model_column/cfl_sweep/cfl_0250_20s.txt \
        --candidate 0.50=outputs/model_column/cfl_sweep/cfl_0500_20s.txt \
        --out outputs/model_column/cfl_sweep/cfl_accuracy_20s.md

Snapshots are matched by physical time (nearest snapshot within --time-tol).
Relative norms use a denominator floor of --floor-frac times the domain maximum
of |reference|, so a cell where the reference field is essentially zero cannot
manufacture a large relative error. The floor is reported for every field.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# gamma_diag column order, from chromo_main.cpp:
#   rho_total v_cm T x_eq n_e n_HI p_total Gamma1 kappa_physical kappa_solver
GAMMA_COLS = ["rho", "v", "T", "x_eq", "n_e", "n_HI", "p", "Gamma1",
              "kappa_phys", "kappa_solver"]

# Conserved-variable column order (chromosphere.hpp namespace cons).
CONS_COLS = ["RHO_I", "RHO_N", "MOM_I", "MOM_N", "E_I", "E_N", "E_E"]

# Transition-region criterion. model_column with ISO_T_TOP=22000 tops out in the
# upper chromosphere, so the usual 1e5 K coronal criterion never fires; the TRAC
# cutoff temperature the solver itself reports (T_c = 20000 K) is used instead.
# Overridable with --tr-temperature.
DEFAULT_TR_T = 20000.0


@dataclass
class Snapshot:
    time: float
    step: int
    data: np.ndarray          # (ns, ncol)


@dataclass
class Run:
    label: str
    path: Path
    height_km: np.ndarray     # (ns,) cumulative cell TOP edges
    ds_km: np.ndarray         # (ns,) cell widths
    cons: list                # list[Snapshot] over CONS_COLS
    gamma: list               # list[Snapshot] over GAMMA_COLS
    log: dict = field(default_factory=dict)

    def cell_center_km(self) -> np.ndarray:
        return self.height_km - 0.5 * self.ds_km


def _parse_snapshot_file(path: Path, ncol_expected: int):
    """Parse a chromo_main snapshot file (main output or .gamma_diag).

    Layout: 'ns ncol' / heights / optional '#' comment lines / repeated
    '# t = T step = S' markers each followed by ns rows of ncol floats.
    """
    with open(path) as fh:
        lines = fh.read().splitlines()
    if not lines:
        raise ValueError(f"{path}: empty")
    head = lines[0].split()
    ns, ncol = int(head[0]), int(head[1])
    if ncol != ncol_expected:
        raise ValueError(f"{path}: expected {ncol_expected} columns, header says {ncol}")
    height = np.asarray([float(x) for x in lines[1].split()], dtype=float)
    if height.size != ns:
        raise ValueError(f"{path}: height row has {height.size} entries, expected {ns}")

    marker = re.compile(r"^#\s*t\s*=\s*(\S+)\s+step\s*=\s*(\S+)")
    snaps = []
    i = 2
    n = len(lines)
    while i < n:
        line = lines[i]
        m = marker.match(line)
        if not m:
            i += 1           # header/comment lines between the grid and the data
            continue
        t, step = float(m.group(1)), int(m.group(2))
        rows = lines[i + 1: i + 1 + ns]
        if len(rows) < ns:
            break            # truncated trailing snapshot (interrupted run)
        block = np.fromstring(" ".join(rows), sep=" ")
        if block.size != ns * ncol:
            break
        snaps.append(Snapshot(t, step, block.reshape(ns, ncol)))
        i += 1 + ns
    if not snaps:
        raise ValueError(f"{path}: no complete snapshots found")
    return ns, height, snaps


LOG_PATTERNS = {
    "wall_s":               r"^real\s+([0-9.]+)",
    "final_step":           r"^final_step=(\d+)",
    "final_time":           r"^final_physical_time=([0-9.eE+-]+)",
    "requested_end_time":   r"^requested_end_time=([0-9.eE+-]+)",
    "termination":          r"^termination=(\w+)",
    "effective_step_cap":   r"^effective_step_cap=(\d+)",
    "time.rhs_s":           r"^time\.rhs_s=([0-9.eE+-]+)",
    "time.conduction_s":    r"^time\.conduction_s=([0-9.eE+-]+)",
    "time.decode_s":        r"^time\.decode_s=([0-9.eE+-]+)",
    "time.projection_s":    r"^time\.projection_s=([0-9.eE+-]+)",
    "time.cfl_s":           r"^time\.cfl_s=([0-9.eE+-]+)",
    "time.boundary_s":      r"^time\.boundary_s=([0-9.eE+-]+)",
    "time.output_s":        r"^time\.output_s=([0-9.eE+-]+)",
    "limiter.acoustic":     r"^limiter\.acoustic=(\d+)",
    "limiter.other":        r"^limiter\.other=(\d+)",
    "limiter.cooling":      r"^limiter\.cooling=(\d+)",
    "limiter.electron_energy": r"^limiter\.electron_energy=(\d+)",
    "limiter.beam_heating": r"^limiter\.beam_heating=(\d+)",
    "limiter.coronal_heating": r"^limiter\.coronal_heating=(\d+)",
    "inversion.calls":      r"^inversion\.calls=(\d+)",
    "inversion.initial_guess_accepts": r"^inversion\.initial_guess_accepts=(\d+)",
    "inversion.one_update_convergences": r"^inversion\.one_update_convergences=(\d+)",
    "inversion.average_iterations": r"^inversion\.average_iterations=([0-9.eE+-]+)",
    "inversion.maximum_iterations": r"^inversion\.maximum_iterations=(\d+)",
    "inversion.bisection_fallbacks": r"^inversion\.bisection_fallbacks=(\d+)",
    "conduction.calls":     r"^conduction\.calls=(\d+)",
    "conduction.average_newton_updates": r"^conduction\.average_newton_updates=([0-9.eE+-]+)",
    "conduction.maximum_newton_updates": r"^conduction\.maximum_newton_updates=(\d+)",
}


def parse_log(path: Path) -> dict:
    out = {}
    if not path.exists():
        return out
    text = path.read_text(errors="replace")
    for key, pat in LOG_PATTERNS.items():
        m = re.search(pat, text, re.MULTILINE)
        if m:
            val = m.group(1)
            try:
                out[key] = int(val) if re.fullmatch(r"\d+", val) else float(val)
            except ValueError:
                out[key] = val
    m = re.search(r"^termination=(\w+)", text, re.MULTILINE)
    if m:
        out["termination"] = m.group(1)
    # Mean/min dt from the periodic progress lines (a cheap, log-only estimate).
    dts = [float(x) for x in re.findall(r"dt = ([0-9.eE+-]+)", text)]
    if dts:
        out["dt_progress_mean"] = float(np.mean(dts))
        out["dt_progress_min"] = float(np.min(dts))
    return out


def load_run(label: str, path: Path) -> Run:
    ns, height, cons = _parse_snapshot_file(path, len(CONS_COLS))
    gpath = Path(str(path) + ".gamma_diag")
    if not gpath.exists():
        raise FileNotFoundError(f"missing gamma diagnostic sidecar: {gpath}")
    ns_g, height_g, gamma = _parse_snapshot_file(gpath, len(GAMMA_COLS))
    if ns_g != ns or not np.allclose(height_g, height):
        raise ValueError(f"{path}: snapshot and gamma_diag grids disagree")
    # Cell widths from the cumulative top-edge column. The first cell's width is
    # not recoverable from the diffs alone (the base offset out_base_km is not in
    # the file); the inner region of this mesh is uniform and only the OUTER end
    # is refined, so ds[0] == ds[1] exactly for these runs.
    ds = np.empty_like(height)
    ds[1:] = np.diff(height)
    ds[0] = ds[1]
    log = parse_log(path.with_suffix(path.suffix + ".log"))
    if not log:
        log = parse_log(Path(str(path).rsplit(".", 1)[0] + ".log"))
    return Run(label, path, height, ds, cons, gamma, log)


def match_snapshots(ref: list, cand: list, tol: float):
    """Pair reference and candidate snapshots by nearest physical time."""
    ref_t = np.array([s.time for s in ref])
    pairs = []
    for cs in cand:
        j = int(np.argmin(np.abs(ref_t - cs.time)))
        if abs(ref_t[j] - cs.time) <= tol:
            pairs.append((ref[j], cs))
    return pairs


def field_errors(ref: np.ndarray, cand: np.ndarray, ds: np.ndarray,
                 height: np.ndarray, floor_frac: float) -> dict:
    """Mesh-weighted relative L1/L2/Linf plus absolute Linf and its location."""
    diff = cand - ref
    scale = float(np.max(np.abs(ref)))
    floor = floor_frac * scale if scale > 0 else 1.0
    denom_cell = np.maximum(np.abs(ref), floor)

    l1 = float(np.sum(np.abs(diff) * ds) / np.sum(denom_cell * ds))
    l2 = float(math.sqrt(np.sum(diff ** 2 * ds) / np.sum(denom_cell ** 2 * ds)))
    rel_cell = np.abs(diff) / denom_cell
    k = int(np.argmax(rel_cell))
    ka = int(np.argmax(np.abs(diff)))
    # Cells whose relative error is set by the floor rather than by the local
    # reference value are flagged separately (see the module docstring).
    floored = np.abs(ref) < floor
    rel_unfloored = np.where(floored, 0.0, rel_cell)
    ku = int(np.argmax(rel_unfloored))
    return {
        "rel_L1": l1,
        "rel_L2": l2,
        "rel_Linf": float(rel_cell[k]),
        "rel_Linf_km": float(height[k]),
        "rel_Linf_floored": bool(floored[k]),
        "rel_Linf_unfloored": float(rel_unfloored[ku]),
        "rel_Linf_unfloored_km": float(height[ku]),
        "abs_Linf": float(np.abs(diff)[ka]),
        "abs_Linf_km": float(height[ka]),
        "floor": floor,
        "ref_max_abs": scale,
    }


def observables(snap_g: np.ndarray, snap_c: np.ndarray, ds_km: np.ndarray,
                tr_T: float, height: np.ndarray) -> dict:
    rho = snap_g[:, GAMMA_COLS.index("rho")]
    v = snap_g[:, GAMMA_COLS.index("v")]
    T = snap_g[:, GAMMA_COLS.index("T")]
    p = snap_g[:, GAMMA_COLS.index("p")]
    ds_m = ds_km * 1.0e3
    mass_flux = rho * v
    e_tot = snap_c[:, CONS_COLS.index("E_I")] + snap_c[:, CONS_COLS.index("E_N")]

    # TR location: first upward crossing of tr_T, linearly interpolated in height.
    tr_km = float("nan")
    above = np.nonzero(T >= tr_T)[0]
    if above.size and above[0] > 0:
        i = above[0]
        t0, t1 = T[i - 1], T[i]
        if t1 != t0:
            f = (tr_T - t0) / (t1 - t0)
            tr_km = float(height[i - 1] + f * (height[i] - height[i - 1]))
        else:
            tr_km = float(height[i])
    elif above.size:
        tr_km = float(height[0])

    grad = np.abs(np.gradient(T, height))
    return {
        "integrated_mass": float(np.sum(rho * ds_m)),
        "integrated_energy": float(np.sum(e_tot * ds_m)),
        "max_up_mass_flux": float(np.max(mass_flux)),
        "max_down_mass_flux": float(np.min(mass_flux)),
        "T_max": float(np.max(T)),
        "T_min": float(np.min(T)),
        "p_top": float(p[-1]),
        "v_top": float(v[-1]),
        "tr_km": tr_km,
        "max_dTds_km": float(height[int(np.argmax(grad))]),
        "max_dTds": float(np.max(grad)),
    }


def rel(a: float, b: float) -> float:
    """Relative difference of a from reference b, guarding a zero reference."""
    if b == 0.0:
        return float("nan") if a != 0.0 else 0.0
    return (a - b) / abs(b)


def fmt(x, nd=4):
    if isinstance(x, str):
        return x
    if x is None:
        return "-"
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return "nan"
    if isinstance(x, int):
        return str(x)
    return f"{x:.{nd}g}"


def table(rows, headers) -> str:
    cols = len(headers)
    widths = [len(h) for h in headers]
    srows = []
    for r in rows:
        s = [fmt(v) for v in r]
        srows.append(s)
        for i in range(cols):
            widths[i] = max(widths[i], len(s[i]))
    out = ["| " + " | ".join(h.ljust(widths[i]) for i, h in enumerate(headers)) + " |",
           "|" + "|".join("-" * (widths[i] + 2) for i in range(cols)) + "|"]
    for s in srows:
        out.append("| " + " | ".join(s[i].ljust(widths[i]) for i in range(cols)) + " |")
    return "\n".join(out)


# Provisional engineering acceptance thresholds (relative, fractional).
THRESHOLDS = {
    "integrated_mass": 0.005,
    "integrated_energy": 0.005,
    "max_up_mass_flux": 0.01,
    "T_max": 0.01,
    "p_top": 0.01,
    "v_top": 0.01,
    "field_rel_L1": 0.01,
    "field_rel_Linf": 0.03,
}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--reference-label", default="ref")
    ap.add_argument("--candidate", action="append", default=[],
                    metavar="LABEL=PATH",
                    help="candidate run; repeatable")
    ap.add_argument("--floor-frac", type=float, default=1.0e-3,
                    help="relative-error denominator floor as a fraction of the "
                         "domain max |reference| (default 1e-3)")
    ap.add_argument("--time-tol", type=float, default=0.05,
                    help="snapshot time-matching tolerance in seconds")
    ap.add_argument("--tr-temperature", type=float, default=DEFAULT_TR_T)
    ap.add_argument("--fields", default="rho,v,T,p,x_eq",
                    help="gamma_diag fields to score")
    ap.add_argument("--out", type=Path, help="write the markdown report here")
    args = ap.parse_args(argv)

    ref = load_run(args.reference_label, args.reference)
    cands = []
    for spec in args.candidate:
        if "=" not in spec:
            ap.error(f"--candidate expects LABEL=PATH, got {spec!r}")
        label, path = spec.split("=", 1)
        cands.append(load_run(label, Path(path)))

    fields = [f.strip() for f in args.fields.split(",") if f.strip()]
    for f in fields:
        if f not in GAMMA_COLS + ["mass_flux", "e_tot"]:
            ap.error(f"unknown field {f!r}")

    out = []
    out.append(f"# CFL sweep accuracy report")
    out.append("")
    out.append(f"Reference: `{args.reference}` (label `{ref.label}`), "
               f"{len(ref.gamma)} snapshots, {ref.height_km.size} cells.")
    out.append(f"Relative-error denominator floor: {args.floor_frac:g} x "
               f"domain max |reference| (per field, per snapshot).")
    out.append(f"TR criterion: first upward crossing of T = {args.tr_temperature:g} K, "
               f"linearly interpolated in height.")
    out.append("")

    verdicts = {}

    # ---- solver behaviour ---------------------------------------------------
    out.append("## Solver counters and performance")
    out.append("")
    keys = ["wall_s", "final_step", "final_time", "termination",
            "dt_progress_mean", "dt_progress_min",
            "limiter.acoustic", "limiter.other",
            "inversion.calls", "inversion.initial_guess_accepts",
            "inversion.one_update_convergences", "inversion.average_iterations",
            "inversion.maximum_iterations", "inversion.bisection_fallbacks",
            "conduction.average_newton_updates", "conduction.maximum_newton_updates"]
    rows = []
    for run in [ref] + cands:
        lg = run.log
        row = [run.label]
        for k in keys:
            row.append(lg.get(k, "-"))
        # derived rates
        calls = lg.get("inversion.calls") or 0
        row.append(lg["inversion.initial_guess_accepts"] / calls if calls else "-")
        row.append(lg["inversion.one_update_convergences"] / calls if calls else "-")
        row.append(lg["inversion.bisection_fallbacks"] / calls if calls else "-")
        rows.append(row)
    out.append(table(rows, ["run"] + keys +
                     ["guess_accept_rate", "one_update_rate", "fallback_rate"]))
    out.append("")

    # ---- per-candidate accuracy --------------------------------------------
    for cand in cands:
        out.append(f"## Candidate `{cand.label}` vs `{ref.label}`")
        out.append("")
        pairs_g = match_snapshots(ref.gamma, cand.gamma, args.time_tol)
        pairs_c = match_snapshots(ref.cons, cand.cons, args.time_tol)
        if not pairs_g:
            out.append("**No snapshots matched within the time tolerance.**")
            verdicts[cand.label] = "NO DATA"
            continue
        if cand.height_km.size != ref.height_km.size or \
           not np.allclose(cand.height_km, ref.height_km):
            out.append("**Grids differ between reference and candidate — "
                       "field comparison aborted.**")
            verdicts[cand.label] = "GRID MISMATCH"
            continue

        # field errors, worst over all matched times
        rows = []
        worst_l1, worst_linf = 0.0, 0.0
        for f in fields:
            worst = None
            for rs, cs in pairs_g:
                if f == "mass_flux":
                    r = rs.data[:, GAMMA_COLS.index("rho")] * rs.data[:, GAMMA_COLS.index("v")]
                    c = cs.data[:, GAMMA_COLS.index("rho")] * cs.data[:, GAMMA_COLS.index("v")]
                else:
                    j = GAMMA_COLS.index(f)
                    r, c = rs.data[:, j], cs.data[:, j]
                e = field_errors(r, c, cand.ds_km, cand.height_km, args.floor_frac)
                e["t"] = cs.time
                if worst is None or e["rel_L1"] > worst["rel_L1"]:
                    worst = e
            rows.append([f, worst["t"], worst["rel_L1"], worst["rel_L2"],
                         worst["rel_Linf"], worst["rel_Linf_km"],
                         worst["rel_Linf_unfloored"], worst["abs_Linf"],
                         worst["abs_Linf_km"], worst["floor"]])
            worst_l1 = max(worst_l1, worst["rel_L1"])
            worst_linf = max(worst_linf, worst["rel_Linf_unfloored"])
        # total conserved energy from the main snapshot file
        worst = None
        for rs, cs in pairs_c:
            r = rs.data[:, CONS_COLS.index("E_I")] + rs.data[:, CONS_COLS.index("E_N")]
            c = cs.data[:, CONS_COLS.index("E_I")] + cs.data[:, CONS_COLS.index("E_N")]
            e = field_errors(r, c, cand.ds_km, cand.height_km, args.floor_frac)
            e["t"] = cs.time
            if worst is None or e["rel_L1"] > worst["rel_L1"]:
                worst = e
        if worst is not None:
            rows.append(["e_tot(E_I+E_N)", worst["t"], worst["rel_L1"], worst["rel_L2"],
                         worst["rel_Linf"], worst["rel_Linf_km"],
                         worst["rel_Linf_unfloored"], worst["abs_Linf"],
                         worst["abs_Linf_km"], worst["floor"]])
            worst_l1 = max(worst_l1, worst["rel_L1"])
            worst_linf = max(worst_linf, worst["rel_Linf_unfloored"])

        out.append("### Field errors (worst matched snapshot per field)")
        out.append("")
        out.append(table(rows, ["field", "t_worst", "rel_L1", "rel_L2", "rel_Linf",
                                "Linf_km", "rel_Linf_nofloor", "abs_Linf",
                                "abs_Linf_km", "floor"]))
        out.append("")

        # observables at every matched time
        obs_rows = []
        worst_obs = {k: 0.0 for k in
                     ["integrated_mass", "integrated_energy", "max_up_mass_flux",
                      "T_max", "p_top", "v_top"]}
        worst_tr = 0.0
        cons_by_time = {round(cs.time, 6): cs for _, cs in pairs_c}
        ref_cons_by_time = {round(rs.time, 6): rs for rs, _ in pairs_c}
        for (rs, cs), (rc, cc) in zip(pairs_g, pairs_c):
            ro = observables(rs.data, rc.data, ref.ds_km, args.tr_temperature, ref.height_km)
            co = observables(cs.data, cc.data, cand.ds_km, args.tr_temperature, cand.height_km)
            row = [cs.time]
            for k in ["integrated_mass", "integrated_energy", "max_up_mass_flux",
                      "T_max", "T_min", "p_top", "v_top"]:
                d = rel(co[k], ro[k])
                row.append(d)
                if k in worst_obs:
                    worst_obs[k] = max(worst_obs[k], abs(d))
            dtr = co["tr_km"] - ro["tr_km"]
            row.append(dtr)
            row.append(co["max_dTds_km"] - ro["max_dTds_km"])
            obs_rows.append(row)
            if not math.isnan(dtr):
                worst_tr = max(worst_tr, abs(dtr))
        out.append("### Physical observables (relative difference vs reference; "
                   "TR/gradient columns in km)")
        out.append("")
        out.append(table(obs_rows, ["t", "d_mass", "d_energy", "d_max_up_flux",
                                    "d_T_max", "d_T_min", "d_p_top", "d_v_top",
                                    "d_TR_km", "d_gradmax_km"]))
        out.append("")

        min_ds_km = float(np.min(cand.ds_km))
        checks = []
        for k, lim in [("integrated_mass", THRESHOLDS["integrated_mass"]),
                       ("integrated_energy", THRESHOLDS["integrated_energy"]),
                       ("max_up_mass_flux", THRESHOLDS["max_up_mass_flux"]),
                       ("T_max", THRESHOLDS["T_max"]),
                       ("p_top", THRESHOLDS["p_top"]),
                       ("v_top", THRESHOLDS["v_top"])]:
            checks.append([k, worst_obs[k], lim, "PASS" if worst_obs[k] <= lim else "FAIL"])
        checks.append(["field_rel_L1 (worst)", worst_l1, THRESHOLDS["field_rel_L1"],
                       "PASS" if worst_l1 <= THRESHOLDS["field_rel_L1"] else "FAIL"])
        checks.append(["field_rel_Linf (worst, unfloored)", worst_linf,
                       THRESHOLDS["field_rel_Linf"],
                       "PASS" if worst_linf <= THRESHOLDS["field_rel_Linf"] else "FAIL"])
        checks.append(["TR shift (km)", worst_tr, min_ds_km,
                       "PASS" if worst_tr <= min_ds_km else "FAIL"])
        out.append("### Acceptance gate (worst over all matched snapshots)")
        out.append("")
        out.append(table(checks, ["metric", "worst", "threshold", "status"]))
        out.append("")
        verdict = "PASS" if all(c[3] == "PASS" for c in checks) else "FAIL"
        verdicts[cand.label] = verdict
        out.append(f"**Candidate `{cand.label}`: {verdict}**")
        out.append("")

    out.append("## Summary")
    out.append("")
    out.append(table([[k, v] for k, v in verdicts.items()], ["candidate", "verdict"]))
    out.append("")

    text = "\n".join(out)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
        print(f"wrote {args.out}")
    print(text)
    return 0 if all(v == "PASS" for v in verdicts.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
