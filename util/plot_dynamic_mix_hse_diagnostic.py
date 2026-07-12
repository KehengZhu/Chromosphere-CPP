#!/usr/bin/env python3
"""Analyze total-mixture dynamic-HSE diagnostic CSVs at a matched time.

The input is sampled every N Stage-D calls. This script streams each large CSV
twice, selects the two samples bracketing ``--time``, and linearly interpolates
the four lower-chromosphere diagnostics to that physical time. It never treats
the charged/neutral pressure partitions as independently hydrostatic.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/chromosphere-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/chromosphere-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_FILES = [
    "outputs/iso_dynamic_mix_hse_ns2000.csv",
    "outputs/iso_dynamic_mix_hse_ns3000.csv",
    "outputs/iso_dynamic_mix_hse_ns4000.csv",
]

METRICS = (
    ("dp_cond_mix_over_dt", r"$\Delta p_{\rm cond,mix}/\Delta t$ [Pa s$^{-1}$]"),
    ("dpi_cond_mix_over_dt", r"$\Delta \pi_{\rm cond,mix}/\Delta t$ [Pa s$^{-1}$]"),
    ("pi_mix_after", r"$\pi_{\rm mix}$ after conduction [Pa]"),
    ("mass_flux_mix", r"$\rho_iV+\rho_nU$ [kg m$^{-2}$ s$^{-1}$]"),
)


def csv_rows(path: Path):
    with path.open(newline="") as stream:
        yield from csv.DictReader(line for line in stream if not line.startswith("#"))


def sample_times(path: Path) -> list[tuple[float, int]]:
    samples: list[tuple[float, int]] = []
    previous_stage: int | None = None
    for row in csv_rows(path):
        stage = int(row["stage"])
        if stage != previous_stage:
            samples.append((float(row["time_s"]), stage))
            previous_stage = stage
    if not samples:
        raise ValueError(f"no data rows in {path}")
    return samples


def time_bracket(samples: list[tuple[float, int]], target_time: float) -> tuple[tuple[float, int], tuple[float, int]]:
    samples = sorted(samples)
    if target_time < samples[0][0] or target_time > samples[-1][0]:
        raise ValueError(
            f"t={target_time:g} s outside [{samples[0][0]:.6g}, {samples[-1][0]:.6g}] s"
        )
    for left, right in zip(samples, samples[1:]):
        if left[0] == target_time:
            return left, left
        if left[0] <= target_time <= right[0]:
            return left, right
    return samples[-1], samples[-1]


def read_profiles(path: Path, stages: set[int]) -> dict[int, dict[str, np.ndarray]]:
    names = ["height_km"] + [metric for metric, _ in METRICS]
    columns = {stage: {name: [] for name in names} for stage in stages}
    for row in csv_rows(path):
        stage = int(row["stage"])
        if stage not in columns:
            continue
        columns[stage]["height_km"].append(float(row["height_km"]))
        for metric, _ in METRICS:
            columns[stage][metric].append(float(row[metric]))

    profiles: dict[int, dict[str, np.ndarray]] = {}
    for stage, values in columns.items():
        if not values["height_km"]:
            raise ValueError(f"selected Stage-D sample {stage} missing from {path}")
        profiles[stage] = {name: np.asarray(value, dtype=float) for name, value in values.items()}
    return profiles


def interpolate_profile(left: dict[str, np.ndarray], right: dict[str, np.ndarray], weight: float) -> dict[str, np.ndarray]:
    h_left = left["height_km"]
    if not np.allclose(h_left, right["height_km"], rtol=0.0, atol=1.0e-7):
        raise ValueError("cell centers differ inside one diagnostic run")
    result = {"height_km": h_left}
    for metric, _ in METRICS:
        result[metric] = (1.0 - weight) * left[metric] + weight * right[metric]
    return result


def ns_from_path(path: Path, fallback: int) -> int:
    match = re.search(r"ns(\d+)", path.stem)
    return int(match.group(1)) if match else fallback


def load_run(path: Path, target_time: float) -> dict[str, object]:
    low, high = time_bracket(sample_times(path), target_time)
    samples = read_profiles(path, {low[1], high[1]})
    weight = 0.0 if low[0] == high[0] else (target_time - low[0]) / (high[0] - low[0])
    profile = interpolate_profile(samples[low[1]], samples[high[1]], weight)
    return {
        "path": path,
        "ns": ns_from_path(path, len(profile["height_km"])),
        "ds_km": float(np.median(np.diff(profile["height_km"]))),
        "t_low": low[0],
        "t_high": high[0],
        "stage_low": low[1],
        "stage_high": high[1],
        "profile": profile,
    }


def linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    slope, intercept = np.polyfit(x, y, 1)
    prediction = intercept + slope * x
    total = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1.0 if total == 0.0 else 1.0 - np.sum((y - prediction) ** 2) / total
    return float(intercept), float(slope), float(r_squared)


def window_metrics(h: np.ndarray, values: np.ndarray, h0: float, h1: float) -> tuple[float, float, float, float]:
    mask = (h >= h0) & (h <= h1)
    h_window, y_window = h[mask], values[mask]
    index = int(np.argmax(np.abs(y_window)))
    return (
        float(np.max(np.abs(y_window))),
        float(np.trapezoid(np.abs(y_window), h_window)),
        float(h_window[index]),
        float(y_window[index]),
    )


def plot(runs: list[dict[str, object]], target_time: float, h0: float, h1: float, out_path: Path) -> None:
    runs = sorted(runs, key=lambda item: float(item["ds_km"]), reverse=True)
    colors = ["C0", "C1", "C2", "C3", "C4"]
    fig, axes = plt.subplots(len(METRICS), 2, figsize=(13.0, 13.6), constrained_layout=True)
    fig.suptitle(
        f"Total-mixture dynamic-HSE diagnostic at t={target_time:g} s (linear time interpolation)",
        fontsize=14,
    )

    for row, (metric, ylabel) in enumerate(METRICS):
        profile_ax, fit_ax = axes[row]
        ds_values: list[float] = []
        values_250: list[float] = []
        print(f"{metric} window [{h0:g}, {h1:g}] km:")
        for color, run in zip(colors, runs):
            profile = run["profile"]
            assert isinstance(profile, dict)
            h = profile["height_km"]
            values = profile[metric]
            ns = int(run["ns"])
            ds_km = float(run["ds_km"])
            mask = (h >= h0) & (h <= h1)
            profile_ax.plot(h[mask], values[mask], color=color, linewidth=1.5,
                            label=fr"ns={ns}, $\Delta s$={ds_km:.3f} km")
            value_250 = float(np.interp(250.0, h, values))
            maxabs, l1, hpeak, peak = window_metrics(h, values, h0, h1)
            ds_values.append(ds_km)
            values_250.append(value_250)
            print(
                f"  ns={ns}: ds={ds_km:.9f} km, at250={value_250:.9e}, "
                f"maxabs={maxabs:.9e}, L1={l1:.9e}, peak_h={hpeak:.6f}, peak={peak:.9e}"
            )

        x, y = np.asarray(ds_values), np.asarray(values_250)
        intercept, slope, r_squared = linear_fit(x, y)
        xfit = np.linspace(0.0, 1.05 * np.max(x), 100)
        fit_ax.plot(x, y, "o", color="C0", markersize=6, label="matched-time value")
        fit_ax.plot(xfit, intercept + slope * xfit, "--", color="C0", linewidth=1.3,
                    label=fr"fit: {intercept:.3e} + {slope:.3e}$\Delta s$")
        fit_ax.scatter([0.0], [intercept], marker="*", s=110, color="crimson", zorder=4,
                       label=fr"$\Delta s\to0$: {intercept:.3e}")
        fit_ax.axhline(0.0, color="0.35", linewidth=0.7)
        fit_ax.set_xlim(left=0.0)
        fit_ax.set_xlabel(r"grid spacing $\Delta s$ [km]")
        fit_ax.set_ylabel(f"value at 250 km\n{ylabel}")
        fit_ax.set_title(fr"250 km linear fit ($R^2$={r_squared:.4f})", fontsize=10)
        fit_ax.legend(fontsize=8, loc="best")
        fit_ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

        profile_ax.axvline(250.0, color="0.35", linewidth=0.8, linestyle=":")
        profile_ax.axhline(0.0, color="0.35", linewidth=0.7)
        profile_ax.set_xlim(h0, h1)
        profile_ax.set_xlabel("height [km]")
        profile_ax.set_ylabel(ylabel)
        profile_ax.set_title("matched-time lower-chromosphere profiles", fontsize=10)
        profile_ax.legend(fontsize=8, loc="best")
        profile_ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        print(f"  fit: intercept={intercept:.9e}, slope={slope:.9e} per km, R2={r_squared:.8f}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    print(f"wrote {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="*", type=Path, default=[Path(path) for path in DEFAULT_FILES])
    parser.add_argument("--time", type=float, default=100.0)
    parser.add_argument("--h0", type=float, default=0.0, help="lower-window start [km]")
    parser.add_argument("--h1", type=float, default=700.0, help="lower-window end [km]")
    parser.add_argument("--out", type=Path, default=Path("visualization/iso_dynamic_mix_hse_t100.png"))
    args = parser.parse_args()

    if len(args.files) < 2:
        parser.error("provide at least two diagnostic CSV files")
    missing = [str(path) for path in args.files if not path.is_file() or path.stat().st_size == 0]
    if missing:
        parser.error("missing or empty diagnostic CSV: " + ", ".join(missing))

    runs = [load_run(path, args.time) for path in args.files]
    for run in sorted(runs, key=lambda item: int(item["ns"])):
        print(
            f"ns={int(run['ns'])}: ds={float(run['ds_km']):.9f} km; "
            f"bracket t=[{float(run['t_low']):.6f}, {float(run['t_high']):.6f}] s "
            f"at stages [{int(run['stage_low'])}, {int(run['stage_high'])}]"
        )
    plot(runs, args.time, args.h0, args.h1, args.out)


if __name__ == "__main__":
    main()
