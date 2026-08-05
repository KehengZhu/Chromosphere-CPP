#!/usr/bin/env python3
"""Compare fixed-gamma and gamma-table physical conductive heat flux."""

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import animate_isentropic as anim
from conductive_flux import conductivity_components, physical_heat_flux


def frame_at(frames, target):
    """Linearly interpolate stored cell-centred fields to an exact target time."""
    times = np.array([frame[0] for frame in frames])
    j = int(np.searchsorted(times, target))
    if j < len(times) and times[j] == target:
        return frames[j][1]
    if j == 0 or j == len(times):
        raise ValueError(f"target time {target:g} s is outside stored frame range")
    weight = (target-times[j-1])/(times[j]-times[j-1])
    return (1.0-weight)*frames[j-1][1] + weight*frames[j][1]


def metrics(h_km, temperature, n_e, n_hi):
    mask = h_km <= 500.0
    dtds = np.gradient(temperature, h_km*1.0e3)
    kappa_e, kappa_n = conductivity_components(n_e, n_hi, temperature)
    q = physical_heat_flux(h_km, temperature, n_e, n_hi)
    return {
        "max_abs_dTds_K_per_m": float(np.max(np.abs(dtds[mask]))),
        "max_n_e_per_m3": float(np.max(n_e[mask])),
        "max_n_HI_per_m3": float(np.max(n_hi[mask])),
        "max_kappa_e_W_per_m_K": float(np.max(kappa_e[mask])),
        "max_kappa_n_W_per_m_K": float(np.max(kappa_n[mask])),
        "max_abs_q_physical_W_per_m2": float(np.max(np.abs(q[mask]))),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixed", required=True)
    parser.add_argument("--gamma", required=True,
                        help="gamma-table primary output; .gamma_diag is appended")
    parser.add_argument("--time", type=float, default=1000.0)
    parser.add_argument("--fixed-gamma", type=float, default=1.05)
    parser.add_argument("--out", required=True)
    parser.add_argument("--metrics-out")
    args = parser.parse_args()

    h_fixed, fixed_frames = anim.load(args.fixed)
    h_gamma, gamma_frames = anim.load(args.gamma + ".gamma_diag")
    if h_fixed.shape != h_gamma.shape or not np.allclose(h_fixed, h_gamma, rtol=0, atol=1e-6):
        raise ValueError("fixed-gamma and gamma-table height coordinates do not match")

    fixed_state = frame_at(fixed_frames, args.time)
    gamma_diag = frame_at(gamma_frames, args.time)
    anim.GM1 = args.fixed_gamma - 1.0
    phi = anim.G * (h_fixed-h_fixed[0])*1.0e3
    _, tf, _, _, nef, nhif = anim.primitives(fixed_state, phi)
    tg, neg, nhig = gamma_diag[:, 2], gamma_diag[:, 4], gamma_diag[:, 5]

    kef, knf = conductivity_components(nef, nhif, tf)
    keg, kng = conductivity_components(neg, nhig, tg)
    qf = physical_heat_flux(h_fixed, tf, nef, nhif)
    qg = physical_heat_flux(h_gamma, tg, neg, nhig)
    result = {
        "snapshot_time_s": args.time,
        "height_grid_identical": True,
        "below_500_km": {
            "fixed_gamma_1p05": metrics(h_fixed, tf, nef, nhif),
            "gamma_table": metrics(h_gamma, tg, neg, nhig),
        },
    }
    print(json.dumps(result, indent=2))
    if args.metrics_out:
        os.makedirs(os.path.dirname(args.metrics_out) or ".", exist_ok=True)
        with open(args.metrics_out, "w") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")

    blue, orange = "#1f4e79", "#c45118"
    h_mm = h_fixed/1.0e3
    low = h_fixed <= 500.0
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.0), constrained_layout=True)
    fig.suptitle(f"Physical conductive heat flux at t = {args.time:g} s\n"
                 r"$q_\mathrm{phys}=-(\kappa_e+\kappa_n)\,dT/ds$; TRAC, numerical diffusion, and boundary flux excluded")

    axes[0, 0].plot(h_mm, tf, color=blue, label=r"fixed $\gamma=1.05$")
    axes[0, 0].plot(h_mm, tg, color=orange, label="gamma table")
    axes[0, 0].set(xlabel="height [Mm]", ylabel="T [K]", title="Temperature")
    axes[0, 0].legend()

    axes[0, 1].plot(h_mm, qf, color=blue)
    axes[0, 1].plot(h_mm, qg, color=orange)
    axes[0, 1].axhline(0, color="0.25", lw=0.7)
    axes[0, 1].set(xlabel="height [Mm]", ylabel=r"$q_\mathrm{phys}$ [W m$^{-2}$]",
                   title="Physical flux — full column")

    axes[1, 0].plot(h_mm[low], kef[low], color=blue, label=r"$\kappa_e$, fixed")
    axes[1, 0].plot(h_mm[low], knf[low], color=blue, ls="--", label=r"$\kappa_n$, fixed")
    axes[1, 0].plot(h_mm[low], keg[low], color=orange, label=r"$\kappa_e$, table")
    axes[1, 0].plot(h_mm[low], kng[low], color=orange, ls="--", label=r"$\kappa_n$, table")
    axes[1, 0].set_yscale("log")
    axes[1, 0].set(xlabel="height [Mm]", ylabel=r"$\kappa$ [W m$^{-1}$ K$^{-1}$]",
                   title="Conductivity components below 500 km")
    axes[1, 0].legend(ncol=2, fontsize=8)

    axes[1, 1].plot(h_mm[low], qf[low], color=blue)
    axes[1, 1].plot(h_mm[low], qg[low], color=orange)
    axes[1, 1].axhline(0, color="0.25", lw=0.7)
    axes[1, 1].set(xlabel="height [Mm]", ylabel=r"$q_\mathrm{phys}$ [W m$^{-2}$]",
                   title="Physical flux below 500 km")

    for ax in axes.flat:
        ax.grid(alpha=0.25)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
