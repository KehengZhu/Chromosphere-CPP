#!/usr/bin/env python3
"""Plot a single frame (default: final) of a model_isentropic run as a static PNG.

Same 6-panel layout as util/animate_isentropic.py, but for a single time. The mass
flux panel ρV is autoscaled to THIS frame — the movie's global ρV range is dominated
by the early conduction-switch-on transient, so on the shared scale the quasi-steady
final frame looks flat near zero; here it is zoomed to its own min/max.

Usage:
    python util/plot_iso_frame.py [input.txt] [frame_index] [output.png]
frame_index: integer, negative allowed (default -1 = final frame).
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from animate_isentropic import load, primitives, heat_flux, G


def main():
    in_path = sys.argv[1] if len(sys.argv) > 1 else "outputs/iso_t22k.txt"
    fidx = int(sys.argv[2]) if len(sys.argv) > 2 else -1
    if len(sys.argv) > 3:
        out_path = sys.argv[3]
    else:
        stem = os.path.splitext(os.path.basename(in_path))[0]
        out_path = os.path.join("visualization", f"{stem}_final_frame.png")
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    h, frames = load(in_path)
    phi = G * (h - h[0]) * 1000.0
    t, arr = frames[fidx]
    V, T, rho, p, n_e, n_n = primitives(arr, phi)
    q = heat_flux(h, T, n_e, n_n)
    mflux = rho * V

    stem = os.path.splitext(os.path.basename(in_path))[0]
    fig, ax = plt.subplots(2, 3, figsize=(17, 8))
    fig.suptitle(f"model_isentropic — {stem}  |  final frame  t = {t:.1f} s", fontsize=13)

    ax[0, 0].plot(h, T, lw=1.4, color="C3")
    ax[0, 0].set_ylabel("T [K]")
    ax[0, 0].set_title("Temperature")

    ax[0, 1].plot(h, V / 1e3, lw=1.4, color="C0")
    ax[0, 1].axhline(0, c="k", lw=0.5)
    ax[0, 1].set_ylabel("V [km/s]")
    ax[0, 1].set_title("Velocity  (>0 up = evaporation, <0 down = condensation)")

    ax[0, 2].plot(h, q / 1e3, lw=1.4, color="C1")
    ax[0, 2].axhline(0, c="k", lw=0.5)
    ax[0, 2].set_ylabel(r"$q$ [kW/m²]")
    ax[0, 2].set_title("Conductive heat flux  (>0 up; <0 = downward into TR)")

    ax[1, 0].semilogy(h, rho, lw=1.4, color="C2")
    ax[1, 0].set_ylabel(r"$\rho$ [kg/m³]")
    ax[1, 0].set_title("Mass density")

    ax[1, 1].semilogy(h, p, lw=1.4, color="C4")
    ax[1, 1].set_ylabel("p [Pa]")
    ax[1, 1].set_title("Pressure")

    ax[1, 2].plot(h, mflux, lw=1.4, color="C5")
    ax[1, 2].axhline(0, c="k", lw=0.5)
    ax[1, 2].set_ylabel(r"$\rho V$ [kg/m²/s]")
    ax[1, 2].set_title("Mass flux  (>0 up = evaporation, <0 down = condensation)  [zoomed]")
    # Zoom to THIS frame, robustly: use 0.5–99.5 percentiles so the couple of dense
    # base-boundary cells (ρ huge ⇒ a tiny residual V gives a large ρV spike) don't
    # dominate the scale and squash the interior evaporation/drainage structure.
    # (The movie's global ρV scale is instead set by the early conduction transient.)
    finite = mflux[np.isfinite(mflux)]
    lo, hi = (float(x) for x in np.percentile(finite, [0.5, 99.5]))
    d = 0.12 * (hi - lo if hi > lo else max(abs(hi), 1e-30))
    ax[1, 2].set_ylim(lo - d, hi + d)
    ax[1, 2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    for a in ax.ravel():
        a.set_xlim(h.min(), h.max())
        a.set_xlabel("height [km]")
        a.grid(True, alpha=0.3)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print(f"wrote {out_path}  (frame {fidx}, t = {t:.1f} s, "
          f"rhoV in [{lo:.3e}, {hi:.3e}] kg/m^2/s)")


if __name__ == "__main__":
    main()
