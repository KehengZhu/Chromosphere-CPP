#!/usr/bin/env python3
"""Preliminary two-fluid decoupling diagnostic — measures how far the neutrals
slip away from the ions in an EXISTING chromo_main run (no recompile needed).

The code is already a two-fluid model: it advances ions (rho_i, V, e_i) and
neutrals (rho_n, U, e_n) separately, coupled by ion-neutral drag (Stage B) and
T_i<->T_n equilibration (Stage C), both scaling with nu_in ~ n_n*sqrt(T).

The single-fluid limit is the strong-coupling limit nu_in -> inf, where V->U and
T_i->T_n at every step. So the size of the *velocity drift* w = V - U and the
*temperature split* dT = T_i - T_n in a baseline run directly measures how far
the real (finite-coupling) two-fluid solution departs from single-fluid.

This reads the 6-component conserved output [rho_i, rho_n, rho_i V, rho_n U,
e_i, e_n] and the .dat (for phi_g), reconstructs primitives EXACTLY as
animate_loop_flare.primitives does, and reports where/when/how-much the fluids
decouple, plus an order-of-magnitude frictional-heating budget.
"""
import os
import sys
import numpy as np

m_i = 1.6726219e-27
m_n = m_i
k_b = 1.380649e-23
bohr_r = 53e-12
GAMMA = 5.0 / 3.0


def parse_dat(path):
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
    s_arc = np.cumsum(ds) - 0.5 * ds
    return ds, s_arc, phi_c


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


def main():
    dat = sys.argv[1]
    out = sys.argv[2]
    tag = sys.argv[3] if len(sys.argv) > 3 else "run"

    ds, s_arc, phi_g = parse_dat(dat)
    frames = read_frames(out)
    ns = len(ds)
    nt = len(frames)
    print(f"[{tag}] {nt} frames x {ns} cells; loop arc = {s_arc[-1]/1e6:.1f} Mm")

    t = np.array([fr[0] for fr in frames])
    # per (time, cell) fields
    W = np.zeros((nt, ns))      # drift V - U  [m/s]
    V = np.zeros((nt, ns))      # ion velocity [m/s]
    CS = np.zeros((nt, ns))     # ion sound speed [m/s]
    DT = np.zeros((nt, ns))     # T_i - T_n [K]
    Ti = np.zeros((nt, ns))
    F = np.zeros((nt, ns))      # ionization fraction
    QF = np.zeros((nt, ns))     # frictional heating rate alpha*w^2 [W/m^3]
    NN = np.zeros((nt, ns))     # neutral number density

    for k, (_, xn) in enumerate(frames):
        rho_i = xn[:, 0]; rho_n = xn[:, 1]
        v = xn[:, 2] / rho_i
        u = xn[:, 3] / rho_n
        n_i = rho_i / m_i; n_n = rho_n / m_n
        p_i = 2/3*xn[:, 4] - 1/3*rho_i*v*v - 2/3*rho_i*phi_g
        p_n = 2/3*xn[:, 5] - 1/3*rho_n*u*u - 2/3*rho_n*phi_g
        T_i = p_i / (2.0 * n_i * k_b)
        T_n = p_n / (n_n * k_b)
        nu = (2*bohr_r)**2 * n_n * np.sqrt(np.abs(8*np.pi*k_b*(T_i+T_n)/m_i))
        alpha = rho_i * nu
        W[k] = v - u
        V[k] = v
        CS[k] = np.sqrt(np.abs(GAMMA * p_i / rho_i))
        DT[k] = T_i - T_n
        Ti[k] = T_i
        F[k] = rho_i / (rho_i + rho_n)
        QF[k] = alpha * (v - u)**2
        NN[k] = n_n

    cs_floor = np.maximum(CS, 1.0e3)              # avoid div by ~0 in cold dense gas
    slip = np.abs(W) / cs_floor                   # slip Mach number |w|/c_s
    dT_rel = np.abs(DT) / np.maximum(Ti, 1.0e3)   # |dT|/T_i

    def loc(arr, label, unit, scale=1.0):
        kmax, imax = np.unravel_index(np.argmax(arr), arr.shape)
        print(f"  max {label:18s} = {arr[kmax,imax]*scale:10.3g} {unit:7s} "
              f"at s={s_arc[imax]/1e6:6.2f} Mm, t={t[kmax]:5.1f}s "
              f"(f={F[kmax,imax]:.2e}, T_i={Ti[kmax,imax]:.2e}K, n_n={NN[kmax,imax]:.2e}m^-3)")

    print("  --- velocity decoupling -------------------------------------------")
    loc(np.abs(W), "|V-U|", "km/s", 1e-3)
    loc(slip, "|V-U|/c_s,ion", "(Mach)")
    print("  --- thermal decoupling --------------------------------------------")
    loc(np.abs(DT), "|T_i-T_n|", "K")
    loc(dT_rel, "|T_i-T_n|/T_i", "")
    print("  --- frictional (drift) heating ------------------------------------")
    # crude energy budget: integrate Q_fric over volume(~area*ds) and time.
    # Use per-area column (drop the cross-section A, same for both channels).
    dt_frame = np.gradient(t)
    E_fric = np.sum(QF * ds[None, :] * dt_frame[:, None])   # J/m^2 (column)
    print(f"  column-integrated frictional heat   = {E_fric:.3e} J/m^2")
    # beam energy per area over the run (if beam flux known via env-default 3e7):
    beam_flux = 3.0e7   # W/m^2 (M8.2 default); over the heating window
    print(f"  (for scale: M8.2 beam 3e7 W/m^2 x 10 s = {beam_flux*10:.2e} J/m^2)")
    print(f"  frictional / beam-input ratio       ~ {E_fric/(beam_flux*10):.2e}")

    # fraction of (cell,time) samples that are 'decoupled' by a 10% threshold
    frac_slip = np.mean(slip > 0.1)
    frac_dT = np.mean(dT_rel > 0.1)
    print("  --- prevalence ----------------------------------------------------")
    print(f"  cells*frames with |w|/c_s > 0.1 : {100*frac_slip:5.2f} %")
    print(f"  cells*frames with |dT|/T  > 0.1 : {100*frac_dT:5.2f} %")

    # optional heatmap
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
        xMm = s_arc / 1e6
        im0 = ax[0].pcolormesh(xMm, t, slip, cmap="inferno", vmin=0, vmax=0.5, shading="auto")
        ax[0].set_ylabel("t [s]"); ax[0].set_title(f"{tag}: slip Mach |V-U|/c_s,ion")
        fig.colorbar(im0, ax=ax[0])
        im1 = ax[1].pcolormesh(xMm, t, dT_rel, cmap="viridis", vmin=0, vmax=0.5, shading="auto")
        ax[1].set_ylabel("t [s]"); ax[1].set_xlabel("arc length s [Mm]")
        ax[1].set_title("temperature split |T_i-T_n|/T_i")
        fig.colorbar(im1, ax=ax[1])
        fig.tight_layout()
        png = f"visualization/two_fluid/two_fluid_diag_{tag}.png"
        os.makedirs(os.path.dirname(png), exist_ok=True)
        fig.savefig(png, dpi=120)
        print(f"  wrote {png}")
    except Exception as e:
        print(f"  (plot skipped: {e})")


if __name__ == "__main__":
    main()
