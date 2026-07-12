#!/usr/bin/env python3
"""Compare a NEUTRALS-ON (two-fluid) run against a NEUTRALS-OFF (single-fluid
limit, SINGLE_FLUID=1) run of the SAME field line. Same IC / beam / BCs, so the
difference is purely the ion-neutral decoupling.

Reconstructs ion temperature T_i and ion velocity V exactly as
animate_loop_flare.primitives, then reports:
  * peak upflow velocity (evaporation) on vs off
  * peak apex temperature on vs off
  * max |T_on - T_off| and |V_on - V_off| with location/time
  * relative RMS difference of the T and V fields over (s,t)
and writes a difference figure: dT(s,t), dV(s,t) heatmaps + overlay profiles at
the epoch of maximum departure.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

m_i = 1.6726219e-27
m_n = m_i
k_b = 1.380649e-23


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


def fields(frames, phi_g):
    nt, ns = len(frames), len(phi_g)
    T = np.zeros((nt, ns)); Vv = np.zeros((nt, ns))
    NE = np.zeros((nt, ns)); t = np.zeros(nt)
    for k, (tt, xn) in enumerate(frames):
        rho_i = xn[:, 0]
        v = xn[:, 2] / rho_i
        n_i = rho_i / m_i
        p_i = 2/3*xn[:, 4] - 1/3*rho_i*v*v - 2/3*rho_i*phi_g
        T[k] = p_i / (2.0 * n_i * k_b)
        Vv[k] = v
        NE[k] = n_i
        t[k] = tt
    return t, T, Vv, NE


def main():
    dat = sys.argv[1]
    on_txt = sys.argv[2]
    off_txt = sys.argv[3]
    tag = sys.argv[4] if len(sys.argv) > 4 else "line"

    ds, s_arc, phi_g = parse_dat(dat)
    fon = read_frames(on_txt)
    foff = read_frames(off_txt)
    nt = min(len(fon), len(foff))
    fon, foff = fon[:nt], foff[:nt]
    t, T_on, V_on, NE_on = fields(fon, phi_g)
    _, T_off, V_off, NE_off = fields(foff, phi_g)
    xMm = s_arc / 1e6
    print(f"\n=== {tag}: neutrals ON (2-fluid) vs OFF (single-fluid) ===")
    print(f"  {nt} frames x {len(ds)} cells; arc = {xMm[-1]:.1f} Mm; t_end = {t[-1]:.1f} s")

    dT = T_on - T_off
    dV = V_on - V_off
    # upflow = positive V into corona; report peak |V| each (evaporation jet)
    print("  --- bulk evaporation observables ----------------------------------")
    print(f"  peak |V| (upflow/jet)   ON = {np.abs(V_on).max()/1e3:8.2f} km/s   "
          f"OFF = {np.abs(V_off).max()/1e3:8.2f} km/s   "
          f"Δ = {100*(np.abs(V_on).max()-np.abs(V_off).max())/np.abs(V_off).max():+5.1f} %")
    print(f"  peak T (apex/loop)      ON = {T_on.max()/1e6:8.3f} MK     "
          f"OFF = {T_off.max()/1e6:8.3f} MK     "
          f"Δ = {100*(T_on.max()-T_off.max())/T_off.max():+5.1f} %")
    # evaporated mass proxy: max coronal column electron density change
    print(f"  peak n_e (any cell)     ON = {NE_on.max():.3e}  OFF = {NE_off.max():.3e}  "
          f"Δ = {100*(NE_on.max()-NE_off.max())/NE_off.max():+5.1f} %")

    print("  --- field-wide departure (single-fluid error) ---------------------")
    kT, iT = np.unravel_index(np.argmax(np.abs(dT)), dT.shape)
    kV, iV = np.unravel_index(np.argmax(np.abs(dV)), dV.shape)
    print(f"  max |T_on - T_off| = {np.abs(dT).max():.3e} K  "
          f"({100*np.abs(dT)[kT,iT]/T_on[kT,iT]:.1f}% of T) at s={xMm[iT]:.2f}Mm t={t[kT]:.1f}s")
    print(f"  max |V_on - V_off| = {np.abs(dV).max()/1e3:.3f} km/s "
          f"at s={xMm[iV]:.2f}Mm t={t[kV]:.1f}s")
    rmsT = np.sqrt(np.mean(dT**2)) / np.sqrt(np.mean(T_on**2))
    rmsV = np.sqrt(np.mean(dV**2)) / (np.sqrt(np.mean(V_on**2)) + 1.0)
    print(f"  relative RMS(T) departure over (s,t) = {100*rmsT:.3f} %")
    print(f"  relative RMS(V) departure over (s,t) = {100*rmsV:.3f} %")

    # ---- figure ----
    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0])
    axTm = fig.add_subplot(gs[0, 0]); axVm = fig.add_subplot(gs[0, 1])
    axTp = fig.add_subplot(gs[1, 0]); axVp = fig.add_subplot(gs[1, 1])

    Tmax = max(1.0, np.abs(dT).max())
    im0 = axTm.pcolormesh(xMm, t, dT, cmap="RdBu_r", vmin=-Tmax, vmax=Tmax, shading="auto")
    axTm.set_title("T_on − T_off  [K]"); axTm.set_ylabel("t [s]"); axTm.set_xlabel("s [Mm]")
    fig.colorbar(im0, ax=axTm)
    Vmax = max(1.0, np.abs(dV).max()) / 1e3
    im1 = axVm.pcolormesh(xMm, t, dV/1e3, cmap="RdBu_r", vmin=-Vmax, vmax=Vmax, shading="auto")
    axVm.set_title("V_on − V_off  [km/s]"); axVm.set_xlabel("s [Mm]")
    fig.colorbar(im1, ax=axVm)

    # overlay profiles at the epoch of max T departure, zoomed to that region
    def zoom(ax, i0):
        lo = max(0, xMm[i0] - 4.0); hi = min(xMm[-1], xMm[i0] + 4.0)
        ax.set_xlim(lo, hi)
    axTp.plot(xMm, T_on[kT]/1e6, "-", color="C0", lw=1.8, label="neutrals ON (2-fluid)")
    axTp.plot(xMm, T_off[kT]/1e6, "--", color="C3", lw=1.5, label="neutrals OFF (1-fluid)")
    axTp.axvline(xMm[iT], color="gray", ls=":", lw=0.8)
    zoom(axTp, iT)
    axTp.set_title(f"T profile @ t={t[kT]:.1f}s (max ΔT epoch)")
    axTp.set_xlabel("s [Mm]"); axTp.set_ylabel("T [MK]"); axTp.legend(fontsize=8)
    axVp.plot(xMm, V_on[kV]/1e3, "-", color="C0", lw=1.8, label="neutrals ON")
    axVp.plot(xMm, V_off[kV]/1e3, "--", color="C3", lw=1.5, label="neutrals OFF")
    axVp.axvline(xMm[iV], color="gray", ls=":", lw=0.8)
    zoom(axVp, iV)
    axVp.set_title(f"V profile @ t={t[kV]:.1f}s (max ΔV epoch)")
    axVp.set_xlabel("s [Mm]"); axVp.set_ylabel("V [km/s]"); axVp.legend(fontsize=8)

    fig.suptitle(f"{tag}: two-fluid (ON) vs single-fluid (OFF) — same beam/IC/BCs", fontsize=12)
    fig.tight_layout()
    png = f"visualization/two_fluid/compare_{tag}.png"
    fig.savefig(png, dpi=120)
    print(f"  wrote {png}")


if __name__ == "__main__":
    main()
