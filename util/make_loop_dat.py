#!/usr/bin/env python3
"""Generate a full chromosphere→corona loop .dat for the pfss_field_line scenario.

This is Phase-1 "Option A": instead of the truncated ~1 Mm chromospheric stub that
extract_field_line.py currently emits, build the WHOLE field line — Model C7
chromosphere + a transition region + a hot, hydrostatic, fully-ionized corona out
to the loop apex (closed) or coronal top (open) — so flare-evaporated plasma has a
real corona to fill and rest against (fixes the truncated-domain artifact).

Geometry:
  * closed: semicircular loop, model one leg footpoint→apex (HALF loop). Field-
    aligned gravity g_∥(s) = g·cos(s/R) with R = 2L/π (full at footpoint, zero at
    apex). The apex is a symmetry point → the C++ scenario applies a reflecting BC
    there. EXACT only for a symmetric flare (both footpoints heated identically);
    for one-sided heating it bounces the jet at the apex instead of letting it
    cross to the far leg — use `full` instead.
  * full:   the WHOLE loop footpoint A→apex→footpoint B (arc 0…2L), built by
    mirroring one leg about the apex. φ_g is a symmetric "hill" (g_∥ = −∇φ_g flips
    sign on the descending leg automatically), and BOTH ends are chromospheric
    footpoints → both get the reflecting wall, the apex is an ordinary interior
    cell. This removes the instant-reflection timing artifact: an asymmetric
    evaporation jet streams THROUGH the apex and down the far leg with the correct
    transit time. The beam (keyed on arc length from cell 0) lands on footpoint A
    only ⇒ asymmetric heating, the realistic flare-ribbon case.
  * open:   vertical/radial rise, g_∥(s) = g. The C++ scenario applies a coronal
    outflow BC at the top.

Corona is built fully ionized (n_n = 1e-3 n_e ⇒ ionization fraction ≈ 0.999 >
F_SLAVE=0.99), so apply_flare_floor_stage collapses it to single fluid (U=V,
T_n=T_e) — the single-fluid corona limit. The chromosphere keeps C7's partially
ionized two-fluid state.

B(s) is held uniform at the footpoint value for now (no flux-tube expansion);
PFSS-derived B(s) is a later refinement. Writes the [META]/[CELLS]/[GHOSTS]
format scenarios/data_file_parser.cpp reads, plus a topology= key in [META].

Usage:
  python util/make_loop_dat.py --output scenarios/data/loop_closed_test.dat \
      --topology closed --half-length-mm 25 --apex-T-MK 2.0 --fp-B-G 100
"""
import argparse
import os

import numpy as np

from extract_field_line import MODEL_C7, interp_c7, G_SI, H_C7_BASE_KM

M_I = 1.6726219e-27
M_N = M_I
K_B = 1.380649e-23

CHROMO_TOP_KM = float(MODEL_C7[-1, 0])      # 1989 km (top of the Python C7 table)
CHROMO_TOP_T = float(MODEL_C7[-1, 3])       # 6674 K
CHROMO_TOP_RISE_M = (CHROMO_TOP_KM - H_C7_BASE_KM) * 1.0e3   # rise above footpoint [m]
F_CORONA = 1.0e-3                            # n_n/n_e in the corona ⇒ f_ion ≈ 0.999


def stretched_faces(L, ds_fine, s_fine_top, growth, ds_max):
    """Face positions [m] from 0 to L: uniform ds_fine up to s_fine_top, then
    geometric growth to ds_max."""
    faces = [0.0]
    d = ds_fine
    while faces[-1] < L - 1e-6:
        s = faces[-1]
        if s >= s_fine_top:
            d = min(d * growth, ds_max)
        faces.append(min(s + d, L))
    return np.array(faces)


def rise_of_s(s, L, topology):
    """Height gained above the footpoint [m] at arc length s [m]."""
    if topology == "closed":
        R = 2.0 * L / np.pi
        return R * np.sin(s / R)
    return s.copy()  # open: vertical


def g_parallel(s, L, topology):
    """Field-aligned gravity [m/s^2] at arc length s."""
    if topology == "closed":
        R = 2.0 * L / np.pi
        return G_SI * np.cos(s / R)
    return np.full_like(s, G_SI)


def tr_corona_temperature(rise, apex_T, tr_thickness_m):
    """T [K] above the chromosphere top: cosine ramp C7-top→apex_T over a TR of
    thickness tr_thickness_m (in rise), then flat at apex_T."""
    x = rise - CHROMO_TOP_RISE_M
    frac = np.clip(x / tr_thickness_m, 0.0, 1.0)
    ramp = 0.5 * (1.0 - np.cos(np.pi * frac))
    return CHROMO_TOP_T + (apex_T - CHROMO_TOP_T) * ramp


def build_profile(faces, L, topology, apex_T, tr_thickness_m):
    """Return per-cell (ds, phi_g_imh, phi_g_iph, ne, nn, T) and cell-center rise."""
    ds = np.diff(faces)
    rise_f = rise_of_s(faces, L, topology)
    phi_g_f = G_SI * rise_f                         # potential, gauge-zeroed at footpoint
    phi_g_imh = phi_g_f[:-1]
    phi_g_iph = phi_g_f[1:]
    s_ctr = 0.5 * (faces[:-1] + faces[1:])
    rise_c = rise_of_s(s_ctr, L, topology)
    g_par = g_parallel(s_ctr, L, topology)

    ns = ds.size
    ne = np.zeros(ns); nn = np.zeros(ns); T = np.zeros(ns)

    chromo = rise_c <= CHROMO_TOP_RISE_M
    # --- chromosphere: Model C7 (two-fluid, partially ionized) ----------------
    h_cell_m = H_C7_BASE_KM * 1.0e3 + rise_c[chromo]
    ne_c, nn_c, T_c = interp_c7(h_cell_m)
    ne[chromo] = ne_c; nn[chromo] = nn_c; T[chromo] = T_c

    # --- TR + corona: prescribed T ramp, hydrostatic p, fully ionized ---------
    T[~chromo] = tr_corona_temperature(rise_c[~chromo], apex_T, tr_thickness_m)

    # Pressure continuity at the chromosphere top: total p (ion 2x + neutral).
    i_top = np.where(chromo)[0][-1]
    p_junction = (2.0 * ne[i_top] + nn[i_top]) * K_B * T[i_top]

    # March hydrostatic pressure upward through the TR+corona cells.
    p = p_junction
    s_prev = s_ctr[i_top]
    for i in range(i_top + 1, ns):
        dln = -(M_I * g_par[i] / (2.0 * K_B * T[i])) * (s_ctr[i] - s_prev)
        p = p * np.exp(dln)
        ne[i] = p / (2.0 * K_B * T[i])
        nn[i] = F_CORONA * ne[i]
        s_prev = s_ctr[i]

    return ds, phi_g_imh, phi_g_iph, ne, nn, T, rise_c


def build_full_loop(half_faces, L, apex_T, tr_thickness_m):
    """Whole loop footpoint A→apex→footpoint B by mirroring one leg about the apex.

    Build the ascending leg [0, L] with build_profile (C7 chromosphere + hydrostatic
    TR/corona), then reflect every cell onto the descending leg [L, 2L]. The result
    is symmetric in (ds, ne, nn, T) and rise(s); φ_g = g·rise is a symmetric hill
    peaking at the apex, so the code's gravity (gb = −∇φ_g) decelerates upflow on
    the ascending leg and pulls toward footpoint B on the descending leg with no
    special-casing. Both end faces are footpoints → reflecting walls; the apex is an
    interior cell."""
    ds_h, _, _, ne_h, nn_h, T_h, _ = build_profile(half_faces, L, "closed",
                                                    apex_T, tr_thickness_m)
    # Mirror cells (the two coarse apex cells, one per leg, straddle the apex).
    ds = np.concatenate([ds_h, ds_h[::-1]])
    ne = np.concatenate([ne_h, ne_h[::-1]])
    nn = np.concatenate([nn_h, nn_h[::-1]])
    T  = np.concatenate([T_h,  T_h[::-1]])

    # Rebuild faces & the symmetric potential hill over the full arc [0, 2L].
    faces = np.concatenate([[0.0], np.cumsum(ds)])     # ends at 2L
    R = 2.0 * L / np.pi
    rise_f = R * np.sin(np.clip(faces / R, 0.0, np.pi))
    phi_g_f = G_SI * rise_f
    phi_g_imh = phi_g_f[:-1]
    phi_g_iph = phi_g_f[1:]
    s_ctr = 0.5 * (faces[:-1] + faces[1:])
    rise_c = R * np.sin(np.clip(s_ctr / R, 0.0, np.pi))
    return ds, phi_g_imh, phi_g_iph, ne, nn, T, rise_c


def write_dat(path, topology, L, fp_B_T, apex_T,
              ds, phi_g_imh, phi_g_iph, ne, nn, T):
    ns = ds.size
    B = float(fp_B_T)
    with open(path, "w") as f:
        f.write("# full-loop .dat — generated by util/make_loop_dat.py\n")
        f.write(f"# topology={topology}, half-length={L/1e6:.1f} Mm, "
                f"apex_T={apex_T/1e6:.2f} MK, footpoint |B|={B*1e4:.1f} G (uniform)\n")
        f.write("# chromosphere: Model C7 (two-fluid); TR+corona: hydrostatic, "
                "fully ionized (single-fluid via F_SLAVE)\n")
        f.write("# Units: SI (m, T, K, m^-3, J/kg)\n")
        f.write("\n[META]\n")
        f.write(f"ns={ns}\n")
        f.write(f"g_si={G_SI}\n")
        f.write(f"B_outer_T={B:.6e}\n")
        f.write("phi_g_offset_Jpkg=0.0\n")
        f.write(f"topology={topology}\n")
        f.write(f"loop_half_length_m={L:.6e}\n")
        f.write("\n[CELLS]\n")
        f.write("# i  ds_m         B_imh_T      B_iph_T      "
                "phi_g_imh    phi_g_iph    ne_im3       nn_im3       T_K\n")
        for i in range(ns):
            f.write(f"{i}  {ds[i]:.6e}  {B:.6e}  {B:.6e}  "
                    f"{phi_g_imh[i]:.6e}  {phi_g_iph[i]:.6e}  "
                    f"{ne[i]:.6e}  {nn[i]:.6e}  {T[i]:.6e}\n")
        # Ghosts. Outer = top cell (apex for closed / coronal top for open);
        # inner = base cell. The C++ BC overrides these each step.
        f.write("\n[GHOSTS]\n")
        f.write("# tag      B_T          phi_g_Jpkg   ne_im3       nn_im3       T_K          T_e_factor\n")
        outer = dict(B=B, phi_g=float(phi_g_iph[-1]), ne=float(ne[-1]), nn=float(nn[-1]),
                     T=float(T[-1]), tef=1.0)
        inner = dict(B=B, phi_g=float(phi_g_imh[0]), ne=float(ne[0]), nn=float(nn[0]),
                     T=float(T[0]), tef=1.0)
        for tag, g in [("outer_0", outer), ("outer_1", outer),
                       ("inner_0", inner), ("inner_1", inner)]:
            f.write(f"{tag}  {g['B']:.6e}  {g['phi_g']:.6e}  "
                    f"{g['ne']:.6e}  {g['nn']:.6e}  {g['T']:.6e}  {g['tef']:.3f}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", required=True)
    ap.add_argument("--topology", choices=["closed", "full", "open"], default="closed",
                    help="closed=half loop (reflecting apex, symmetric flares only); "
                         "full=whole loop A→apex→B (asymmetry-capable); open=coronal outflow")
    ap.add_argument("--half-length-mm", type=float, default=25.0,
                    help="footpoint→apex arc length [Mm] (closed/full) or total length (open)")
    ap.add_argument("--apex-T-MK", type=float, default=2.0, help="coronal temperature [MK]")
    ap.add_argument("--fp-B-G", type=float, default=100.0, help="footpoint |B| [G], uniform")
    ap.add_argument("--tr-thickness-km", type=float, default=400.0,
                    help="transition-region thickness in rise [km]")
    ap.add_argument("--ds-fine-km", type=float, default=11.0)
    ap.add_argument("--s-fine-top-mm", type=float, default=1.5,
                    help="arc length below which the grid stays fine [Mm]")
    ap.add_argument("--ds-max-km", type=float, default=400.0)
    ap.add_argument("--growth", type=float, default=1.06)
    args = ap.parse_args()

    L = args.half_length_mm * 1.0e6
    apex_T = args.apex_T_MK * 1.0e6
    fp_B_T = args.fp_B_G * 1.0e-4

    faces = stretched_faces(L, args.ds_fine_km * 1e3, args.s_fine_top_mm * 1e6,
                            args.growth, args.ds_max_km * 1e3)
    if args.topology == "full":
        ds, phi_g_imh, phi_g_iph, ne, nn, T, rise_c = build_full_loop(
            faces, L, apex_T, args.tr_thickness_km * 1e3)
    else:
        ds, phi_g_imh, phi_g_iph, ne, nn, T, rise_c = build_profile(
            faces, L, args.topology, apex_T, args.tr_thickness_km * 1e3)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    write_dat(args.output, args.topology, L, fp_B_T, apex_T,
              ds, phi_g_imh, phi_g_iph, ne, nn, T)

    f_ion = ne / (ne + nn)
    print(f"[make_loop] {args.output}")
    print(f"  topology={args.topology}  ns={ds.size}  L={L/1e6:.1f} Mm  "
          f"total arc={ds.sum()/1e6:.1f} Mm  apex height={rise_c.max()/1e6:.2f} Mm")
    print(f"  grid: ds {ds.min()/1e3:.1f}–{ds.max()/1e3:.1f} km, "
          f"fine cells (ds<20km)={int((ds<20e3).sum())}, coronal cells={int((ds>=20e3).sum())}")
    ia = int(np.argmax(rise_c))                 # true apex (max rise), not last cell
    print(f"  chromosphere base: ne={ne[0]:.2e} nn={nn[0]:.2e} T={T[0]:.0f} K  f={f_ion[0]:.3f}")
    print(f"  corona apex:       ne={ne[ia]:.2e} nn={nn[ia]:.2e} T={T[ia]:.3e} K  f={f_ion[ia]:.3f}")


if __name__ == "__main__":
    main()
