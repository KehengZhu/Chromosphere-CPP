/*!
 * @file two_fluid/rhs.cpp
 * @brief Explicit (TVD-MUSCL + Rusanov) and implicit (drag, frictional heating,
 *        field-aligned conduction) right-hand sides of the two-fluid solver.
 * @ingroup two_fluid_solver
 */
#include "two_fluid/two_fluid.hpp"
#include "physics.hpp"
#include "profiling.hpp"
#include "parallel.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace chromosphere {

// Well-balanced gravity correction (chromosphere.hpp::well_balanced).
// cons2prim of a SPATIALLY-SHIFTED conserved state recovers the wrong thermal
// pressure: it subtracts the LOCAL-index φ_g, but the shifted cell's total energy
// carries the SHIFTED cell's ρφ_g, so each neighbour pressure is biased by
// ⅔ρ·Δφ_g = ⅔ρg·Δs. That flattens the reconstructed hydrostatic slope dp/ds = −ρg
// to −⅓ρg, giving a hydrostatic atmosphere a spurious resolution-independent
// (γ−1)g downforce. Undo it: subtract ⅔ρ(φ_src − φ_local) from p_i and p_n, where
// φ_src is the cell-centred potential actually baked into the shifted data's
// energy. p_e carries no φ_g (shift-invariant) and is left untouched. Stage-D
// conduction already does the analogous φ_g shift (rhs_implicit_state).
static void wb_correct_shifted(const Grid& grid, Vec& prim,
                               const Vec& phi_src, const Vec& phi_cell) {
    const Vec d     = grid.gm1() * (phi_src - phi_cell);
    const Vec rho_i = get_scalar(grid, prim, prim::RHO_I);
    const Vec rho_n = get_scalar(grid, prim, prim::RHO_N);
    prim -= scalar_to(grid, rho_i % d, prim::P_I);
    prim -= scalar_to(grid, rho_n % d, prim::P_N);
}

// Log-space reconstruction (chromosphere.hpp::log_reconstruct). The exponentially-
// stratified, strictly-positive primitive slots — the densities and pressures — are
// transformed to log before the limiter/extrapolation arithmetic and back (exp) on
// every reconstructed face before prim2cons, so an isothermal-hydrostatic column is
// piecewise-linear and minmod no longer clips its smooth gradient to first order
// (the lower-boundary diffusion source). The SIGNED velocities {V, U} are left linear
// (log undefined for them, and they carry no large dynamic range). A scalar block for
// primitive variable v occupies the contiguous subvec [v·ns, v·ns+ns-1] (cons2prim /
// scalar_to packing). The pre-log clamp guards log(0) at a vacuum cell; exp is safe.
static const arma::uword LOG_PRIM_VARS[5] =
    {prim::RHO_I, prim::RHO_N, prim::P_I, prim::P_N, prim::P_E};

static void prim_to_log(const Grid& grid, Vec& prim) {
    const arma::uword ns = grid.ns;
    for (arma::uword v : LOG_PRIM_VARS) {
        const arma::uword lo = v * ns, hi = lo + ns - 1;
        const Vec blk = prim.subvec(lo, hi);
        prim.subvec(lo, hi) = arma::log(arma::clamp(blk, 1.0e-30f, arma::datum::inf));
    }
}

static void prim_from_log(const Grid& grid, Vec& prim) {
    const arma::uword ns = grid.ns;
    for (arma::uword v : LOG_PRIM_VARS) {
        const arma::uword lo = v * ns, hi = lo + ns - 1;
        const Vec blk = prim.subvec(lo, hi);
        prim.subvec(lo, hi) = arma::exp(blk);
    }
}

// Broadcast a length-ns cell field across all num_of_eq equation slots — the
// same packing scalar_to/get_scalar use — so a per-cell metric vector can
// multiply a packed primitive state elementwise. Used by the non-uniform MUSCL
// reconstruction below.
static Vec bcast_state(const Grid& grid, const Vec& v_i) {
    Vec s(grid.n_state, arma::fill::zeros);
    for (arma::uword k = 0; k < num_of_eq; ++k) s += scalar_to(grid, v_i, k);
    return s;
}

// ============================================================================
// Explicit RHS: TVD-MUSCL reconstruction + predictor-corrector +
// Rusanov flux differencing, minus the explicit source.
// (1.5D field-aligned; z-direction not considered. CoMFi Re_MUSCL analogue.)
// ============================================================================

Vec rhs_explicit_state(const Grid& grid, const Vec& xn_state) {
    if (!grid.eos_gamma_table.empty())
        throw std::logic_error(
            "the legacy two-fluid RHS cannot advance a Gamma1-table Grid; "
            "use the release solver in mixture.hpp");
    // ---- shifted neighbours (conserved) ---------------------------------
    Vec cons_xn_state_ip1 = ip1(grid, xn_state);
    Vec cons_xn_state_im1 = im1(grid, xn_state);
    Vec cons_xn_state_ip2 = ip2(grid, xn_state);
    Vec cons_xn_state_im2 = im2(grid, xn_state);

    // ---- limit on primitive variables -----------------------------------
    Vec prim_xn_state     = cons2prim(grid, xn_state);
    Vec prim_xn_state_ip1 = cons2prim(grid, cons_xn_state_ip1);
    Vec prim_xn_state_im1 = cons2prim(grid, cons_xn_state_im1);
    Vec prim_xn_state_ip2 = cons2prim(grid, cons_xn_state_ip2);
    Vec prim_xn_state_im2 = cons2prim(grid, cons_xn_state_im2);

    // Well-balanced gravity: correct the φ_g inconsistency in every shifted-state
    // pressure (wb_correct_shifted). φ_src is the cell-centred potential baked into
    // the shifted data's energy, matching what the ip*/im* shifts pull in: interior
    // → neighbour φ_cell; the outer/inner ghosts → the face potentials the BCs used
    // when packing the ghost energies (as rhs_implicit_state does for conduction).
    // The same correction is reused for the predictor-shifted states below.
    const Vec phi_cell = 0.5 * (grid.phi_g_imh + grid.phi_g_iph);
    Vec phi_ip1, phi_im1;
    if (grid.well_balanced) {
        const arma::uword ns = grid.ns;
        phi_ip1 = ip1(grid, phi_cell, SLICE); phi_ip1[ns - 1] = grid.phi_g_iph[ns - 1];
        phi_im1 = im1(grid, phi_cell, SLICE); phi_im1[0]      = grid.phi_g_imh[0];
        Vec phi_ip2 = ip1(grid, phi_ip1, SLICE); phi_ip2[ns - 1] = grid.phi_g_iph[ns - 1];
        Vec phi_im2 = im1(grid, phi_im1, SLICE); phi_im2[0]      = grid.phi_g_imh[0];
        wb_correct_shifted(grid, prim_xn_state_ip1, phi_ip1, phi_cell);
        wb_correct_shifted(grid, prim_xn_state_im1, phi_im1, phi_cell);
        wb_correct_shifted(grid, prim_xn_state_ip2, phi_ip2, phi_cell);
        wb_correct_shifted(grid, prim_xn_state_im2, phi_im2, phi_cell);
    }

    // Log-space reconstruction: transform the stratified ρ/p slots to log (AFTER the
    // well-balanced pressure correction above, which needs linear p). All the slope
    // differences, limiter ratios and edge extrapolations below then run in log space
    // for those slots; the faces are exp'd back before prim2cons. V, U stay linear.
    if (grid.log_reconstruct) {
        prim_to_log(grid, prim_xn_state);
        prim_to_log(grid, prim_xn_state_ip1);
        prim_to_log(grid, prim_xn_state_im1);
        prim_to_log(grid, prim_xn_state_ip2);
        prim_to_log(grid, prim_xn_state_im2);
    }

    const Vec dxn_state_iph = prim_xn_state_ip1 - prim_xn_state;
    const Vec dxn_state_imh = prim_xn_state     - prim_xn_state_im1;

    if ((dxn_state_iph.has_nan() + dxn_state_imh.has_nan()) > 0) {
        std::cout << "XN" << std::endl;
        print_xn(grid, xn_state);
        std::cout << "XNPRIM" << std::endl;
        print_xn(grid, prim_xn_state);
        std::cout << "dxn_iph has nan: " << dxn_state_iph.has_nan() << std::endl
                  << "dxn_imh has nan: " << dxn_state_imh.has_nan() << std::endl;
        std::exit(1);
    }

    // limiter inputs
    Vec r_state     = dxn_state_imh / dxn_state_iph;
    Vec r_state_ip1 = dxn_state_iph / (prim_xn_state_ip2 - prim_xn_state_ip1);
    Vec r_state_im1 = (prim_xn_state_im1 - prim_xn_state_im2) / dxn_state_imh;

    // ---- non-uniform-mesh metric corrections ----------------------------
    // The MUSCL slope ratios above compare raw primitive DIFFERENCES (correct
    // only when the cells are equally spaced). On a refined mesh the limiter must
    // instead compare CENTER-TO-CENTER GRADIENTS, and each face is reconstructed
    // with the owning cell's HALF-WIDTH — so the difference terms carry the
    // factor 0.5·ds_i/ds_{i±½}. Both reduce to the legacy form on a uniform mesh
    // (ratio metric = 1, weight = 0.5), so the reconstruction stays byte-for-byte
    // identical there (the uniform branch below runs the original literals).
    //   r_i        gradient ratio ×= ds_{i+½}/ds_{i-½}
    //   r_{i+1}    ×= ds_{i+3/2}/ds_{i+½}   (mirrored-ghost width shift)
    //   r_{i-1}    ×= ds_{i-½}/ds_{i-3/2}
    // Reconstruction weights (state-packed):
    //   W1 = 0.5·ds_i/ds_{i+½}         (owner cell i: Lxn_iph, Rxn_imh)
    //   W3 = ip1(W1) = 0.5·ds_{i+1}/ds_{i+3/2}   (owner i+1: Rxn_iph)
    //   W4 = im1(W1) = 0.5·ds_{i-1}/ds_{i-½}      (owner i-1: Lxn_imh)
    Vec W1, W3, W4;
    if (!grid.uniform_mesh) {
        const Vec w1_i = 0.5f * grid.ds_i / grid.ds_iph_i;
        W1 = bcast_state(grid, w1_i);
        W3 = bcast_state(grid, ip1(grid, w1_i, SLICE));
        W4 = bcast_state(grid, im1(grid, w1_i, SLICE));
        r_state     %= bcast_state(grid, grid.ds_iph_i / grid.ds_imh_i);
        r_state_ip1 %= bcast_state(grid, ip1(grid, grid.ds_iph_i, SLICE) / grid.ds_iph_i);
        r_state_im1 %= bcast_state(grid, grid.ds_imh_i / im1(grid, grid.ds_imh_i, SLICE));
    }

    // Flare scenario: the positivity floors can pin several adjacent cells to the
    // SAME density/pressure value (flat plateaus in the evaporated column), making
    // the slope ratio r = 0/0 = NaN, which propagates through flux_lim into the
    // reconstruction and the Rusanov spectral radius. Replace any non-finite ratio
    // with 0 → flux_lim(0)=0 → first-order (no reconstructed slope) at those flat
    // cells, which is the correct TVD behavior at extrema. Gated on the floors
    // (beam heating or the C7 vacuum floor) so untouched scenarios / the test
    // suite are unchanged.
    if (grid.floors_active()) {
        r_state.elem(arma::find_nonfinite(r_state)).zeros();
        r_state_ip1.elem(arma::find_nonfinite(r_state_ip1)).zeros();
        r_state_im1.elem(arma::find_nonfinite(r_state_im1)).zeros();
    }

    // Slope limiter. Default: symmetric minmod (slimL=slimR ⇒ one φ per cell).
    // grid.mc3_limiter ⇒ the MC3 / Koren limiter, which is ASYMMETRIC: the '+'
    // (right-face, Lxn = u+½φΔ₊) and '−' (left-face, Rxn = u−½φΔ₊) reconstruction
    // terms take different limited slopes, so they get different limiter functions.
    // With mc3 off both lambdas are flux_lim ⇒ byte-identical to the minmod baseline.
    const float beta = grid.limiter_beta;
    auto lim_plus  = [&](const Vec& r) { return grid.mc3_limiter ? flux_lim_mc3_plus (r, beta) : flux_lim(r); };
    auto lim_minus = [&](const Vec& r) { return grid.mc3_limiter ? flux_lim_mc3_minus(r, beta) : flux_lim(r); };

    // extrapolated cell-edge variables (MUSCL paper eq 4.5). Uniform mesh keeps
    // the exact legacy 0.5 half-cell weight; the refined mesh uses the owning
    // cell's half-width factor (W1/W3/W4, = 0.5 when uniform).
    Vec Rxn_state_iph, Lxn_state_iph, Rxn_state_imh, Lxn_state_imh;
    if (grid.uniform_mesh) {
        Rxn_state_iph = prim_xn_state_ip1 - 0.5 * lim_minus(r_state_ip1) % (prim_xn_state_ip2 - prim_xn_state_ip1);
        Lxn_state_iph = prim_xn_state     + 0.5 * lim_plus (r_state)     % (prim_xn_state_ip1 - prim_xn_state);
        Rxn_state_imh = prim_xn_state     - 0.5 * lim_minus(r_state)     % (prim_xn_state_ip1 - prim_xn_state);
        Lxn_state_imh = prim_xn_state_im1 + 0.5 * lim_plus (r_state_im1) % (prim_xn_state     - prim_xn_state_im1);
    } else {
        Rxn_state_iph = prim_xn_state_ip1 - W3 % lim_minus(r_state_ip1) % (prim_xn_state_ip2 - prim_xn_state_ip1);
        Lxn_state_iph = prim_xn_state     + W1 % lim_plus (r_state)     % (prim_xn_state_ip1 - prim_xn_state);
        Rxn_state_imh = prim_xn_state     - W1 % lim_minus(r_state)     % (prim_xn_state_ip1 - prim_xn_state);
        Lxn_state_imh = prim_xn_state_im1 + W4 % lim_plus (r_state_im1) % (prim_xn_state     - prim_xn_state_im1);
    }

    if ((Lxn_state_iph.has_nan() + Rxn_state_iph.has_nan() +
         Lxn_state_imh.has_nan() + Rxn_state_imh.has_nan()) > 0) {
        std::cout << "Lxn_iph has nan: " << Lxn_state_iph.has_nan() << std::endl
                  << "Rxn_iph has nan: " << Rxn_state_iph.has_nan() << std::endl
                  << "Lxn_imh has nan: " << Lxn_state_imh.has_nan() << std::endl
                  << "Rxn_imh has nan: " << Rxn_state_imh.has_nan() << std::endl;
        std::exit(1);
    }

    // Undo the log transform on the reconstructed faces before mapping to conserved.
    if (grid.log_reconstruct) {
        prim_from_log(grid, Rxn_state_iph);
        prim_from_log(grid, Lxn_state_iph);
        prim_from_log(grid, Rxn_state_imh);
        prim_from_log(grid, Lxn_state_imh);
    }

    Rxn_state_iph = prim2cons(grid, Rxn_state_iph);
    Lxn_state_iph = prim2cons(grid, Lxn_state_iph);
    Rxn_state_imh = prim2cons(grid, Rxn_state_imh);
    Lxn_state_imh = prim2cons(grid, Lxn_state_imh);

    // ---- prediction step (half-step) ------------------------------------
    Vec prim_xt_state     = cons2prim(grid,
        xn_state - grid.dt_state / grid.ds_state %
        (cal_flux_state(grid, Lxn_state_iph) - cal_flux_state(grid, Rxn_state_imh)));
    // Positivity safeguard (flare scenario only): the half-step predictor can
    // undershoot density/pressure to ≤0 at the steep evaporation front on the
    // coarse grid, which then poisons the 2nd-order reconstruction and the
    // Rusanov spectral radius (c_s = √(γp/ρ) → NaN). Floor the predicted density
    // and pressure slots to small positive values. Gated on floors_active() so
    // untouched scenarios and the existing test suite are byte-for-byte unchanged
    // (a floor that only clips negatives never fires on those runs anyway).
    if (grid.floors_active()) {
        const auto  sz = arma::size(grid.ns, num_of_eq);
        const float RHO_FLOOR = grid.m_i * 1.0e10f;
        const float P_FLOOR   = 1.0e-8f;
        for (arma::uword i = 0; i < grid.ns; ++i) {
            float& ri = prim_xt_state(arma::sub2ind(sz, i, prim::RHO_I));
            float& rn = prim_xt_state(arma::sub2ind(sz, i, prim::RHO_N));
            float& pi = prim_xt_state(arma::sub2ind(sz, i, prim::P_I));
            float& pn = prim_xt_state(arma::sub2ind(sz, i, prim::P_N));
            if (ri < RHO_FLOOR) ri = RHO_FLOOR;
            if (rn < RHO_FLOOR) rn = RHO_FLOOR;
            if (pi < P_FLOOR)   pi = P_FLOOR;
            if (pn < P_FLOOR)   pn = P_FLOOR;
        }
    }

    Vec prim_xt_state_ip1 = cons2prim(grid, ip1(grid, prim2cons(grid, prim_xt_state)));
    Vec prim_xt_state_im1 = cons2prim(grid, im1(grid, prim2cons(grid, prim_xt_state)));
    // Same well-balanced φ_g correction for the predictor-shifted neighbours, which
    // feed the 2nd-order corrector below (otherwise they re-introduce the ⅓-of-truth
    // hydrostatic slope and a residual (γ−1)/2 g downforce).
    if (grid.well_balanced) {
        wb_correct_shifted(grid, prim_xt_state_ip1, phi_ip1, phi_cell);
        wb_correct_shifted(grid, prim_xt_state_im1, phi_im1, phi_cell);
    }

    // Log the predictor states too, so the 2nd-order blend below (which mixes them
    // with the still-logged prim_xn_state_* and reuses the log-space limiter ratios)
    // is consistent. After the positivity floor / wb correction, which need linear p.
    if (grid.log_reconstruct) {
        prim_to_log(grid, prim_xt_state);
        prim_to_log(grid, prim_xt_state_ip1);
        prim_to_log(grid, prim_xt_state_im1);
    }

    // ---- 2nd-order-in-time reconstruction -------------------------------
    // The 0.5·(xn+xt) prefactor is the predictor-corrector TIME average (mesh
    // independent); only the SPATIAL slope term carries the half-cell weight,
    // so on a refined mesh it becomes W1/W3/W4 (= 0.5 when uniform).
    if (grid.uniform_mesh) {
        Rxn_state_iph = 0.5 * (prim_xn_state_ip1 + prim_xt_state_ip1) - 0.5 * lim_minus(r_state_ip1) % (prim_xn_state_ip2 - prim_xn_state_ip1);
        Lxn_state_iph = 0.5 * (prim_xn_state     + prim_xt_state)     + 0.5 * lim_plus (r_state)     % (prim_xn_state_ip1 - prim_xn_state);
        Rxn_state_imh = 0.5 * (prim_xn_state     + prim_xt_state)     - 0.5 * lim_minus(r_state)     % (prim_xn_state_ip1 - prim_xn_state);
        Lxn_state_imh = 0.5 * (prim_xn_state_im1 + prim_xt_state_im1) + 0.5 * lim_plus (r_state_im1) % (prim_xn_state     - prim_xn_state_im1);
    } else {
        Rxn_state_iph = 0.5 * (prim_xn_state_ip1 + prim_xt_state_ip1) - W3 % lim_minus(r_state_ip1) % (prim_xn_state_ip2 - prim_xn_state_ip1);
        Lxn_state_iph = 0.5 * (prim_xn_state     + prim_xt_state)     + W1 % lim_plus (r_state)     % (prim_xn_state_ip1 - prim_xn_state);
        Rxn_state_imh = 0.5 * (prim_xn_state     + prim_xt_state)     - W1 % lim_minus(r_state)     % (prim_xn_state_ip1 - prim_xn_state);
        Lxn_state_imh = 0.5 * (prim_xn_state_im1 + prim_xt_state_im1) + W4 % lim_plus (r_state_im1) % (prim_xn_state     - prim_xn_state_im1);
    }

    // Undo the log transform before the positivity floor (linear thresholds) and
    // prim2cons. With log reconstruction the ρ/p faces are exp(·) > 0 by construction,
    // so the floor only ever fires on the linear V/U slots' companions — harmless.
    if (grid.log_reconstruct) {
        prim_from_log(grid, Rxn_state_iph);
        prim_from_log(grid, Lxn_state_iph);
        prim_from_log(grid, Rxn_state_imh);
        prim_from_log(grid, Lxn_state_imh);
    }

    // Positivity safeguard (flare scenario only): floor the reconstructed face
    // primitives' density/pressure before prim2cons, so the Rusanov spectral
    // radius c_s = √(γp/ρ) at the steep evaporation front stays real. Gated on
    // floors_active() → untouched scenarios / test suite unchanged.
    if (grid.floors_active()) {
        const auto  sz = arma::size(grid.ns, num_of_eq);
        const float RHO_FLOOR = grid.m_i * 1.0e10f;
        const float P_FLOOR   = 1.0e-8f;
        Vec* faces[4] = {&Rxn_state_iph, &Lxn_state_iph, &Rxn_state_imh, &Lxn_state_imh};
        for (Vec* fp : faces) {
            Vec& F = *fp;
            for (arma::uword i = 0; i < grid.ns; ++i) {
                float& ri = F(arma::sub2ind(sz, i, prim::RHO_I));
                float& rn = F(arma::sub2ind(sz, i, prim::RHO_N));
                float& pi = F(arma::sub2ind(sz, i, prim::P_I));
                float& pn = F(arma::sub2ind(sz, i, prim::P_N));
                if (ri < RHO_FLOOR) ri = RHO_FLOOR;
                if (rn < RHO_FLOOR) rn = RHO_FLOOR;
                if (pi < P_FLOOR)   pi = P_FLOOR;
                if (pn < P_FLOOR)   pn = P_FLOOR;
            }
        }
    }

    Rxn_state_iph = prim2cons(grid, Rxn_state_iph);
    Lxn_state_iph = prim2cons(grid, Lxn_state_iph);
    Rxn_state_imh = prim2cons(grid, Rxn_state_imh);
    Lxn_state_imh = prim2cons(grid, Lxn_state_imh);

    // ---- Rusanov / local Lax–Friedrichs flux ----------------------------
    const Vec a_state_imh = arma::max(cal_spectral_radius_state(grid, Lxn_state_imh),
                                      cal_spectral_radius_state(grid, Rxn_state_imh));
    const Vec a_state_iph = arma::max(cal_spectral_radius_state(grid, Lxn_state_iph),
                                      cal_spectral_radius_state(grid, Rxn_state_iph));

    const Vec F_state_imh = 0.5 * (cal_flux_state(grid, Lxn_state_imh) + cal_flux_state(grid, Rxn_state_imh)
                                   - a_state_imh % (Rxn_state_imh - Lxn_state_imh));
    const Vec F_state_iph = 0.5 * (cal_flux_state(grid, Lxn_state_iph) + cal_flux_state(grid, Rxn_state_iph)
                                   - a_state_iph % (Rxn_state_iph - Lxn_state_iph));

    if ((a_state_imh.has_nan() + a_state_iph.has_nan() +
         F_state_imh.has_nan() + F_state_iph.has_nan()) > 0) {
        std::cout << "a_imh has nan: " << a_state_imh.has_nan() << std::endl
                  << "a_iph has nan: " << a_state_iph.has_nan() << std::endl
                  << "Fximh has nan: " << F_state_imh.has_nan() << std::endl
                  << "Fxiph has nan: " << F_state_iph.has_nan() << std::endl;
        std::exit(1);
    }

    Vec rhs = -1.0 * grid.B_state % (F_state_iph / grid.B_state_iph - F_state_imh / grid.B_state_imh) / grid.ds_state
              + cal_source_state(grid, xn_state);

    // Equilibrium-reference well-balancing (chromosphere.hpp::eq_wb): subtract the
    // frozen scheme residual at the reference equilibrium so eq_state is an exact
    // discrete steady state. The empty() guard means the one-time computation of
    // eq_residual itself (integrators::ensure_eq_residual, which calls this on
    // eq_state while eq_residual is still empty) is NOT self-corrected, and the
    // default-off path (eq_residual never populated) is byte-identical.
    if (grid.eq_wb && !grid.eq_residual.is_empty()) {
        rhs -= grid.eq_residual;
    }
    return rhs;
}

// ============================================================================
// Implicit RHS (writeup §4.3, R_I): ion-neutral drag, collisional + frictional
// heating, conservative field-aligned heat conduction. Pressure-area and
// gravity live in cal_source_state (R_E).
// ============================================================================

Vec rhs_implicit_state(const Grid& grid, const Vec& xn_state) {
    Vec RI_state(arma::size(xn_state), arma::fill::zeros);

    // ---- cell-centered primitives ---------------------------------------
    const Vec rho_i  = get_scalar(grid, xn_state, cons::RHO_I);
    const Vec rho_n  = get_scalar(grid, xn_state, cons::RHO_N);
    const Vec rhoV_i = get_scalar(grid, xn_state, cons::MOM_I);
    const Vec rhoU_n = get_scalar(grid, xn_state, cons::MOM_N);
    const Vec e_i    = get_scalar(grid, xn_state, cons::E_I);
    const Vec e_n    = get_scalar(grid, xn_state, cons::E_N);

    const Vec n_i   = rho_i / grid.m_i;
    const Vec n_n   = rho_n / grid.m_n;
    const Vec V     = rhoV_i / rho_i;
    const Vec U     = rhoU_n / rho_n;
    const Vec phi_g = 0.5 * (grid.phi_g_imh + grid.phi_g_iph);
    const Vec p_i   = grid.gm1() * e_i - grid.half_gm1() * rho_i % V % V - grid.gm1() * rho_i % phi_g;
    const Vec p_n   = grid.gm1() * e_n - grid.half_gm1() * rho_n % U % U - grid.gm1() * rho_n % phi_g;
    const Vec T_i   = p_i / (2.0 * n_i * grid.k_b);
    const Vec T_n   = p_n / (n_n * grid.k_b);

    const Vec Ke = kappa_e(n_i, n_n, T_i);
    const Vec Ki(arma::size(T_i), arma::fill::zeros);
    const Vec Kn = kappa_n(n_i, n_n, T_i, T_n);

    // ---- i+1 neighbour (uses outer ghost at i = ns-1) -------------------
    const Vec xn_state_ip1 = ip1(grid, xn_state);
    const Vec rho_i_ip1  = get_scalar(grid, xn_state_ip1, cons::RHO_I);
    const Vec rho_n_ip1  = get_scalar(grid, xn_state_ip1, cons::RHO_N);
    const Vec rhoV_i_ip1 = get_scalar(grid, xn_state_ip1, cons::MOM_I);
    const Vec rhoU_n_ip1 = get_scalar(grid, xn_state_ip1, cons::MOM_N);
    const Vec e_i_ip1    = get_scalar(grid, xn_state_ip1, cons::E_I);
    const Vec e_n_ip1    = get_scalar(grid, xn_state_ip1, cons::E_N);

    const Vec n_i_ip1 = rho_i_ip1 / grid.m_i;
    const Vec n_n_ip1 = rho_n_ip1 / grid.m_n;
    const Vec V_ip1   = rhoV_i_ip1 / rho_i_ip1;
    const Vec U_ip1   = rhoU_n_ip1 / rho_n_ip1;
    Vec phi_g_ip1 = ip1(grid, 0.5 * (grid.phi_g_iph + grid.phi_g_imh), SLICE);
    phi_g_ip1[grid.ns - 1] = grid.phi_g_iph[grid.ns - 1];
    const Vec p_i_ip1 = grid.gm1() * e_i_ip1 - grid.half_gm1() * rho_i_ip1 % V_ip1 % V_ip1 - grid.gm1() * rho_i_ip1 % phi_g_ip1;
    const Vec p_n_ip1 = grid.gm1() * e_n_ip1 - grid.half_gm1() * rho_n_ip1 % U_ip1 % U_ip1 - grid.gm1() * rho_n_ip1 % phi_g_ip1;
    const Vec T_i_ip1 = p_i_ip1 / (2.0 * n_i_ip1 * grid.k_b);
    const Vec T_n_ip1 = p_n_ip1 / (n_n_ip1 * grid.k_b);

    const Vec Ke_ip1 = kappa_e(n_i_ip1, n_n_ip1, T_i_ip1);
    const Vec Ki_ip1(arma::size(T_i), arma::fill::zeros);
    const Vec Kn_ip1 = kappa_n(n_i_ip1, n_n_ip1, T_i_ip1, T_n_ip1);

    // ---- i-1 neighbour (uses inner ghost at i = 0) ----------------------
    const Vec xn_state_im1 = im1(grid, xn_state);
    const Vec rho_i_im1  = get_scalar(grid, xn_state_im1, cons::RHO_I);
    const Vec rho_n_im1  = get_scalar(grid, xn_state_im1, cons::RHO_N);
    const Vec rhoV_i_im1 = get_scalar(grid, xn_state_im1, cons::MOM_I);
    const Vec rhoU_n_im1 = get_scalar(grid, xn_state_im1, cons::MOM_N);
    const Vec e_i_im1    = get_scalar(grid, xn_state_im1, cons::E_I);
    const Vec e_n_im1    = get_scalar(grid, xn_state_im1, cons::E_N);

    const Vec n_i_im1 = rho_i_im1 / grid.m_i;
    const Vec n_n_im1 = rho_n_im1 / grid.m_n;
    const Vec V_im1   = rhoV_i_im1 / rho_i_im1;
    const Vec U_im1   = rhoU_n_im1 / rho_n_im1;
    Vec phi_g_im1 = im1(grid, 0.5 * (grid.phi_g_iph + grid.phi_g_imh), SLICE);
    phi_g_im1[0] = grid.phi_g_imh[0];
    const Vec p_i_im1 = grid.gm1() * e_i_im1 - grid.half_gm1() * rho_i_im1 % V_im1 % V_im1 - grid.gm1() * rho_i_im1 % phi_g_im1;
    const Vec p_n_im1 = grid.gm1() * e_n_im1 - grid.half_gm1() * rho_n_im1 % U_im1 % U_im1 - grid.gm1() * rho_n_im1 % phi_g_im1;
    const Vec T_i_im1 = p_i_im1 / (2.0 * n_i_im1 * grid.k_b);
    const Vec T_n_im1 = p_n_im1 / (n_n_im1 * grid.k_b);

    const Vec Ke_im1 = kappa_e(n_i_im1, n_n_im1, T_i_im1);
    const Vec Ki_im1(arma::size(T_i), arma::fill::zeros);
    const Vec Kn_im1 = kappa_n(n_i_im1, n_n_im1, T_i_im1, T_n_im1);

    // ---- ion-neutral drag (α = ρ_i ν_in) --------------------------------
    const Vec nu_in_v = nu_in(grid, n_n, T_i, T_n);
    const Vec alpha   = rho_i % nu_in_v;
    const Vec w       = V - U;
    RI_state += scalar_to(grid, -alpha % w, cons::MOM_I);
    RI_state += scalar_to(grid,  alpha % w, cons::MOM_N);

    // ---- collisional + frictional heating (R_I energy rows) -------------
    RI_state += scalar_to(grid,
        alpha / (grid.m_i + grid.m_n) % (3.0 * grid.k_b * (T_n - T_i) + grid.m_n * w % w) - alpha % V % w,
        cons::E_I);
    RI_state += scalar_to(grid,
        alpha / (grid.m_i + grid.m_n) % (3.0 * grid.k_b * (T_i - T_n) + grid.m_i * w % w) + alpha % U % w,
        cons::E_N);

    // ---- conservative field-aligned heat conduction (writeup §4.4) ----
    //   C = (B_i/Δs_i) [ K_{i+1/2}/B_{i+1/2} · (T_{i+1}-T_i)/Δs_{i+1/2}
    //                  - K_{i-1/2}/B_{i-1/2} · (T_i-T_{i-1})/Δs_{i-1/2} ]
    // Face spacings use Neumann BC at both ends (mirror Δs_i across the
    // outermost/innermost interior cell).
    const Vec ds_i_ip1 = ip1(grid, grid.ds_i, SLICE);
    const Vec ds_i_im1 = im1(grid, grid.ds_i, SLICE);
    const Vec ds_iph   = 0.5 * (grid.ds_i + ds_i_ip1);
    const Vec ds_imh   = 0.5 * (ds_i_im1 + grid.ds_i);

    // Face conductivities. Uniform mesh keeps the legacy arithmetic average;
    // a refined mesh uses the width-weighted series-resistance combination
    // (physics.hpp::face_conductivity_series) so unequal half-cells conduct
    // correctly. Ki ≡ 0, so (Ke+Ki) is the charged-row conductivity.
    Vec Kei_iph, Kei_imh, Kn_iph, Kn_imh;
    if (grid.uniform_mesh) {
        Kei_iph = 0.5 * ((Ke + Ki) + (Ke_ip1 + Ki_ip1));
        Kei_imh = 0.5 * ((Ke + Ki) + (Ke_im1 + Ki_im1));
        Kn_iph  = 0.5 * (Kn + Kn_ip1);
        Kn_imh  = 0.5 * (Kn + Kn_im1);
    } else {
        Kei_iph = face_conductivity_series(Ke + Ki, Ke_ip1 + Ki_ip1, grid.ds_i, ds_i_ip1);
        Kei_imh = face_conductivity_series(Ke + Ki, Ke_im1 + Ki_im1, grid.ds_i, ds_i_im1);
        Kn_iph  = face_conductivity_series(Kn, Kn_ip1, grid.ds_i, ds_i_ip1);
        Kn_imh  = face_conductivity_series(Kn, Kn_im1, grid.ds_i, ds_i_im1);
    }

    // Same face-local artificial conduction as the production implicit solver:
    // chi_num,f = C_num*ds_face and K_num,f = chi_num,f*C_V,f. The charged
    // single-temperature capacity is 3 n_i k_B; the neutral capacity is
    // 1.5 n_n k_B, both face-averaged with the existing ghost geometry.
    const Vec chi_iph = numerical_diffusivity_at_face(grid, ds_iph);
    const Vec chi_imh = numerical_diffusivity_at_face(grid, ds_imh);
    Kei_iph += (1.5f*grid.k_b)*(n_i+n_i_ip1)%chi_iph;
    Kei_imh += (1.5f*grid.k_b)*(n_i+n_i_im1)%chi_imh;
    Kn_iph  += (0.75f*grid.k_b)*(n_n+n_n_ip1)%chi_iph;
    Kn_imh  += (0.75f*grid.k_b)*(n_n+n_n_im1)%chi_imh;

    const Vec q_i_iph = Kei_iph / grid.B_iph % (T_i_ip1 - T_i)     / ds_iph;
    const Vec q_i_imh = Kei_imh / grid.B_imh % (T_i     - T_i_im1) / ds_imh;
    const Vec q_n_iph = Kn_iph  / grid.B_iph % (T_n_ip1 - T_n)     / ds_iph;
    const Vec q_n_imh = Kn_imh  / grid.B_imh % (T_n     - T_n_im1) / ds_imh;

    const Vec C_i = grid.B_i % (q_i_iph - q_i_imh) / grid.ds_i;
    const Vec C_n = grid.B_i % (q_n_iph - q_n_imh) / grid.ds_i;
    RI_state += scalar_to(grid, C_i, cons::E_I);
    RI_state += scalar_to(grid, C_n, cons::E_N);

    return RI_state;
}

} // namespace chromosphere
