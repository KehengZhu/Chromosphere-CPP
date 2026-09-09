/*!
 * @file two_temp/conduction.cpp
 * @brief EXPERIMENTAL two-temperature solver: coupled implicit conduction on
 *        (T_e, T_i) with the collisional exchange, and the boundary conditions
 *        the characteristic budget of the 2T system implies.
 * @ingroup two_temp_solver
 *
 * The second and last stage of the two-temperature timestep, replacing the
 * release's single scalar conduction solve with a coupled one:
 *
 *   E_e(rho, T_e) - E_e*  = dt [ div(kappa_e grad T_e) + g_ex (T_i - T_e) ]
 *   e_h(rho, T_i) - e_h*  = dt [ div(kappa_n grad T_i) - g_ex (T_i - T_e) ]
 *
 * solved together by quasi-Newton on the cell temperatures with a 2x2-BLOCK
 * tridiagonal Jacobian: the two conduction operators never couple across cells
 * (each channel conducts on its own temperature), and the only coupling is the
 * pointwise collisional exchange, which lands entirely in the diagonal block.
 * kappa and g_ex are lagged, exactly as the release conduction solve lags
 * d(kappa)/dT.
 *
 * The exchange is ANTISYMMETRIC by construction and both conductive divergences
 * are re-summed into the total-energy row, so the stage conserves total energy
 * exactly and conserves mass and momentum trivially (it never touches them).
 *
 * ## Why this stage must be implicit
 *
 * The measured electron-heavy equilibration time in this column is 0.092 ms at
 * the 1600 km base, 0.163 ms mid-domain and 3.92 ms at the 2153 km top
 * (`util/tau_ei_equilibration.py`), against a MEASURED CFL-0.50 step of
 * 5.53-5.62 ms on this configuration. So tau_eq/dt is about 0.016 at the base,
 * 0.029 mid-domain and 0.70 in the fully ionized top cells: the exchange is
 * stiff by one to two orders of magnitude through the bulk, and in the top cells
 * it is only MARGINALLY resolved rather than stiff. (An earlier version of this
 * comment said "tens of milliseconds" and "three to four orders of magnitude
 * everywhere"; both were wrong.) Treating it implicitly is still clearly right —
 * an explicit exchange would need dt <~ 90 us, about 60x smaller — and it is
 * what keeps the acoustic CFL condition the only timestep constraint, the same
 * property the release relies on for its conduction stage. The top-cell result
 * is not merely a caveat: those are exactly the cells where a real decoupling is
 * physically plausible, which is where this study looks first.
 *
 * ## Boundary conditions, from the characteristics
 *
 * The 4x4 hyperbolic system in (rho, rho u, E, E_e) has eigenvalues
 * {u - c, u, u, u + c}: the acoustic pair, and a DOUBLE eigenvalue at u carrying
 * the total entropy and the electron/heavy energy partition. Both temperatures
 * therefore ride on lambda = u, not on an acoustic characteristic. At the outer
 * face this column is a slow subsonic OUTFLOW (0 < u << c), so three
 * characteristics leave and exactly one (u - c) enters: exactly ONE condition may
 * be imposed from outside, and the release already spends it on the reservoir
 * back-pressure. Both temperatures must therefore be EXTRAPOLATED by the
 * hydrodynamic boundary closure, never imposed — which is what
 * `model_column_2t` does. (Were the top flow to reverse to subsonic inflow, three
 * characteristics would enter and both temperatures would have to be supplied
 * with the pressure; the Mach-capped outflow is what keeps that from happening.)
 *
 * The CONDUCTION operator is parabolic, not hyperbolic, and its boundary data is
 * a separate budget: one condition per temperature per face. Here the physics
 * chooses them:
 *
 *   * OUTER, electron: DIRICHLET at the imposed reservoir temperature. The
 *     downward conductive flux from the transition region and corona is carried
 *     by ELECTRONS — above the domain top the hydrogen is fully ionized, so
 *     kappa_n collapses (it is proportional to n_HI) while kappa_e is the Spitzer
 *     T^{5/2} channel. Assigning the release's 22 000 K wall to the electron
 *     channel is what makes the heating chain explicit:
 *     TR conduction -> electrons -> collisional exchange -> heavies -> pressure
 *     -> evaporation.
 *   * OUTER, heavy: NEUMANN, zero conductive flux. There is no heavy-particle
 *     conductive flux to import from a fully ionized corona, and Spitzer PROTON
 *     conduction, which we do not model, is only ~2.3 % of kappa_e. Imposing a
 *     heavy Dirichlet datum instead would manufacture a spurious neutral-channel
 *     flux at a face where n_HI has already gone to zero. `TT_OUTER_TI=dirichlet`
 *     selects that alternative for the controlled comparison.
 *   * INNER, both: DIRICHLET at the truncation reservoir temperature. The 1600 km
 *     reservoir is dense and cool, tau_eq there is 0.09 ms, so T_e = T_i is a
 *     statement about the reservoir rather than an extra assumption.
 *     `Grid::inner_conduction_neumann` makes the base insulating for BOTH
 *     channels, keeping the release knob's meaning.
 */

#include "two_temp/two_temp.hpp"
#include "profiling.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <string>

namespace chromosphere {

namespace {

inline arma::SizeMat state_size(const Grid& grid) {
    return arma::size(grid.ns, num_of_two_temp_eq);
}

/// Thermal state of one boundary layer as the CONDUCTION rows see it.
struct ConductionWall {
    double T_e, T_i, n_e, n_HI;
};

/// Decode a two-temperature ghost buffer into its wall thermal state.
ConductionWall wall_from_ghost(const Grid& grid, const Vec& ghost, double phi) {
    const TwoTempThermo th = two_temp_decode(
        grid.eos_gamma_table, ghost(tt::RHO), ghost(tt::MOM), ghost(tt::ENERGY),
        ghost(tt::E_ELEC), phi, std::numeric_limits<double>::quiet_NaN());
    return {th.T_e, th.T_i, th.n_e, th.n_HI};
}

/// The OUTER wall with the imposed conduction temperature substituted in. The
/// material state is rebuilt at the imposed face temperature from the externally
/// imposed back-pressure the ghost carries, the same construction the release
/// conduction stage uses. At 22 000 K the gas is fully ionized and T_e = T_i
/// there, so the equilibrium pressure inverse is exact at that state.
///
/// The imposed temperature always enters the ELECTRON slot: the conductive
/// reservoir above the domain is electron-conducted. What goes into the HEAVY
/// slot depends on which comparison is running, and the distinction is the whole
/// point of `TT_OUTER_TI`:
///
///   * `tt_outer_Ti_neumann` (default) — the heavy channel is insulating, so its
///     outer face conductance is zeroed by the caller and this slot is never
///     read. It is left at the ghost value rather than at the wall so that a
///     stray read would carry the physically meaningful extrapolated value.
///   * otherwise — the CONTROLLED COMPARISON deliberately imposes the SAME wall
///     temperature on the heavy channel that the electron channel gets. That is
///     the point of the test: it asks what the naive closure, which lets the
///     heavy species feel the coronal reservoir directly, would have produced.
///     Handing the heavy channel the extrapolated ghost temperature instead
///     would make this leg a near-copy of the Neumann leg and the comparison
///     would be worthless.
ConductionWall outer_conduction_wall(const Grid& grid,
                                     const ConductionWall& ghost,
                                     double ghost_pressure) {
    if (!grid.outer_conduction_temperature_override) return ghost;
    const double T = static_cast<double>(grid.outer_conduction_temperature);
    const double rho = equilibrium_density_from_pressure(ghost_pressure, T);
    const double n_h = rho/eos_constants::m_h;
    const double x = saha_ionization_fraction_n_h(n_h, T);
    const double T_i_wall = grid.tt_outer_Ti_neumann ? ghost.T_i : T;
    return {T, T_i_wall, x*n_h, (1.0 - x)*n_h};
}

/// Block-tridiagonal Thomas algorithm with 2x2 blocks. `a` and `c` hold the two
/// DIAGONAL entries of the sub/super blocks (the two conduction channels never
/// couple across cells); `b` holds the full 2x2 diagonal block row-major.
/// Destroys `b` and `d`.
void block_thomas_solve(std::vector<double>& a, std::vector<double>& b,
                        std::vector<double>& c, std::vector<double>& d,
                        std::vector<double>& x, std::size_t n) {
    auto invert = [](const double* m, double* out) {
        const double det = m[0]*m[3] - m[1]*m[2];
        if (!(std::abs(det) > 0.0) || !std::isfinite(det))
            throw std::runtime_error(
                "two-temperature conduction Jacobian block is singular");
        const double inv = 1.0/det;
        out[0] =  m[3]*inv; out[1] = -m[1]*inv;
        out[2] = -m[2]*inv; out[3] =  m[0]*inv;
    };
    for (std::size_t i = 1; i < n; ++i) {
        double b_inv[4];
        invert(&b[4*(i-1)], b_inv);
        // m = A_i * B_{i-1}^{-1}, with A_i diagonal.
        const double m[4] = {a[2*i]*b_inv[0],     a[2*i]*b_inv[1],
                             a[2*i+1]*b_inv[2],   a[2*i+1]*b_inv[3]};
        // B_i -= m * C_{i-1}, with C_{i-1} diagonal.
        b[4*i+0] -= m[0]*c[2*(i-1)];
        b[4*i+1] -= m[1]*c[2*(i-1)+1];
        b[4*i+2] -= m[2]*c[2*(i-1)];
        b[4*i+3] -= m[3]*c[2*(i-1)+1];
        // d_i -= m * d_{i-1}
        const double d0 = d[2*(i-1)], d1 = d[2*(i-1)+1];
        d[2*i]   -= m[0]*d0 + m[1]*d1;
        d[2*i+1] -= m[2]*d0 + m[3]*d1;
    }
    {
        double b_inv[4];
        invert(&b[4*(n-1)], b_inv);
        x[2*(n-1)]   = b_inv[0]*d[2*(n-1)] + b_inv[1]*d[2*(n-1)+1];
        x[2*(n-1)+1] = b_inv[2]*d[2*(n-1)] + b_inv[3]*d[2*(n-1)+1];
    }
    for (std::size_t i = n - 1; i > 0; --i) {
        const std::size_t j = i - 1;
        double b_inv[4];
        invert(&b[4*j], b_inv);
        const double r0 = d[2*j]     - c[2*j]*x[2*i];
        const double r1 = d[2*j+1]   - c[2*j+1]*x[2*i+1];
        x[2*j]   = b_inv[0]*r0 + b_inv[1]*r1;
        x[2*j+1] = b_inv[2]*r0 + b_inv[3]*r1;
    }
}

} // namespace

Vec two_temp_apply_conduction(const Grid& grid, const Vec& state, double dt,
                              TwoTempWorkspace& s) {
    if (!(dt > 0.0)) return state;
    ProfileScope timer(ProfileRegion::Conduction);
    const arma::uword ns = grid.ns;
    const std::size_t n = ns;
    const auto sz = state_size(grid);
    const double t_min = grid.eos_gamma_table.min_temperature();
    const double t_max = grid.eos_gamma_table.max_temperature();

    // ---- seed the solve from the post-hydro state -------------------------
    for (arma::uword i = 0; i < ns; ++i) {
        auto at = [&](arma::uword row) {
            return static_cast<double>(state(arma::sub2ind(sz, i, row)));
        };
        const TwoTempThermo th = two_temp_decode(
            grid.eos_gamma_table, at(tt::RHO), at(tt::MOM), at(tt::ENERGY),
            at(tt::E_ELEC), two_temp_cell_phi(grid, i),
            i < s.c_Te.size() && s.c_Te[i] > 0.0
                ? s.c_Te[i] : std::numeric_limits<double>::quiet_NaN());
        s.c_rho[i] = th.rho;
        s.c_nH[i] = th.n_H;
        s.c_Te[i] = th.T_e;
        s.c_Ti[i] = th.T_i;
        s.c_Ee_old[i] = th.e_elec;
        s.c_eh_old[i] = th.e_heavy;
        // The heavy heat capacity is exactly (3/2) n_H k_B, independent of T_i:
        // the heavy pool carries no ionization energy.
        s.c_cap_i[i] = 1.5*th.n_H*eos_constants::k_b;
    }

    // ---- boundary thermal data --------------------------------------------
    const ConductionWall inner = wall_from_ghost(
        grid, grid.tt_inner_boundary0, grid.phi_g_imh(0));
    const ConductionWall outer_ghost = wall_from_ghost(
        grid, grid.tt_outer_boundary0, grid.phi_g_iph(ns - 1));
    double outer_ghost_pressure = 0.0;
    {
        const TwoTempThermo th = two_temp_decode(
            grid.eos_gamma_table, grid.tt_outer_boundary0(tt::RHO),
            grid.tt_outer_boundary0(tt::MOM), grid.tt_outer_boundary0(tt::ENERGY),
            grid.tt_outer_boundary0(tt::E_ELEC), grid.phi_g_iph(ns - 1),
            std::numeric_limits<double>::quiet_NaN());
        outer_ghost_pressure = th.p;
    }
    const ConductionWall outer =
        outer_conduction_wall(grid, outer_ghost, outer_ghost_pressure);

    const double k_inner_e = two_temp_kappa_e(inner.n_e, inner.n_HI, inner.T_e);
    const double k_inner_i = two_temp_kappa_i(inner.n_e, inner.n_HI, inner.T_i);
    const double k_outer_e = two_temp_kappa_e(outer.n_e, outer.n_HI, outer.T_e);
    const double k_outer_i = two_temp_kappa_i(outer.n_e, outer.n_HI, outer.T_i);

    auto face_k = [&](double kh, double kt, double dh, double dtw) {
        if (grid.uniform_mesh) return 0.5*(kh + kt);
        kh = std::max(kh, 1.0e-30);
        kt = std::max(kt, 1.0e-30);
        return (dh + dtw)/(dh/kh + dtw/kt);
    };
    // The imposed outer electron temperature is a PHYSICAL-FACE datum, so its
    // conductance uses the half top cell and the wall conductivity directly —
    // the same relocation the release conduction stage performs.
    const bool physical_outer_face = grid.outer_conduction_temperature_override;
    const double ds_right_top = physical_outer_face
        ? 0.5*static_cast<double>(grid.ds_i(ns - 1))
        : static_cast<double>(grid.ds_iph_i(ns - 1));

    bool converged = false;
    int passes = 0;
    for (int iteration = 0; iteration <= 40; ++iteration) {
        // --- lagged coefficients at the current iterate --------------------
        for (arma::uword i = 0; i < ns; ++i) {
            const CaloricState cs = equilibrium_caloric_state(s.c_rho[i], s.c_Te[i]);
            s.c_ne[i] = cs.n_e;
            s.c_nHI[i] = cs.n_hi;
            // E_e at the iterate temperature, and its derivative, as exact
            // differences of the release caloric closure.
            s.c_Ee_at_T[i] = cs.internal_energy
                - 1.5*s.c_nH[i]*eos_constants::k_b*s.c_Te[i];
            s.c_cap_e[i] = cs.heat_capacity - 1.5*s.c_nH[i]*eos_constants::k_b;
            s.c_ke[i] = two_temp_kappa_e(cs.n_e, cs.n_hi, s.c_Te[i]);
            s.c_ki[i] = two_temp_kappa_i(cs.n_e, cs.n_hi, s.c_Ti[i]);
            s.c_gex[i] = two_temp_exchange_conductance(cs.n_e, cs.n_hi, s.c_Te[i]);
        }

        // --- face conductances --------------------------------------------
        for (arma::uword i = 0; i < ns; ++i) {
            const double ds_i = grid.ds_i(i);
            const double kl_e = i == 0
                ? face_k(s.c_ke[i], k_inner_e, ds_i, ds_i)
                : face_k(s.c_ke[i], s.c_ke[i-1], ds_i, grid.ds_i(i-1));
            const double kl_i = i == 0
                ? face_k(s.c_ki[i], k_inner_i, ds_i, ds_i)
                : face_k(s.c_ki[i], s.c_ki[i-1], ds_i, grid.ds_i(i-1));
            double kr_e, kr_i, ds_right;
            if (i + 1 == ns) {
                kr_e = physical_outer_face
                    ? k_outer_e : face_k(s.c_ke[i], k_outer_e, ds_i, ds_i);
                kr_i = physical_outer_face
                    ? k_outer_i : face_k(s.c_ki[i], k_outer_i, ds_i, ds_i);
                ds_right = ds_right_top;
            } else {
                kr_e = face_k(s.c_ke[i], s.c_ke[i+1], ds_i, grid.ds_i(i+1));
                kr_i = face_k(s.c_ki[i], s.c_ki[i+1], ds_i, grid.ds_i(i+1));
                ds_right = grid.ds_iph_i(i);
            }
            const double geo_l = grid.B_i(i)/ds_i/grid.B_imh(i)/grid.ds_imh_i(i);
            const double geo_r = grid.B_i(i)/ds_i/grid.B_iph(i)/ds_right;
            s.c_gE_L[i] = kl_e*geo_l;
            s.c_gI_L[i] = kl_i*geo_l;
            s.c_gE_R[i] = kr_e*geo_r;
            s.c_gI_R[i] = kr_i*geo_r;
            // Inner Neumann: drop the base coupling for BOTH channels.
            if (i == 0 && grid.inner_conduction_neumann) {
                s.c_gE_L[i] = 0.0;
                s.c_gI_L[i] = 0.0;
            }
            // Outer heavy Neumann: no heavy conductive flux is imported from the
            // fully ionized plasma above the domain (the default; see the file
            // header for the alternative).
            if (i + 1 == ns && grid.tt_outer_Ti_neumann) s.c_gI_R[i] = 0.0;
        }

        // --- residual and Jacobian ----------------------------------------
        double residual_max = 0.0;
        arma::uword residual_cell = 0;
        double residual_e_worst = 0.0, residual_i_worst = 0.0;
        for (arma::uword i = 0; i < ns; ++i) {
            const double te_l = i == 0 ? inner.T_e : s.c_Te[i-1];
            const double ti_l = i == 0 ? inner.T_i : s.c_Ti[i-1];
            const double te_r = i + 1 == ns ? outer.T_e : s.c_Te[i+1];
            const double ti_r = i + 1 == ns ? outer.T_i : s.c_Ti[i+1];
            const double div_e = s.c_gE_R[i]*(te_r - s.c_Te[i])
                               - s.c_gE_L[i]*(s.c_Te[i] - te_l);
            const double div_i = s.c_gI_R[i]*(ti_r - s.c_Ti[i])
                               - s.c_gI_L[i]*(s.c_Ti[i] - ti_l);
            const double exchange = s.c_gex[i]*(s.c_Ti[i] - s.c_Te[i]);
            s.c_Ee_target[i] = s.c_Ee_old[i] + dt*(div_e + exchange);
            s.c_eh_target[i] = s.c_eh_old[i] + dt*(div_i - exchange);
            const double r_e = s.c_Ee_at_T[i] - s.c_Ee_target[i];
            const double r_i = s.c_cap_i[i]*s.c_Ti[i] - s.c_eh_target[i];
            // Both residuals are normalized on the same energy scale (the heavy
            // pool, which is the larger of the two everywhere in this column), so
            // the tolerance means the same thing for both rows.
            const double scale = std::max(s.c_eh_old[i], 1.0e-30);
            const double scaled = std::max(std::abs(r_e), std::abs(r_i))/scale;
            if (scaled > residual_max) {
                residual_max = scaled;
                residual_cell = i;
                residual_e_worst = r_e/scale;
                residual_i_worst = r_i/scale;
            }

            s.blk_a[2*i]   = -dt*s.c_gE_L[i];
            s.blk_a[2*i+1] = -dt*s.c_gI_L[i];
            s.blk_c[2*i]   = -dt*s.c_gE_R[i];
            s.blk_c[2*i+1] = -dt*s.c_gI_R[i];
            s.blk_b[4*i+0] = s.c_cap_e[i]
                + dt*(s.c_gE_L[i] + s.c_gE_R[i] + s.c_gex[i]);
            s.blk_b[4*i+1] = -dt*s.c_gex[i];
            s.blk_b[4*i+2] = -dt*s.c_gex[i];
            s.blk_b[4*i+3] = s.c_cap_i[i]
                + dt*(s.c_gI_L[i] + s.c_gI_R[i] + s.c_gex[i]);
            s.blk_rhs[2*i]   = -r_e;
            s.blk_rhs[2*i+1] = -r_i;
        }
        s.blk_a[0] = 0.0; s.blk_a[1] = 0.0;
        s.blk_c[2*(n-1)] = 0.0; s.blk_c[2*(n-1)+1] = 0.0;

        if (residual_max < 2.0e-11) { converged = true; passes = iteration; break; }
        if (iteration == 40) break;

        block_thomas_solve(s.blk_a, s.blk_b, s.blk_c, s.blk_rhs, s.blk_delta, n);
        double step_max = 0.0;
        for (arma::uword i = 0; i < ns; ++i) {
            auto accept = [&](double current, double delta) {
                double candidate = current + delta;
                if (!std::isfinite(candidate))
                    throw std::runtime_error(
                        "two-temperature conduction Newton produced a non-finite T");
                if (candidate <= t_min) candidate = 0.5*(current + t_min);
                if (candidate >= t_max) candidate = 0.5*(current + t_max);
                step_max = std::max(step_max,
                    std::abs(candidate - current)/std::max(current, 1.0));
                return candidate;
            };
            s.c_Te[i] = accept(s.c_Te[i], s.blk_delta[2*i]);
            s.c_Ti[i] = accept(s.c_Ti[i], s.blk_delta[2*i+1]);
        }
        // A vanishing step with a live residual means the iteration has stalled,
        // which must be an error rather than a silently accepted state. The
        // threshold sits two decades BELOW the residual tolerance on purpose: for
        // the heavy row C_i T_i is identically e_h, so a relative residual and a
        // relative Newton step are the SAME number, and equal thresholds would
        // make the stall detector fire on an ordinary final Newton step.
        if (step_max < 1.0e-13 && iteration > 0) {
            std::ostringstream message;
            message.setf(std::ios::scientific);
            message.precision(6);
            message << "two-temperature conduction stagnated before residual "
                       "convergence: scaled_residual=" << residual_max
                    << " (cell " << residual_cell
                    << ", r_e=" << residual_e_worst
                    << ", r_i=" << residual_i_worst
                    << ", T_e=" << s.c_Te[residual_cell]
                    << ", T_i=" << s.c_Ti[residual_cell]
                    << "), relative_step=" << step_max
                    << ", iteration=" << iteration;
            throw std::runtime_error(message.str());
        }
    }
    if (!converged)
        throw std::runtime_error(
            "two-temperature conduction Newton did not converge");
    s.stats.max_conduction_iterations = std::max(
        s.stats.max_conduction_iterations, static_cast<std::uint64_t>(passes));

    // ---- write back -------------------------------------------------------
    // The ACCEPTED energies are the backward-Euler targets, not an EOS
    // round-trip of the converged temperatures, so the discrete energy balance
    // holds exactly and the exchange cancels between the two pools.
    Vec updated = state;
    for (arma::uword i = 0; i < ns; ++i) {
        const double rho = s.c_rho[i];
        const double momentum =
            static_cast<double>(state(arma::sub2ind(sz, i, tt::MOM)));
        const double phi = two_temp_cell_phi(grid, i);
        double e_elec = s.c_Ee_target[i];
        double e_heavy = s.c_eh_target[i];
        if (!std::isfinite(e_elec) || !std::isfinite(e_heavy))
            throw std::domain_error(
                "non-finite two-temperature conduction target at cell "
                + std::to_string(i) + ", time=" + std::to_string(grid.sim_time));
        const double e_elec_min = two_temp_electron_energy(rho, t_min);
        const double e_elec_max = two_temp_electron_energy(rho, t_max);
        if (e_elec < e_elec_min) { e_elec = e_elec_min; ++s.stats.decode_Te_clamps; }
        if (e_elec > e_elec_max) { e_elec = e_elec_max; ++s.stats.decode_Te_clamps; }
        const double e_heavy_min = 1.5*s.c_nH[i]*eos_constants::k_b*t_min;
        if (e_heavy < e_heavy_min) {
            e_heavy = e_heavy_min;
            ++s.stats.decode_Ti_clamps;
        }
        updated(arma::sub2ind(sz, i, tt::E_ELEC)) = static_cast<Real>(e_elec);
        updated(arma::sub2ind(sz, i, tt::ENERGY)) = static_cast<Real>(
            two_temp_total_energy(rho, momentum, e_heavy, e_elec, phi));
    }
    return updated;
}

} // namespace chromosphere
