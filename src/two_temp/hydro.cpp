/*!
 * @file two_temp/hydro.cpp
 * @brief EXPERIMENTAL two-temperature solver: MUSCL-Hancock reconstruction, the
 *        SWMF exact-Riemann Godunov flux on the three release rows, contact
 *        upwinding of the electron energy, the sources, and the CFL timestep.
 * @ingroup two_temp_solver
 *
 * One explicit RHS evaluation, in order:
 *
 *  1. reconstruct `(log rho, u, log p, log T_e)` with the asymmetric MC3/Koren
 *     limiter (beta = 2) on the non-uniform-mesh-aware slope metrics — the
 *     release primitive set plus ONE extra slot for the electron temperature.
 *     The TOTAL pressure remains the authoritative mechanical variable, so the
 *     new slot changes only the electron/heavy energy PARTITION at a face;
 *  2. take a MUSCL-Hancock predictor half step applying the same sources as the
 *     corrector (flux-tube pressure term, gravity, electron compression work);
 *  3. build every one-sided face state algebraically — no temperature inversion
 *     anywhere, because (rho, p, T_e) closes the state in closed form;
 *  4. evaluate the release SWMF-style exact-Riemann Godunov flux at the frozen
 *     index gamma = 5/3 on (rho, rho u, E), and upwind the specific electron
 *     energy E_e/rho on the contact with the exact mass flux;
 *  5. difference the fluxes with the flux-tube area factor and add the sources.
 *
 * ## Why the release flux is valid here unchanged
 *
 * Both species are monatomic, so at frozen composition
 * `e_int = e_h + e_e = p/(5/3 - 1) + e_ion` regardless of how the thermal energy
 * is partitioned. The exact ideal-gas Riemann solve at gamma = 5/3 with `e_ion`
 * in SWMF's passive offset E0 therefore describes the acoustic waves of the 2T
 * system exactly as it describes those of the release system, and the frozen
 * sound speed sqrt(5/3 p/rho) is unchanged. The extra physics of the 2T system
 * lives on the SECOND lambda = u characteristic (the electron/heavy partition),
 * which is precisely what step 4's contact upwinding transports. The side that
 * supplies E0 and the side that supplies E_e/rho are the SAME side, so the total
 * energy flux and the electron energy flux never disagree about which state's
 * composition crossed the face.
 *
 * The electron equation is advanced in the ADVECTIVE form
 *
 *     d_t E_e + (1/A) d_s(A E_e u) = -p_e (1/A) d_s(A u) + ...
 *
 * i.e. a mass-consistent flux plus a compression-work source, which is the
 * standard non-conservative electron-pressure treatment of AWSoM/BATS-R-US. The
 * TOTAL energy row stays in exact conservation form and the heavy internal energy
 * is obtained by difference, so total energy is conserved exactly and any
 * splitting error is absorbed by the heavy pool rather than by the total.
 */

#include "two_temp/two_temp.hpp"
#include "single_fluid/exact_rs.hpp"
#include "parallel.hpp"
#include "profiling.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace chromosphere {

void TwoTempFaceArrays::resize(arma::uword ns) {
    conserved.set_size(ns*num_of_two_temp_eq);
    flux.set_size(ns*num_of_two_temp_eq);
    velocity.set_size(ns);
    pressure.set_size(ns);
    pressure_e.set_size(ns);
    sound_speed.set_size(ns);
}

void TwoTempWorkspace::resize(arma::uword ns) {
    for (arma::Col<double>& v : work) v.set_size(ns*num_of_two_temp_eq);
    for (TwoTempFaceArrays& f : faces) f.resize(ns);
    flux_iph.set_size(ns*num_of_two_temp_eq);
    flux_imh.set_size(ns*num_of_two_temp_eq);
    source.set_size(ns*num_of_two_temp_eq);
    u_face_iph.set_size(ns);
    u_face_imh.set_size(ns);
    predicted_state.set_size(ns*num_of_two_temp_eq);
    rhs.set_size(ns*num_of_two_temp_eq);
    const std::size_t n = ns;
    std::vector<double>* cell_arrays[] = {
        &c_rho, &c_Te, &c_Ti, &c_ne, &c_nHI, &c_nH,
        &c_Ee_old, &c_eh_old, &c_Ee_at_T, &c_cap_e, &c_cap_i,
        &c_ke, &c_ki, &c_gE_L, &c_gE_R, &c_gI_L, &c_gI_R, &c_gex,
        &c_Ee_target, &c_eh_target};
    for (std::vector<double>* v : cell_arrays) v->assign(n, 0.0);
    blk_a.assign(2*n, 0.0);
    blk_c.assign(2*n, 0.0);
    blk_b.assign(4*n, 0.0);
    blk_rhs.assign(2*n, 0.0);
    blk_delta.assign(2*n, 0.0);
}

namespace {

using Work = arma::Col<double>;

/// Frozen-composition index of the local Riemann problem. Identical to the
/// release: both species are monatomic, so the translational gas is gamma = 5/3
/// whatever the electron/heavy energy partition is.
constexpr double kGodunovGamma = 5.0/3.0;

inline arma::SizeMat state_size(const Grid& grid) {
    return arma::size(grid.ns, num_of_two_temp_eq);
}

/// Neumann-mirrored neighbour shift of a cell-centred STATIC mesh metric. Only
/// mesh metrics are shifted this way, so no ghost data is involved.
Vec shift_up(const Vec& v) {
    Vec out(arma::size(v));
    const arma::uword n = v.n_elem;
    for (arma::uword i = 0; i + 1 < n; ++i) out(i) = v(i + 1);
    out(n - 1) = v(n - 1);
    return out;
}

Vec shift_down(const Vec& v) {
    Vec out(arma::size(v));
    const arma::uword n = v.n_elem;
    for (arma::uword i = 1; i < n; ++i) out(i) = v(i - 1);
    out(0) = v(0);
    return out;
}

void stencil_view(const Grid& grid, const TwoTempField& decoded, int offset,
                  Work& result) {
    result.set_size(grid.ns*num_of_two_temp_eq);
    const auto ext_size = arma::size(grid.ns + 4, num_of_two_temp_eq);
    const auto packed = state_size(grid);
    for (arma::uword slot = 0; slot < num_of_two_temp_eq; ++slot)
        for (arma::uword i = 0; i < grid.ns; ++i) {
            const arma::uword ext_i =
                static_cast<arma::uword>(static_cast<int>(i) + 2 + offset);
            result(arma::sub2ind(packed, i, slot)) =
                decoded.extended_primitive(arma::sub2ind(ext_size, ext_i, slot));
        }
}

void broadcast_slots(const Grid& grid, const Vec& values, Work& packed) {
    packed.set_size(grid.ns*num_of_two_temp_eq);
    for (arma::uword k = 0; k < num_of_two_temp_eq; ++k)
        for (arma::uword i = 0; i < grid.ns; ++i)
            packed(arma::sub2ind(state_size(grid), i, k)) = values(i);
}

void limit_minmod(const Work& ratio, Work& result) {
    result.set_size(ratio.n_elem);
    for (arma::uword i = 0; i < ratio.n_elem; ++i) {
        const double r = ratio(i);
        result(i) = std::isfinite(r) ? std::max(0.0, std::min(1.0, r)) : 0.0;
    }
}

void limit_mc3(const Work& ratio, double beta, bool plus, Work& result) {
    result.set_size(ratio.n_elem);
    for (arma::uword i = 0; i < ratio.n_elem; ++i) {
        const double r = ratio(i);
        if (!std::isfinite(r)) { result(i) = 0.0; continue; }
        const double third_order = plus ? (2.0*r + 1.0)/3.0 : (r + 2.0)/3.0;
        result(i) = std::max(0.0, std::min(std::min(beta*r, beta), third_order));
    }
}

/// Build one one-sided face state per cell from the reconstructed variables and
/// fill both the conserved rows and the physical flux at the face potential.
void build_face_cell(const Grid& grid, const Work& w, const Vec& phi_face,
                     TwoTempFaceArrays& out, arma::uword i,
                     std::uint64_t& clamps) {
    const auto sz = state_size(grid);
    const double rho = std::exp(w(arma::sub2ind(sz, i, TT_LOG_RHO)));
    const double velocity = w(arma::sub2ind(sz, i, TT_U));
    const double pressure = std::exp(w(arma::sub2ind(sz, i, TT_LOG_P)));
    const double T_e = std::exp(w(arma::sub2ind(sz, i, TT_LOG_TE)));
    bool clamped = false;
    const TwoTempThermo face = two_temp_face_state(
        grid.eos_gamma_table, rho, velocity, pressure, T_e, &clamped);
    if (clamped) ++clamps;

    const double phi = static_cast<double>(phi_face(i));
    const double momentum = face.rho*velocity;
    const double energy = two_temp_total_energy(
        face.rho, momentum, face.e_heavy, face.e_elec, phi);
    const double conserved[num_of_two_temp_eq] = {
        face.rho, momentum, energy, face.e_elec};
    const double flux[num_of_two_temp_eq] = {
        momentum,
        momentum*velocity + face.p,
        (energy + face.p)*velocity,
        face.e_elec*velocity};
    for (arma::uword k = 0; k < num_of_two_temp_eq; ++k) {
        out.conserved(arma::sub2ind(sz, i, k)) = conserved[k];
        out.flux(arma::sub2ind(sz, i, k)) = flux[k];
    }
    out.velocity(i) = velocity;
    out.pressure(i) = face.p;
    out.pressure_e(i) = face.p_e;
    out.sound_speed(i) = face.sound_speed;
}

struct FaceRequest {
    const Work* w;
    const Vec* phi;
    TwoTempFaceArrays* out;
};

void build_faces(const Grid& grid, const FaceRequest* requests, int count,
                 TwoTempStats& stats) {
    for (int r = 0; r < count; ++r) requests[r].out->resize(grid.ns);
    std::array<std::uint64_t, kMaximumParallelThreads> clamp_slots{};
    ParallelFailure failure;
    parallel_for_cells(grid.ns, [&](std::size_t raw_i) {
        const arma::uword i = static_cast<arma::uword>(raw_i);
        const std::size_t tid = static_cast<std::size_t>(parallel_thread_index());
        try {
            for (int r = 0; r < count; ++r)
                build_face_cell(grid, *requests[r].w, *requests[r].phi,
                                *requests[r].out, i, clamp_slots[tid]);
        } catch (...) {
            failure.capture(i, std::current_exception());
        }
    });
    failure.rethrow_lowest();
    for (std::uint64_t count_slot : clamp_slots) stats.face_pe_clamps += count_slot;
}

/// Face-local Rusanov (local Lax-Friedrichs) numerical flux on all four rows,
/// with the arithmetic-mean face velocity for the compression-work source. This
/// is the per-face fallback whenever the exact Riemann solve is not usable.
inline void rusanov_face(const Grid& grid, const TwoTempFaceArrays& left,
                         const TwoTempFaceArrays& right, arma::uword i,
                         double a, Work& output, Work& u_face) {
    const auto sz = state_size(grid);
    for (arma::uword row = 0; row < num_of_two_temp_eq; ++row) {
        const arma::uword k = arma::sub2ind(sz, i, row);
        output(k) = 0.5*(left.flux(k) + right.flux(k)
                         - a*(right.conserved(k) - left.conserved(k)));
    }
    u_face(i) = 0.5*(left.velocity(i) + right.velocity(i));
}

/// The release SWMF-style Godunov flux on (rho, rho u, E), plus contact
/// upwinding of the specific electron energy on the exact mass flux.
void build_godunov_flux(const Grid& grid, const TwoTempFaceArrays& left,
                        const TwoTempFaceArrays& right, const Work& a_face,
                        Work& output, Work& u_face, TwoTempStats& stats) {
    const auto sz = state_size(grid);
    constexpr double gm1_inv = 1.0/(kGodunovGamma - 1.0);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const arma::uword k_rho = arma::sub2ind(sz, i, tt::RHO);
        const arma::uword k_mom = arma::sub2ind(sz, i, tt::MOM);
        const arma::uword k_e   = arma::sub2ind(sz, i, tt::ENERGY);
        const arma::uword k_ee  = arma::sub2ind(sz, i, tt::E_ELEC);
        ++stats.godunov_faces;

        const swmf_rs::ExactRsState state_l{left.conserved(k_rho),
                                            left.velocity(i), left.pressure(i)};
        const swmf_rs::ExactRsState state_r{right.conserved(k_rho),
                                            right.velocity(i), right.pressure(i)};
        const double energy_l = left.conserved(k_e);
        const double energy_r = right.conserved(k_e);
        auto fall_back = [&]() {
            ++stats.godunov_fallbacks;
            rusanov_face(grid, left, right, i, a_face(i), output, u_face);
        };
        if (!std::isfinite(energy_l) || !std::isfinite(energy_r)) {
            fall_back();
            continue;
        }

        const swmf_rs::ExactRsSolution solution = swmf_rs::exact_rs_pu_star(
            state_l, state_r, kGodunovGamma, kGodunovGamma);
        if (solution.status != swmf_rs::ExactRsStatus::Ok) { fall_back(); continue; }

        // ONE upwind decision serves both the passive energy offset and the
        // electron energy, so the total-energy flux and the electron-energy flux
        // always carry the same side's composition. SWMF's tie-break at
        // u_star == 0 exactly is the right state.
        const bool from_left = solution.u_star > 0.0;
        const double side_rho = from_left ? state_l.rho : state_r.rho;
        const double side_u   = from_left ? state_l.u   : state_r.u;
        const double side_p   = from_left ? state_l.p   : state_r.p;
        const double side_energy = from_left ? energy_l : energy_r;
        const double side_e_elec = from_left ? left.conserved(k_ee)
                                             : right.conserved(k_ee);
        const double e0 =
            (side_energy - 0.5*side_rho*side_u*side_u - side_p*gm1_inv)/side_rho;
        const double eps_e = side_e_elec/side_rho;   // specific electron energy [J/kg]

        const swmf_rs::ExactRsState face = swmf_rs::exact_rs_sample(
            0.0, solution, state_l, state_r, kGodunovGamma, kGodunovGamma);
        if (!(face.rho > 0.0) || !(face.p > 0.0) || !std::isfinite(face.u)
            || !std::isfinite(e0) || !std::isfinite(eps_e)) {
            fall_back();
            continue;
        }

        const double energy_star =
            face.p*gm1_inv + 0.5*face.rho*face.u*face.u + face.rho*e0;
        const double mass_flux = face.rho*face.u;
        const double flux[num_of_two_temp_eq] = {
            mass_flux,
            mass_flux*face.u + face.p,
            (energy_star + face.p)*face.u,
            mass_flux*eps_e};
        bool admissible = true;
        for (double value : flux) admissible = admissible && std::isfinite(value);
        if (!admissible) { fall_back(); continue; }
        output(k_rho) = flux[0];
        output(k_mom) = flux[1];
        output(k_e)   = flux[2];
        output(k_ee)  = flux[3];
        u_face(i) = face.u;
    }
}

/// Momentum source of the field-aligned PDE (flux-tube pressure term plus
/// gravity), plus the electron compression work -p_e (1/A) d_s(A u) built from
/// the face velocities of the stage that is calling. The energy row has no source
/// because E carries the gravitational potential.
void build_source(const Grid& grid, const Vec& state,
                  const TwoTempField& decoded, const Work& u_iph,
                  const Work& u_imh, Work& source) {
    source.zeros(grid.ns*num_of_two_temp_eq);
    const auto sz = state_size(grid);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const double g_par = -static_cast<double>(
            (grid.phi_g_iph(i) - grid.phi_g_imh(i))/grid.ds_i(i));
        source(arma::sub2ind(sz, i, tt::MOM)) =
            decoded.cells[i].p*static_cast<double>(grid.B_i(i))
                *static_cast<double>(grid.dinvB_ds_i(i))
            + static_cast<double>(state(arma::sub2ind(sz, i, tt::RHO)))*g_par;
        const double div_u = static_cast<double>(grid.B_i(i))
            * (u_iph(i)/static_cast<double>(grid.B_iph(i))
             - u_imh(i)/static_cast<double>(grid.B_imh(i)))
            / static_cast<double>(grid.ds_i(i));
        source(arma::sub2ind(sz, i, tt::E_ELEC)) = -decoded.cells[i].p_e*div_u;
    }
}

} // namespace

// ============================================================================
// Explicit RHS
// ============================================================================

Vec two_temp_rhs_explicit(const Grid& grid, const Vec& state,
                          const TwoTempField& decoded, double dt_predictor,
                          TwoTempWorkspace& scratch) {
    ProfileScope timer(ProfileRegion::Rhs);
    if (decoded.ns != grid.ns || decoded.cells.size() != grid.ns)
        throw std::logic_error("two-temperature field does not match the grid");

    auto& m = scratch.work;
    Work& w=m[0]; Work& w_ip1=m[1]; Work& w_im1=m[2];
    Work& w_ip2=m[3]; Work& w_im2=m[4];
    Work& r=m[5]; Work& r_ip1=m[6]; Work& r_im1=m[7];
    Work& wr_iph=m[8]; Work& wl_iph=m[9];
    Work& wr_imh=m[10]; Work& wl_imh=m[11];
    Work& wt=m[12]; Work& wt_ip1=m[13]; Work& wt_im1=m[14];
    Work& lp_r=m[15]; Work& lm_r=m[16];
    Work& lm_rip1=m[17]; Work& lp_rim1=m[18];
    Work& dw_iph=m[19]; Work& dw_imh=m[20];
    Work& a_imh=m[21]; Work& a_iph=m[22];

    stencil_view(grid, decoded,  0, w);
    stencil_view(grid, decoded,  1, w_ip1);
    stencil_view(grid, decoded, -1, w_im1);
    stencil_view(grid, decoded,  2, w_ip2);
    stencil_view(grid, decoded, -2, w_im2);

    dw_iph = w_ip1 - w;
    dw_imh = w - w_im1;
    r = dw_imh/dw_iph;
    r_ip1 = dw_iph/(w_ip2 - w_ip1);
    r_im1 = (w_im1 - w_im2)/dw_imh;

    if (!grid.uniform_mesh) {
        if (!scratch.weights_valid
            || scratch.weights_generation != grid.metrics_generation()
            || scratch.W1.n_elem != grid.ns*num_of_two_temp_eq) {
            const Vec w1_i = 0.5*grid.ds_i/grid.ds_iph_i;
            broadcast_slots(grid, w1_i, scratch.W1);
            broadcast_slots(grid, shift_up(w1_i), scratch.W3);
            broadcast_slots(grid, shift_down(w1_i), scratch.W4);
            broadcast_slots(grid, grid.ds_iph_i/grid.ds_imh_i, scratch.metric_r);
            broadcast_slots(grid, shift_up(grid.ds_iph_i)/grid.ds_iph_i,
                            scratch.metric_r_ip1);
            broadcast_slots(grid, grid.ds_imh_i/shift_down(grid.ds_imh_i),
                            scratch.metric_r_im1);
            scratch.weights_generation = grid.metrics_generation();
            scratch.weights_valid = true;
        }
        r %= scratch.metric_r;
        r_ip1 %= scratch.metric_r_ip1;
        r_im1 %= scratch.metric_r_im1;
    }
    r.elem(arma::find_nonfinite(r)).zeros();
    r_ip1.elem(arma::find_nonfinite(r_ip1)).zeros();
    r_im1.elem(arma::find_nonfinite(r_im1)).zeros();

    const double beta = grid.limiter_beta;
    auto lim_plus = [&](const Work& q, Work& out) {
        if (grid.mc3_limiter) limit_mc3(q, beta, true, out);
        else limit_minmod(q, out);
    };
    auto lim_minus = [&](const Work& q, Work& out) {
        if (grid.mc3_limiter) limit_mc3(q, beta, false, out);
        else limit_minmod(q, out);
    };
    lim_plus(r, lp_r);         lim_minus(r, lm_r);
    lim_minus(r_ip1, lm_rip1); lim_plus(r_im1, lp_rim1);

    const Work& W1 = scratch.W1;
    const Work& W3 = scratch.W3;
    const Work& W4 = scratch.W4;

    // ---- Predictor pass ---------------------------------------------------
    // The Hancock half step differences only cell i's OWN two extrapolated
    // fluxes, so only those two reconstructions are needed.
    if (grid.uniform_mesh) {
        wl_iph = w + 0.5*lp_r%(w_ip1 - w);
        wr_imh = w - 0.5*lm_r%(w_ip1 - w);
    } else {
        wl_iph = w + W1%lp_r%(w_ip1 - w);
        wr_imh = w - W1%lm_r%(w_ip1 - w);
    }

    TwoTempFaceArrays& fr_iph = scratch.faces[0];
    TwoTempFaceArrays& fl_iph = scratch.faces[1];
    TwoTempFaceArrays& fr_imh = scratch.faces[2];
    TwoTempFaceArrays& fl_imh = scratch.faces[3];
    {
        const FaceRequest requests[2] = {
            {&wl_iph, &grid.phi_g_iph, &fl_iph},
            {&wr_imh, &grid.phi_g_imh, &fr_imh}};
        build_faces(grid, requests, 2, scratch.stats);
    }

    // Both stages apply the same sources; the predictor uses cell i's own
    // one-sided face velocities, exactly as it differences cell i's own fluxes.
    Work& u_iph = scratch.u_face_iph;
    Work& u_imh = scratch.u_face_imh;
    u_iph.set_size(grid.ns);
    u_imh.set_size(grid.ns);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        u_iph(i) = fl_iph.velocity(i);
        u_imh(i) = fr_imh.velocity(i);
    }
    Work& source = scratch.source;
    build_source(grid, state, decoded, u_iph, u_imh, source);

    const auto sz = state_size(grid);
    Vec& predicted = scratch.predicted_state;
    predicted.set_size(grid.ns*num_of_two_temp_eq);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const double scale = dt_predictor/static_cast<double>(grid.ds_i(i));
        for (arma::uword row = 0; row < num_of_two_temp_eq; ++row) {
            const arma::uword k = arma::sub2ind(sz, i, row);
            predicted(k) = static_cast<Real>(
                static_cast<double>(state(k))
                - scale*(fl_iph.flux(k) - fr_imh.flux(k))
                + dt_predictor*source(k));
        }
    }
    two_temp_decode_into(grid, predicted, scratch.predicted, scratch.stats,
                         &decoded);
    stencil_view(grid, scratch.predicted,  0, wt);
    stencil_view(grid, scratch.predicted,  1, wt_ip1);
    stencil_view(grid, scratch.predicted, -1, wt_im1);

    // ---- Corrector pass ---------------------------------------------------
    if (grid.uniform_mesh) {
        wr_iph = 0.5*(w_ip1 + wt_ip1) - 0.5*lm_rip1%(w_ip2 - w_ip1);
        wl_iph = 0.5*(w + wt)         + 0.5*lp_r%(w_ip1 - w);
        wr_imh = 0.5*(w + wt)         - 0.5*lm_r%(w_ip1 - w);
        wl_imh = 0.5*(w_im1 + wt_im1) + 0.5*lp_rim1%(w - w_im1);
    } else {
        wr_iph = 0.5*(w_ip1 + wt_ip1) - W3%lm_rip1%(w_ip2 - w_ip1);
        wl_iph = 0.5*(w + wt)         + W1%lp_r%(w_ip1 - w);
        wr_imh = 0.5*(w + wt)         - W1%lm_r%(w_ip1 - w);
        wl_imh = 0.5*(w_im1 + wt_im1) + W4%lp_rim1%(w - w_im1);
    }
    {
        const FaceRequest requests[4] = {
            {&wr_iph, &grid.phi_g_iph, &fr_iph},
            {&wl_iph, &grid.phi_g_iph, &fl_iph},
            {&wr_imh, &grid.phi_g_imh, &fr_imh},
            {&wl_imh, &grid.phi_g_imh, &fl_imh}};
        build_faces(grid, requests, 4, scratch.stats);
    }

    a_imh.set_size(grid.ns);
    a_iph.set_size(grid.ns);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        a_imh(i) = std::max(
            std::abs(fl_imh.velocity(i)) + fl_imh.sound_speed(i),
            std::abs(fr_imh.velocity(i)) + fr_imh.sound_speed(i));
        a_iph(i) = std::max(
            std::abs(fl_iph.velocity(i)) + fl_iph.sound_speed(i),
            std::abs(fr_iph.velocity(i)) + fr_iph.sound_speed(i));
    }
    Work& flux_imh = scratch.flux_imh;
    Work& flux_iph = scratch.flux_iph;
    flux_imh.set_size(grid.ns*num_of_two_temp_eq);
    flux_iph.set_size(grid.ns*num_of_two_temp_eq);
    build_godunov_flux(grid, fl_imh, fr_imh, a_imh, flux_imh, u_imh,
                       scratch.stats);
    build_godunov_flux(grid, fl_iph, fr_iph, a_iph, flux_iph, u_iph,
                       scratch.stats);

    // Rebuild the sources on the corrector's Godunov face velocities.
    build_source(grid, state, decoded, u_iph, u_imh, source);

    Vec& rhs = scratch.rhs;
    rhs.set_size(grid.ns*num_of_two_temp_eq);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const double B_i = grid.B_i(i);
        const double inv_B_iph = 1.0/static_cast<double>(grid.B_iph(i));
        const double inv_B_imh = 1.0/static_cast<double>(grid.B_imh(i));
        const double inv_ds = 1.0/static_cast<double>(grid.ds_i(i));
        for (arma::uword row = 0; row < num_of_two_temp_eq; ++row) {
            const arma::uword k = arma::sub2ind(sz, i, row);
            rhs(k) = static_cast<Real>(
                -B_i*(flux_iph(k)*inv_B_iph - flux_imh(k)*inv_B_imh)*inv_ds
                + source(k));
        }
    }
    return Vec(rhs);
}

// ============================================================================
// Timestep
// ============================================================================

Vec two_temp_timestep(const Grid& grid, const Vec& state,
                      const TwoTempField& decoded) {
    ProfileScope timer(ProfileRegion::Cfl);
    if (decoded.ns != grid.ns)
        throw std::logic_error("two-temperature field does not match the grid");
    const auto sz = state_size(grid);
    Vec dt_i(grid.ns);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const TwoTempThermo& th = decoded.cells[i];
        const double momentum =
            static_cast<double>(state(arma::sub2ind(sz, i, tt::MOM)));
        // The frozen acoustic speed depends only on the TOTAL pressure, so the
        // CFL rule is byte-for-byte the release Godunov rule. The stiff
        // collisional exchange (tau_eq of order 0.1-4 ms here, far below dt)
        // imposes no explicit limit because it is solved implicitly.
        const double speed = std::abs(momentum/th.rho) + th.sound_speed;
        dt_i(i) = static_cast<Real>(grid.CFL*grid.ds_i(i)/speed);
    }
    profile_note_timestep_limiter(TimestepLimiter::Acoustic);
    return arma::min(dt_i)*arma::ones<Vec>(grid.ns);
}

// ============================================================================
// One complete step
// ============================================================================

Vec two_temp_advance(Grid& grid, const Vec& state, const Vec& dt_i,
                     const TwoTempField& decoded, TwoTempWorkspace& work) {
    if (grid.eos_gamma_table.empty())
        throw std::logic_error("the two-temperature solver requires a Gamma1 table");
    const double dt = dt_i(0);   // uniform by two_temp_timestep construction

    Vec next = state;
    {
        const Vec rhs = two_temp_rhs_explicit(grid, state, decoded, dt, work);
        const auto sz = state_size(grid);
        for (arma::uword i = 0; i < grid.ns; ++i)
            for (arma::uword row = 0; row < num_of_two_temp_eq; ++row) {
                const arma::uword k = arma::sub2ind(sz, i, row);
                next(k) = state(k) + dt_i(i)*rhs(k);
            }
    }
    if (grid.enable_conduction)
        next = two_temp_apply_conduction(grid, next, dt, work);
    return next;
}

} // namespace chromosphere
