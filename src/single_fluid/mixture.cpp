// Release solver: decode, MUSCL reconstruction, mixture Roe flux, source.
// See single_fluid/mixture.hpp for the governing equations and the software
// contract. Nothing here reaches into the two-fluid research solver.

#include "single_fluid/mixture.hpp"
#include "profiling.hpp"
#include "parallel.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace chromosphere {

void MixtureFaceArrays::resize(arma::uword ns) {
    conserved.set_size(ns*num_of_mixture_eq);
    flux.set_size(ns*num_of_mixture_eq);
    velocity.set_size(ns);
    pressure.set_size(ns);
    sound_speed.set_size(ns);
    dp_deint.set_size(ns);
}

namespace {

// Reconstruction slots. Slot 2 carries the THERMAL variable: log(p) with the
// release reconstruction, log(T) with the reference one. The packing, limiter
// arithmetic and extrapolation are identical either way; only the decode (which
// quantity is written) and the face builder (how the state is closed) branch.
enum MixtureSlot : arma::uword { MIX_LOG_RHO = 0, MIX_U = 1, MIX_THERMAL = 2 };

using Work = arma::Col<double>;

inline arma::SizeMat state_size(const Grid& grid) {
    return arma::size(grid.ns, num_of_mixture_eq);
}

// Neumann-mirrored neighbour shift of a cell-centred STATIC mesh metric:
// shift_up(v)(i) = v(i+1) in the interior and v(ns-1) at the top face;
// shift_down mirrors at the base. Only mesh metrics are ever shifted this way,
// so no ghost data is involved and the release solver needs none of the
// two-fluid packed-state helpers.
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

} // namespace

// ============================================================================
// Decode
// ============================================================================

void MixtureField::require_matches(const Grid& grid, const Vec& state) const {
    if (source_grid != &grid || source_state != &state
        || source_memory != state.memptr()
        || source_elements != state.n_elem)
        throw std::logic_error(
            "mixture field cache does not match the immutable conserved state");
    const Vec* boundaries[] = {&grid.mix_inner_boundary0, &grid.mix_inner_boundary1,
                               &grid.mix_outer_boundary0, &grid.mix_outer_boundary1};
    for (arma::uword block = 0; block < 4; ++block)
        for (arma::uword row = 0; row < num_of_mixture_eq; ++row)
            if (boundary_signature[block*num_of_mixture_eq+row] != (*boundaries[block])(row))
                throw std::logic_error(
                    "mixture field cache was invalidated by a boundary update");
}

MixtureField mixture_decode(const Grid& grid, const Vec& state,
                            std::uint64_t generation,
                            const MixtureField* previous) {
    MixtureField out;
    mixture_decode_into(grid, state, out, generation, previous);
    return out;
}

void mixture_decode_into(const Grid& grid, const Vec& state, MixtureField& out,
                         std::uint64_t generation, const MixtureField* previous) {
    ProfileScope timer(ProfileRegion::Decode);
    if (grid.eos_gamma_table.empty())
        throw std::logic_error("the mixture solver requires a Gamma1 table");
    if (state.n_elem != grid.n_mixture_state)
        throw std::invalid_argument("mixture state size mismatch");

    out.source_grid = &grid;
    out.source_state = &state;
    out.source_memory = state.memptr();
    out.source_elements = state.n_elem;
    out.state_generation = generation;
    const Vec* boundaries[] = {&grid.mix_inner_boundary0, &grid.mix_inner_boundary1,
                               &grid.mix_outer_boundary0, &grid.mix_outer_boundary1};
    for (arma::uword block = 0; block < 4; ++block)
        for (arma::uword row = 0; row < num_of_mixture_eq; ++row)
            out.boundary_signature[block*num_of_mixture_eq+row] = (*boundaries[block])(row);
    out.cells.resize(grid.ns);
    const arma::uword ext_n = grid.ns+4;
    out.extended_primitive.set_size(ext_n*num_of_mixture_eq);
    out.extended_temperature.set_size(ext_n);
    const auto sz = state_size(grid);
    const auto ext_size = arma::size(ext_n, num_of_mixture_eq);
    const bool log_p = grid.pressure_reconstruct;
    auto put = [&](arma::uword ext_i, double rho, double momentum,
                   double pressure, double temperature) {
        out.extended_primitive(arma::sub2ind(ext_size, ext_i, MIX_LOG_RHO)) =
            std::log(rho);
        out.extended_primitive(arma::sub2ind(ext_size, ext_i, MIX_U)) =
            momentum/rho;
        out.extended_primitive(arma::sub2ind(ext_size, ext_i, MIX_THERMAL)) =
            log_p ? std::log(pressure) : std::log(temperature);
        out.extended_temperature(ext_i) = temperature;
    };
    const bool has_guesses = previous && previous->source_grid == &grid
        && previous->cells.size() == grid.ns;
    parallel_for_cells(grid.ns, [&](std::size_t raw_i) {
        const arma::uword i = static_cast<arma::uword>(raw_i);
        auto at = [&](arma::uword row) {
            return static_cast<double>(state(arma::sub2ind(sz, i, row)));
        };
        const double guess = has_guesses ? previous->cells[i].T
            : std::numeric_limits<double>::quiet_NaN();
        out.cells[i] = decode_equilibrium_mixture(
            grid.eos_gamma_table, at(mix::RHO), at(mix::MOM), at(mix::ENERGY),
            mixture_cell_phi(grid, i), guess, grid.eos_gamma_debug_clamp);
        // Refresh the shared hint, but do NOT read it here: this site has its own
        // guess (`previous`), validated against source_grid and cell count. Each
        // worker writes only its own pre-sized hint slot.
        grid.store_eos_temperature_hint(i, out.cells[i].T);
        const MixtureThermo& th = out.cells[i];
        put(i+2, th.rho, at(mix::MOM), th.p, th.T);
    });

    auto decode_ghost = [&](const Vec& ghost, double phi, arma::uword ext_i) {
        const MixtureThermo th = decode_equilibrium_mixture(
            grid.eos_gamma_table, ghost(mix::RHO), ghost(mix::MOM),
            ghost(mix::ENERGY), phi, std::numeric_limits<double>::quiet_NaN(),
            grid.eos_gamma_debug_clamp);
        put(ext_i, th.rho, static_cast<double>(ghost(mix::MOM)), th.p, th.T);
    };
    const double phi_inner = grid.phi_g_imh(0);
    const double phi_outer = grid.phi_g_iph(grid.ns-1);
    decode_ghost(grid.mix_inner_boundary1, phi_inner, 0);
    decode_ghost(grid.mix_inner_boundary0, phi_inner, 1);
    decode_ghost(grid.mix_outer_boundary0, phi_outer, grid.ns+2);
    decode_ghost(grid.mix_outer_boundary1, phi_outer, grid.ns+3);
}

void mixture_pack_ghost(const Grid& grid, Vec& ghost, double rho,
                        double velocity, double temperature, double phi) {
    if (ghost.n_elem != num_of_mixture_eq) ghost.zeros(num_of_mixture_eq);
    grid.eos_gamma_table.require_n_h_in_bounds(
        rho/eos_constants::m_h, grid.eos_gamma_debug_clamp);
    const double e_int = equilibrium_internal_energy(rho, temperature);
    ghost(mix::RHO)    = static_cast<float>(rho);
    ghost(mix::MOM)    = static_cast<float>(rho*velocity);
    ghost(mix::ENERGY) = static_cast<float>(
        mixture_total_energy(rho, rho*velocity, e_int, phi));
}

void mixture_pack_cell(const Grid& grid, Vec& state, arma::uword cell,
                       double rho, double velocity, double temperature,
                       double phi) {
    Vec packed(num_of_mixture_eq, arma::fill::zeros);
    mixture_pack_ghost(grid, packed, rho, velocity, temperature, phi);
    const auto sz = state_size(grid);
    for (arma::uword row = 0; row < num_of_mixture_eq; ++row)
        state(arma::sub2ind(sz, cell, row)) = packed(row);
}

// ============================================================================
// MUSCL reconstruction + numerical flux
// ============================================================================

namespace {

void stencil_view(const Grid& grid, const MixtureField& decoded, int offset,
                  Work& result) {
    result.set_size(grid.ns*num_of_mixture_eq);
    const arma::uword ext_n = grid.ns+4;
    const auto ext_size = arma::size(ext_n, num_of_mixture_eq);
    const auto packed = state_size(grid);
    for (arma::uword slot = 0; slot < num_of_mixture_eq; ++slot)
        for (arma::uword i = 0; i < grid.ns; ++i) {
            const arma::uword ext_i =
                static_cast<arma::uword>(static_cast<int>(i)+2+offset);
            result(arma::sub2ind(packed, i, slot)) =
                decoded.extended_primitive(arma::sub2ind(ext_size, ext_i, slot));
        }
}

// Parent-cell temperature hints, shifted by exactly the same stencil offset the
// reconstruction of the corresponding one-sided face state used. This is solver
// scratch, not a reconstructed variable: it never enters the limiter, the
// extrapolation, or any flux.
void temperature_stencil_view(const Grid& grid, const MixtureField& decoded,
                              int offset, Work& result) {
    result.set_size(grid.ns);
    for (arma::uword i = 0; i < grid.ns; ++i)
        result(i) = decoded.extended_temperature(
            static_cast<arma::uword>(static_cast<int>(i)+2+offset));
}

void broadcast_slots(const Grid& grid, const Vec& values, Work& packed) {
    packed.set_size(grid.ns*num_of_mixture_eq);
    for (arma::uword k = 0; k < num_of_mixture_eq; ++k)
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

// Build one one-sided face state per cell from the reconstructed variables.
void build_face_cell(const Grid& grid, const Work& w, const Vec& phi_face,
                     MixtureFaceArrays& out, arma::uword i,
                     const Work* temperature_hint) {
    const auto sz = state_size(grid);
    const double log_rho = w(arma::sub2ind(sz, i, MIX_LOG_RHO));
    const double velocity = w(arma::sub2ind(sz, i, MIX_U));
    const double thermal = w(arma::sub2ind(sz, i, MIX_THERMAL));
    const MixtureFaceState face = grid.pressure_reconstruct
        ? equilibrium_mixture_face_state_from_log_pressure(
              grid.eos_gamma_table, log_rho, velocity, thermal,
              static_cast<double>(phi_face(i)), grid.eos_gamma_debug_clamp,
              temperature_hint ? (*temperature_hint)(i)
                               : std::numeric_limits<double>::quiet_NaN())
        : equilibrium_mixture_face_state_from_logs(
              grid.eos_gamma_table, log_rho, velocity, thermal,
              static_cast<double>(phi_face(i)), grid.eos_gamma_debug_clamp);
    const std::array<double, num_of_mixture_eq> flux = equilibrium_mixture_flux(face);
    const double conserved[num_of_mixture_eq] = {
        face.rho, face.rho*face.velocity, face.energy};
    for (arma::uword k = 0; k < num_of_mixture_eq; ++k) {
        out.conserved(arma::sub2ind(sz, i, k)) = conserved[k];
        out.flux(arma::sub2ind(sz, i, k)) = flux[k];
    }
    out.velocity(i) = face.velocity;
    out.pressure(i) = face.pressure;
    out.sound_speed(i) = face.sound_speed;
    out.dp_deint(i) = face.dp_deint_rho;
}

struct FaceRequest {
    const Work* w;
    const Vec* phi;
    MixtureFaceArrays* out;
    const Work* hint;
};

void build_faces(const Grid& grid, const FaceRequest* requests, int count) {
    for (int r = 0; r < count; ++r) requests[r].out->resize(grid.ns);
    parallel_for_cells(grid.ns, [&](std::size_t raw_i) {
        const arma::uword i = static_cast<arma::uword>(raw_i);
        for (int r = 0; r < count; ++r)
            build_face_cell(grid, *requests[r].w, *requests[r].phi,
                            *requests[r].out, i, requests[r].hint);
    });
}

// Face-local Rusanov (local Lax-Friedrichs) numerical flux.
inline void rusanov_face(const Grid& grid, const MixtureFaceArrays& left,
                         const MixtureFaceArrays& right, arma::uword i,
                         double a, Work& output) {
    const auto sz = state_size(grid);
    for (arma::uword row = 0; row < num_of_mixture_eq; ++row) {
        const arma::uword k = arma::sub2ind(sz, i, row);
        output(k) = 0.5*(left.flux(k) + right.flux(k)
                         - a*(right.conserved(k) - left.conserved(k)));
    }
}

// Mixture Roe characteristic flux on U = (rho, rho u, E).
//
// With the Roe-averaged velocity u~, total enthalpy H~ and the general-EOS
// coefficients b = (dp/de_int)_rho and c~ (the density-weighted acoustic speed),
// the wave strengths of the three-field decomposition are
//
//   theta       = (b/c~^2)[dE - u~ d(rho u) + (u~^2 - H~) drho]
//   alpha_+/-   = 1/2 [drho + theta +/- (d(rho u) - u~ drho)/c~]
//   alpha_0     = -theta
//
// and the dissipation is sum_k |lambda_k| alpha_k r_k with
// r_-/+ = (1, u~ -/+ c~, H~ -/+ u~ c~) and r_0 = (1, u~, H~ - c~^2/b).
// This is a local characteristic average, not a claim of exact general-EOS Roe
// Property U; any inadmissible average falls back face-locally to Rusanov.
void build_roe_flux(const Grid& grid, const MixtureFaceArrays& left,
                    const MixtureFaceArrays& right, const Work& a_face,
                    Work& output) {
    const auto sz = state_size(grid);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const arma::uword k_rho = arma::sub2ind(sz, i, mix::RHO);
        const arma::uword k_mom = arma::sub2ind(sz, i, mix::MOM);
        const arma::uword k_e   = arma::sub2ind(sz, i, mix::ENERGY);

        const double rho_l = left.conserved(k_rho), rho_r = right.conserved(k_rho);
        const double mom_l = left.conserved(k_mom), mom_r = right.conserved(k_mom);
        const double energy_l = left.conserved(k_e), energy_r = right.conserved(k_e);
        const double p_l = left.pressure(i), p_r = right.pressure(i);
        const double c_l = left.sound_speed(i), c_r = right.sound_speed(i);
        const double b_l = left.dp_deint(i), b_r = right.dp_deint(i);
        if (!(rho_l > 0.0) || !(rho_r > 0.0) || !(p_l > 0.0) || !(p_r > 0.0)
            || !(c_l > 0.0) || !(c_r > 0.0) || !(b_l > 0.0) || !(b_r > 0.0)
            || !std::isfinite(energy_l) || !std::isfinite(energy_r)) {
            rusanov_face(grid, left, right, i, a_face(i), output);
            continue;
        }

        const double v_l = mom_l/rho_l, v_r = mom_r/rho_r;
        const double h_l = (energy_l + p_l)/rho_l;
        const double h_r = (energy_r + p_r)/rho_r;
        const double sr_l = std::sqrt(rho_l), sr_r = std::sqrt(rho_r);
        const double denom = sr_l + sr_r;
        const double velocity = (sr_l*v_l + sr_r*v_r)/denom;
        const double enthalpy = (sr_l*h_l + sr_r*h_r)/denom;
        const double b = (sr_l*b_l + sr_r*b_r)/denom;
        const double c2 = (sr_l*c_l*c_l + sr_r*c_r*c_r)/denom;
        if (!(b > 0.0) || !(c2 > 0.0) || !std::isfinite(enthalpy)) {
            rusanov_face(grid, left, right, i, a_face(i), output);
            continue;
        }
        const double c = std::sqrt(c2);

        const double d_rho = rho_r - rho_l;
        const double d_mom = mom_r - mom_l;
        const double d_energy = energy_r - energy_l;
        const double theta = b/c2*(d_energy - velocity*d_mom
                                   + (velocity*velocity - enthalpy)*d_rho);
        const double acoustic = (d_mom - velocity*d_rho)/c;
        const double alpha_minus = 0.5*(d_rho + theta - acoustic);
        const double alpha_zero = -theta;
        const double alpha_plus = 0.5*(d_rho + theta + acoustic);
        const double lambda_minus = std::abs(velocity - c);
        const double lambda_zero = std::abs(velocity);
        const double lambda_plus = std::abs(velocity + c);

        const double dissipation[num_of_mixture_eq] = {
            lambda_minus*alpha_minus + lambda_zero*alpha_zero
                + lambda_plus*alpha_plus,
            lambda_minus*alpha_minus*(velocity - c) + lambda_zero*alpha_zero*velocity
                + lambda_plus*alpha_plus*(velocity + c),
            lambda_minus*alpha_minus*(enthalpy - velocity*c)
                + lambda_zero*alpha_zero*(enthalpy - c2/b)
                + lambda_plus*alpha_plus*(enthalpy + velocity*c)};
        if (!std::isfinite(dissipation[0]) || !std::isfinite(dissipation[1])
            || !std::isfinite(dissipation[2])) {
            rusanov_face(grid, left, right, i, a_face(i), output);
            continue;
        }
        for (arma::uword row = 0; row < num_of_mixture_eq; ++row) {
            const arma::uword k = arma::sub2ind(sz, i, row);
            output(k) = 0.5*(left.flux(k) + right.flux(k) - dissipation[row]);
        }
    }
}

// Geometric momentum source: the flux-tube pressure term p d_s(ln A) written in
// the code's B form (A B = const => d_s(ln A) = -d_s(ln B) = B d_s(1/B)) plus
// the gravitational body force -rho d_s(phi_g).
void build_source(const Grid& grid, const Vec& state,
                  const MixtureField& decoded, Work& source) {
    source.zeros(grid.ns*num_of_mixture_eq);
    const auto sz = state_size(grid);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const double g_par = -static_cast<double>(
            (grid.phi_g_iph(i) - grid.phi_g_imh(i))/grid.ds_i(i));
        source(arma::sub2ind(sz, i, mix::MOM)) =
            decoded.cells[i].p * static_cast<double>(grid.B_i(i))
                * static_cast<double>(grid.dinvB_ds_i(i))
            + static_cast<double>(state(arma::sub2ind(sz, i, mix::RHO)))*g_par;
    }
}

// Read-only copy-out of the production face reconstruction and the numerical
// mass flux actually used at the i+1/2 faces. Reads the same arrays the
// continuity row is differenced from; writes nothing the solver consumes.
void capture_face_flux(const Grid& grid, const Work& w, const Work& wl,
                       const Work& wr, const MixtureFaceArrays& fl,
                       const MixtureFaceArrays& fr, const Work& a_face,
                       const Work& flux, const Work& r, const Work& r_ip1,
                       const Work& lp_r, const Work& lm_rip1) {
    MixtureFaceFluxCapture& c = grid.face_flux_capture;
    c.resize(grid.ns);
    const auto sz = state_size(grid);
    const bool has_eq = grid.eq_wb && !grid.eq_residual.is_empty();
    for (arma::uword i = 0; i < grid.ns; ++i) {
        auto W = [&](const Work& v, arma::uword slot) {
            return v(arma::sub2ind(sz, i, slot));
        };
        c.rho_cell[i] = std::exp(W(w, MIX_LOG_RHO));
        c.v_cell[i]   = W(w, MIX_U);
        c.rho_L[i]    = std::exp(W(wl, MIX_LOG_RHO));
        c.rho_R[i]    = std::exp(W(wr, MIX_LOG_RHO));
        c.v_L[i]      = W(wl, MIX_U);
        c.v_R[i]      = W(wr, MIX_U);
        // Both (T, p) pairs, whichever slot 2 actually carried. The missing member
        // of the pair comes from the SAME authoritative closure the solver used for
        // this face, so the mismatch diagnostic is directly comparable between the
        // two reconstructions. Capture steps only.
        auto thermal_pair = [&](const Work& v, double rho, double& temperature,
                                double& pressure) {
            const double slot = W(v, MIX_THERMAL);
            if (grid.pressure_reconstruct) {
                pressure = std::exp(slot);
                temperature = equilibrium_temperature_from_density_pressure(
                    rho, pressure);
            } else {
                temperature = std::exp(slot);
                pressure = equilibrium_caloric_state(rho, temperature).pressure;
            }
        };
        thermal_pair(w,  c.rho_cell[i], c.T_cell[i], c.p_cell[i]);
        thermal_pair(wl, c.rho_L[i],    c.T_L[i],    c.p_L[i]);
        thermal_pair(wr, c.rho_R[i],    c.T_R[i],    c.p_R[i]);
        c.cs_L[i]   = fl.sound_speed(i);
        c.cs_R[i]   = fr.sound_speed(i);
        c.a_face[i] = a_face(i);
        const arma::uword k_rho = arma::sub2ind(sz, i, mix::RHO);
        c.f_central[i] = 0.5*(fl.flux(k_rho) + fr.flux(k_rho));
        c.f_total[i]   = flux(k_rho);
        c.f_diff[i]    = grid.roe_characteristic_flux
            ? c.f_total[i] - c.f_central[i]
            : -0.5*c.a_face[i]*(fr.conserved(k_rho) - fl.conserved(k_rho));
        c.eq_residual_mass[i] = has_eq
            ? static_cast<double>(grid.eq_residual(k_rho)) : 0.0;
        c.r_rho[i]         = W(r,       MIX_LOG_RHO);
        c.phi_plus_rho[i]  = W(lp_r,    MIX_LOG_RHO);
        c.r_ip1_rho[i]     = W(r_ip1,   MIX_LOG_RHO);
        c.phi_minus_rho[i] = W(lm_rip1, MIX_LOG_RHO);
        c.r_v[i]           = W(r,       MIX_U);
        c.phi_plus_v[i]    = W(lp_r,    MIX_U);
        // Slot-2 limiter inputs/outputs: log(p) in the release reconstruction.
        c.r_T[i]           = W(r,       MIX_THERMAL);
        c.phi_plus_T[i]    = W(lp_r,    MIX_THERMAL);
    }
    c.valid = true;
}

// Predictor decode: caloric only (the half step never needs the acoustic Gamma1).
void decode_predictor_into(const Grid& grid, const Vec& state,
                           const MixtureField& previous, MixtureField& output) {
    if (state.n_elem != grid.n_mixture_state)
        throw std::invalid_argument("mixture predictor state size mismatch");
    output.cells.clear();
    const arma::uword ext_n = grid.ns+4;
    output.extended_primitive.set_size(ext_n*num_of_mixture_eq);
    output.extended_temperature.set_size(ext_n);
    const auto sz = state_size(grid);
    const auto ext_size = arma::size(ext_n, num_of_mixture_eq);
    const bool log_p = grid.pressure_reconstruct;
    auto put = [&](arma::uword ext_i, const CaloricMixtureThermo& th,
                   double momentum) {
        output.extended_primitive(arma::sub2ind(ext_size, ext_i, MIX_LOG_RHO)) =
            std::log(th.rho);
        output.extended_primitive(arma::sub2ind(ext_size, ext_i, MIX_U)) =
            momentum/th.rho;
        output.extended_primitive(arma::sub2ind(ext_size, ext_i, MIX_THERMAL)) =
            log_p ? std::log(th.p) : std::log(th.T);
        output.extended_temperature(ext_i) = th.T;
    };
    parallel_for_cells(grid.ns, [&](std::size_t raw_i) {
        const arma::uword i = static_cast<arma::uword>(raw_i);
        auto at = [&](arma::uword row) {
            return static_cast<double>(state(arma::sub2ind(sz, i, row)));
        };
        const CaloricMixtureThermo th = decode_equilibrium_caloric_mixture(
            grid.eos_gamma_table, at(mix::RHO), at(mix::MOM), at(mix::ENERGY),
            mixture_cell_phi(grid, i), previous.cells[i].T,
            grid.eos_gamma_debug_clamp);
        grid.store_eos_temperature_hint(i, th.T);
        put(i+2, th, at(mix::MOM));
    });
    auto decode_ghost = [&](const Vec& ghost, double phi, arma::uword ext_i) {
        const CaloricMixtureThermo th = decode_equilibrium_caloric_mixture(
            grid.eos_gamma_table, ghost(mix::RHO), ghost(mix::MOM),
            ghost(mix::ENERGY), phi, std::numeric_limits<double>::quiet_NaN(),
            grid.eos_gamma_debug_clamp);
        put(ext_i, th, static_cast<double>(ghost(mix::MOM)));
    };
    const double phi_inner = grid.phi_g_imh(0);
    const double phi_outer = grid.phi_g_iph(grid.ns-1);
    decode_ghost(grid.mix_inner_boundary1, phi_inner, 0);
    decode_ghost(grid.mix_inner_boundary0, phi_inner, 1);
    decode_ghost(grid.mix_outer_boundary0, phi_outer, grid.ns+2);
    decode_ghost(grid.mix_outer_boundary1, phi_outer, grid.ns+3);
}

} // namespace

Vec mixture_rhs_explicit(const Grid& grid, const Vec& state,
                         const MixtureField& decoded, double dt_predictor) {
    ProfileScope timer(ProfileRegion::Rhs);
    decoded.require_matches(grid, state);
    MixtureRhsScratch& scratch = grid.mixture_rhs_scratch;
    auto& m = scratch.work;
    Work& w=m[0]; Work& w_ip1=m[1]; Work& w_im1=m[2];
    Work& w_ip2=m[3]; Work& w_im2=m[4];
    Work& dw_iph=m[5]; Work& dw_imh=m[6];
    Work& r=m[7]; Work& r_ip1=m[8]; Work& r_im1=m[9];
    Work& t_hint_0=m[10]; Work& t_hint_ip1=m[11]; Work& t_hint_im1=m[12];
    Work& wr_iph=m[13]; Work& wl_iph=m[14];
    Work& wr_imh=m[15]; Work& wl_imh=m[16];
    Work& wt=m[17]; Work& wt_ip1=m[18]; Work& wt_im1=m[19];
    Work& lp_r=m[20]; Work& lm_r=m[21];
    Work& lm_rip1=m[22]; Work& lp_rim1=m[23];

    // Parent-cell temperature hints for the (rho,p)->T face inversion. Built only
    // when the pressure reconstruction is active; the log(T) path never inverts.
    const Work* h0 = nullptr;
    const Work* h_ip1 = nullptr;
    const Work* h_im1 = nullptr;
    if (grid.pressure_reconstruct) {
        temperature_stencil_view(grid, decoded,  0, t_hint_0);
        temperature_stencil_view(grid, decoded,  1, t_hint_ip1);
        temperature_stencil_view(grid, decoded, -1, t_hint_im1);
        h0 = &t_hint_0; h_ip1 = &t_hint_ip1; h_im1 = &t_hint_im1;
    }
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

    // Non-uniform reconstruction weights and slope-ratio metrics depend only on
    // the static mesh, so they are built once per mesh generation.
    const Work& W1 = scratch.W1;
    const Work& W3 = scratch.W3;
    const Work& W4 = scratch.W4;
    if (!grid.uniform_mesh) {
        if (!scratch.weights_valid
            || scratch.weights_generation != grid.metrics_generation()
            || scratch.W1.n_elem != grid.ns*num_of_mixture_eq) {
            const Vec w1_i = 0.5f*grid.ds_i/grid.ds_iph_i;
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
    lim_plus(r, lp_r);       lim_minus(r, lm_r);
    lim_minus(r_ip1, lm_rip1); lim_plus(r_im1, lp_rim1);

    // First (predictor) pass. The predictor update differences ONLY
    // fl_iph.flux - fr_imh.flux, so only those two reconstructions are needed.
    if (grid.uniform_mesh) {
        wl_iph = w + 0.5*lp_r%(w_ip1 - w);
        wr_imh = w - 0.5*lm_r%(w_ip1 - w);
    } else {
        wl_iph = w + W1%lp_r%(w_ip1 - w);
        wr_imh = w - W1%lm_r%(w_ip1 - w);
    }

    MixtureFaceArrays& fr_iph = scratch.faces[0];
    MixtureFaceArrays& fl_iph = scratch.faces[1];
    MixtureFaceArrays& fr_imh = scratch.faces[2];
    MixtureFaceArrays& fl_imh = scratch.faces[3];
    {
        // Both predictor faces are extrapolated from cell i, so cell i's
        // temperature is the natural Newton seed for both.
        const FaceRequest requests[2] = {
            {&wl_iph, &grid.phi_g_iph, &fl_iph, h0},
            {&wr_imh, &grid.phi_g_imh, &fr_imh, h0}};
        build_faces(grid, requests, 2);
    }

    const auto sz = state_size(grid);
    Vec& predicted = scratch.predicted_state;
    predicted.set_size(grid.n_mixture_state);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const double scale = dt_predictor/static_cast<double>(grid.ds_i(i));
        for (arma::uword row = 0; row < num_of_mixture_eq; ++row) {
            const arma::uword k = arma::sub2ind(sz, i, row);
            predicted(k) = static_cast<float>(
                static_cast<double>(state(k))
                - scale*(fl_iph.flux(k) - fr_imh.flux(k)));
        }
    }
    MixtureField& decoded_predicted = scratch.predicted;
    decode_predictor_into(grid, predicted, decoded, decoded_predicted);
    stencil_view(grid, decoded_predicted,  0, wt);
    stencil_view(grid, decoded_predicted,  1, wt_ip1);
    stencil_view(grid, decoded_predicted, -1, wt_im1);

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

    // Corrector parentage, matching the four reconstructions just above:
    //   wr_iph[i] from cell i+1, wl_iph[i] and wr_imh[i] from cell i,
    //   wl_imh[i] from cell i-1.
    {
        const FaceRequest requests[4] = {
            {&wr_iph, &grid.phi_g_iph, &fr_iph, h_ip1},
            {&wl_iph, &grid.phi_g_iph, &fl_iph, h0},
            {&wr_imh, &grid.phi_g_imh, &fr_imh, h0},
            {&wl_imh, &grid.phi_g_imh, &fl_imh, h_im1}};
        build_faces(grid, requests, 4);
    }

    Work a_imh(grid.ns), a_iph(grid.ns);
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
    flux_imh.set_size(grid.n_mixture_state);
    flux_iph.set_size(grid.n_mixture_state);
    if (grid.roe_characteristic_flux) {
        build_roe_flux(grid, fl_imh, fr_imh, a_imh, flux_imh);
        build_roe_flux(grid, fl_iph, fr_iph, a_iph, flux_iph);
    } else {
        for (arma::uword i = 0; i < grid.ns; ++i) {
            rusanov_face(grid, fl_imh, fr_imh, i, a_imh(i), flux_imh);
            rusanov_face(grid, fl_iph, fr_iph, i, a_iph(i), flux_iph);
        }
    }

    Work& source = scratch.source;
    build_source(grid, state, decoded, source);
    Vec& rhs = scratch.rhs;
    rhs.set_size(grid.n_mixture_state);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const double B_i = grid.B_i(i);
        const double inv_B_iph = 1.0/static_cast<double>(grid.B_iph(i));
        const double inv_B_imh = 1.0/static_cast<double>(grid.B_imh(i));
        const double inv_ds = 1.0/static_cast<double>(grid.ds_i(i));
        for (arma::uword row = 0; row < num_of_mixture_eq; ++row) {
            const arma::uword k = arma::sub2ind(sz, i, row);
            rhs(k) = static_cast<float>(
                -B_i*(flux_iph(k)*inv_B_iph - flux_imh(k)*inv_B_imh)*inv_ds
                + source(k));
        }
    }
    if (grid.eq_wb && !grid.eq_residual.is_empty()) rhs -= grid.eq_residual;
    if (grid.capture_face_flux)
        capture_face_flux(grid, w, wl_iph, wr_iph, fl_iph, fr_iph, a_iph,
                          flux_iph, r, r_ip1, lp_r, lm_rip1);
    return Vec(rhs);
}

// ============================================================================
// Timestep
// ============================================================================

// The release has no volumetric energy source, so the acoustic CFL condition is
// the only timestep constraint: the implicit conduction stage is unconditionally
// stable and imposes none.
Vec mixture_timestep(const Grid& grid, const Vec& state,
                     const MixtureField& decoded) {
    ProfileScope timer(ProfileRegion::Cfl);
    decoded.require_matches(grid, state);
    const auto sz = state_size(grid);
    Vec dt_i(grid.ns);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const MixtureThermo& th = decoded.cells[i];
        const double momentum =
            static_cast<double>(state(arma::sub2ind(sz, i, mix::MOM)));
        const double speed = std::abs(momentum/th.rho)
                           + std::sqrt(th.gamma1*th.p/th.rho);
        dt_i(i) = static_cast<float>(grid.CFL*grid.ds_i(i)/speed);
    }
    profile_note_timestep_limiter(TimestepLimiter::Acoustic);
    return arma::min(dt_i)*arma::ones<Vec>(grid.ns);
}

} // namespace chromosphere
