#include "chromosphere.hpp"
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

void DecodedMixtureField::require_matches(const Grid& grid,
                                          const Vec& state) const {
    if (source_grid != &grid || source_state != &state
        || source_memory != state.memptr()
        || source_elements != state.n_elem)
        throw std::logic_error(
            "decoded mixture cache does not match the immutable conserved state");
    const Vec* boundaries[] = {&grid.inner_boundary0_i,&grid.inner_boundary1_i,
                               &grid.outer_boundary0_i,&grid.outer_boundary1_i};
    for (arma::uword block=0; block<4; ++block)
        for (arma::uword row=0; row<num_of_eq; ++row)
            if (boundary_signature[block*num_of_eq+row] != (*boundaries[block])(row))
                throw std::logic_error(
                    "decoded mixture cache was invalidated by a boundary update");
}

DecodedMixtureField decode_mixture_field(
    const Grid& grid, const Vec& state, std::uint64_t generation,
    const DecodedMixtureField* previous) {
    DecodedMixtureField out;
    decode_mixture_field_into(grid,state,out,generation,previous);
    return out;
}

void decode_mixture_field_into(
    const Grid& grid, const Vec& state, DecodedMixtureField& out,
    std::uint64_t generation, const DecodedMixtureField* previous) {
    ProfileScope timer(ProfileRegion::Decode);
    if (grid.eos_gamma_table.empty())
        throw std::logic_error("decoded mixture field requires a Gamma1 table");
    if (state.n_elem != grid.n_state)
        throw std::invalid_argument("decoded mixture field state size mismatch");

    out.source_grid = &grid;
    out.source_state = &state;
    out.source_memory = state.memptr();
    out.source_elements = state.n_elem;
    out.state_generation = generation;
    const Vec* boundaries[] = {&grid.inner_boundary0_i,&grid.inner_boundary1_i,
                               &grid.outer_boundary0_i,&grid.outer_boundary1_i};
    for (arma::uword block=0; block<4; ++block)
        for (arma::uword row=0; row<num_of_eq; ++row)
            out.boundary_signature[block*num_of_eq+row] = (*boundaries[block])(row);
    out.cells.resize(grid.ns);
    const arma::uword ext_n = grid.ns+4;
    out.extended_primitive.set_size(ext_n*3);
    const auto state_size = arma::size(grid.ns, num_of_eq);
    const auto ext_size = arma::size(ext_n, 3);
    // Slot 2 carries the THERMAL reconstruction variable: log(T) by default, or
    // log(p_total) when grid.pressure_reconstruct selects the pressure-based
    // primitive set. p_total = p_i + p_n is the same total the face builder and
    // the momentum flux use (the electron pressure is already inside p_i).
    const bool log_p = grid.pressure_reconstruct;
    auto put = [&](arma::uword ext_i, const MixtureThermo& th,
                   double momentum) {
        out.extended_primitive(arma::sub2ind(ext_size, ext_i, 0)) = std::log(th.rho);
        out.extended_primitive(arma::sub2ind(ext_size, ext_i, 1)) = momentum/th.rho;
        out.extended_primitive(arma::sub2ind(ext_size, ext_i, 2)) =
            log_p ? std::log(th.p_i + th.p_n) : std::log(th.T);
    };
    const bool has_guesses = previous && previous->source_grid == &grid
        && previous->cells.size() == grid.ns;
    parallel_for_cells(grid.ns, [&](std::size_t raw_i) {
        const arma::uword i = static_cast<arma::uword>(raw_i);
        auto at = [&](arma::uword row) {
            return static_cast<double>(state(arma::sub2ind(state_size, i, row)));
        };
        const double phi = 0.5*static_cast<double>(
            grid.phi_g_imh(i)+grid.phi_g_iph(i));
        const double guess = has_guesses ? previous->cells[i].T
            : std::numeric_limits<double>::quiet_NaN();
        out.cells[i] = decode_equilibrium_mixture(
            grid.eos_gamma_table, at(cons::RHO_I), at(cons::RHO_N),
            at(cons::MOM_I), at(cons::MOM_N), at(cons::E_I), at(cons::E_N),
            phi, guess, grid.eos_gamma_debug_clamp);
        // Refresh the shared hint, but do NOT read it here: this site has its own
        // guess (`previous`), validated against source_grid and cell count. Each
        // worker writes only its own pre-sized hint slot.
        grid.store_eos_temperature_hint(i, out.cells[i].T);
        put(i+2, out.cells[i], at(cons::MOM_I)+at(cons::MOM_N));
    });

    if (USE_NEUMANN_BC) {
        for (arma::uword slot = 0; slot < 3; ++slot) {
            out.extended_primitive(arma::sub2ind(ext_size, 0, slot)) =
                out.extended_primitive(arma::sub2ind(ext_size, 2, slot));
            out.extended_primitive(arma::sub2ind(ext_size, 1, slot)) =
                out.extended_primitive(arma::sub2ind(ext_size, 2, slot));
            out.extended_primitive(arma::sub2ind(ext_size, grid.ns+2, slot)) =
                out.extended_primitive(arma::sub2ind(ext_size, grid.ns+1, slot));
            out.extended_primitive(arma::sub2ind(ext_size, grid.ns+3, slot)) =
                out.extended_primitive(arma::sub2ind(ext_size, grid.ns+1, slot));
        }
    } else {
        auto decode_ghost = [&](const Vec& ghost, double phi,
                                arma::uword ext_i) {
            const MixtureThermo th = decode_equilibrium_mixture(
                grid.eos_gamma_table, ghost(cons::RHO_I), ghost(cons::RHO_N),
                ghost(cons::MOM_I), ghost(cons::MOM_N), ghost(cons::E_I),
                ghost(cons::E_N), phi,
                std::numeric_limits<double>::quiet_NaN(),
                grid.eos_gamma_debug_clamp);
            put(ext_i, th, static_cast<double>(ghost(cons::MOM_I)+ghost(cons::MOM_N)));
        };
        const double phi_inner = grid.phi_g_imh(0);
        const double phi_outer = grid.phi_g_iph(grid.ns-1);
        decode_ghost(grid.inner_boundary1_i, phi_inner, 0);
        decode_ghost(grid.inner_boundary0_i, phi_inner, 1);
        decode_ghost(grid.outer_boundary0_i, phi_outer, grid.ns+2);
        decode_ghost(grid.outer_boundary1_i, phi_outer, grid.ns+3);
    }
}

namespace {

// Slot 2 is the thermal reconstruction variable: log(T) with the default
// reconstruction, log(p_total) when Grid::pressure_reconstruct is set. The
// packing, limiter arithmetic and extrapolation are identical either way; only
// decode (which quantity is written) and build_mixture_face_cell (how the face
// state is closed) branch on the flag.
enum MixtureSlot : arma::uword { MIX_LOG_RHO = 0, MIX_V = 1, MIX_THERMAL = 2 };

struct MixtureFaceBundle {
    Vec& conserved;
    Vec& flux;
    Vec& spectral_radius;
    Vec& dp_deint_rho;
};

using MixtureVec = arma::Col<double>;

// The gamma reconstruction limits exactly THREE variables — log(rho), the
// mixture velocity and log(T) — but its scratch used to be packed at the
// num_of_eq-row conserved-state stride, so every difference, ratio, limiter and
// extrapolation moved 7 rows per cell and discarded 4 of them. All of these
// operations are element-wise, so shrinking the stride to 3 leaves every
// surviving element bit-for-bit unchanged while cutting the reconstruction
// arithmetic and its scratch footprint by 7/3.
constexpr arma::uword MIX_ROWS = 3;

inline arma::SizeMat mixture_size(const Grid& grid) {
    return arma::size(grid.ns, MIX_ROWS);
}

void mixture_stencil_view(const Grid& grid,
                          const DecodedMixtureField& decoded,
                          int offset, MixtureVec& result) {
    result.set_size(grid.ns*MIX_ROWS);
    const arma::uword ext_n = grid.ns+4;
    const auto ext_size = arma::size(ext_n, 3);
    const auto packed_size = mixture_size(grid);
    for (arma::uword slot = 0; slot < MIX_ROWS; ++slot) {
        for (arma::uword i = 0; i < grid.ns; ++i) {
            const arma::uword ext_i = static_cast<arma::uword>(
                static_cast<int>(i)+2+offset);
            result(arma::sub2ind(packed_size, i, slot)) =
                decoded.extended_primitive(arma::sub2ind(ext_size, ext_i, slot));
        }
    }
}

void decode_predicted_caloric_primitives_into(
    const Grid& grid, const Vec& state, const DecodedMixtureField& previous,
    DecodedMixtureField& output) {
    if (state.n_elem != grid.n_state)
        throw std::invalid_argument("predicted caloric decode state size mismatch");
    output.cells.clear();
    const arma::uword ext_n = grid.ns+4;
    output.extended_primitive.set_size(ext_n*MIX_ROWS);
    const auto state_size = arma::size(grid.ns, num_of_eq);
    const auto ext_size = arma::size(ext_n, MIX_ROWS);
    // Same thermal-slot convention as decode_mixture_field_into: the predictor
    // must carry the SAME variable the corrector reconstructs.
    const bool log_p = grid.pressure_reconstruct;
    auto put = [&](arma::uword ext_i, const CaloricMixtureThermo& th,
                   double momentum) {
        output.extended_primitive(arma::sub2ind(ext_size, ext_i, MIX_LOG_RHO)) =
            std::log(th.rho);
        output.extended_primitive(arma::sub2ind(ext_size, ext_i, MIX_V)) =
            momentum/th.rho;
        output.extended_primitive(arma::sub2ind(ext_size, ext_i, MIX_THERMAL)) =
            log_p ? std::log(th.p_i + th.p_n) : std::log(th.T);
    };
    parallel_for_cells(grid.ns, [&](std::size_t raw_i) {
        const arma::uword i = static_cast<arma::uword>(raw_i);
        auto at = [&](arma::uword row) {
            return static_cast<double>(state(arma::sub2ind(state_size, i, row)));
        };
        const double phi = 0.5*static_cast<double>(
            grid.phi_g_imh(i)+grid.phi_g_iph(i));
        const double guess = previous.cells[i].T;
        const CaloricMixtureThermo th = decode_equilibrium_caloric_mixture(
            grid.eos_gamma_table, at(cons::RHO_I), at(cons::RHO_N),
            at(cons::MOM_I), at(cons::MOM_N), at(cons::E_I), at(cons::E_N),
            phi, guess, grid.eos_gamma_debug_clamp);
        grid.store_eos_temperature_hint(i, th.T);
        put(i+2, th, at(cons::MOM_I)+at(cons::MOM_N));
    });
    if (USE_NEUMANN_BC) {
        for (arma::uword slot = 0; slot < MIX_ROWS; ++slot) {
            output.extended_primitive(arma::sub2ind(ext_size,0,slot)) =
                output.extended_primitive(arma::sub2ind(ext_size,2,slot));
            output.extended_primitive(arma::sub2ind(ext_size,1,slot)) =
                output.extended_primitive(arma::sub2ind(ext_size,2,slot));
            output.extended_primitive(arma::sub2ind(ext_size,grid.ns+2,slot)) =
                output.extended_primitive(arma::sub2ind(ext_size,grid.ns+1,slot));
            output.extended_primitive(arma::sub2ind(ext_size,grid.ns+3,slot)) =
                output.extended_primitive(arma::sub2ind(ext_size,grid.ns+1,slot));
        }
    } else {
        auto decode_ghost = [&](const Vec& ghost, double phi, arma::uword ext_i) {
            const CaloricMixtureThermo th = decode_equilibrium_caloric_mixture(
                grid.eos_gamma_table, ghost(cons::RHO_I), ghost(cons::RHO_N),
                ghost(cons::MOM_I), ghost(cons::MOM_N), ghost(cons::E_I),
                ghost(cons::E_N), phi, std::numeric_limits<double>::quiet_NaN(),
                grid.eos_gamma_debug_clamp);
            put(ext_i, th, static_cast<double>(
                ghost(cons::MOM_I)+ghost(cons::MOM_N)));
        };
        const double phi_inner = grid.phi_g_imh(0);
        const double phi_outer = grid.phi_g_iph(grid.ns-1);
        decode_ghost(grid.inner_boundary1_i,phi_inner,0);
        decode_ghost(grid.inner_boundary0_i,phi_inner,1);
        decode_ghost(grid.outer_boundary0_i,phi_outer,grid.ns+2);
        decode_ghost(grid.outer_boundary1_i,phi_outer,grid.ns+3);
    }
}

void bcast_mixture(const Grid& grid, const Vec& values, MixtureVec& packed) {
    packed.set_size(grid.ns*MIX_ROWS);
    for (arma::uword k = 0; k < MIX_ROWS; ++k) {
        for (arma::uword i = 0; i < grid.ns; ++i)
            packed(arma::sub2ind(mixture_size(grid), i, k)) = values(i);
    }
}

void mixture_minmod(const MixtureVec& ratio, MixtureVec& result) {
    result.set_size(ratio.n_elem);
    for (arma::uword i = 0; i < ratio.n_elem; ++i) {
        const double r = ratio(i);
        result(i) = std::isfinite(r) ? std::max(0.0, std::min(1.0, r)) : 0.0;
    }
}

void mixture_mc3(const MixtureVec& ratio, double beta, bool plus,
                 MixtureVec& result) {
    result.set_size(ratio.n_elem);
    for (arma::uword i = 0; i < ratio.n_elem; ++i) {
        const double r = ratio(i);
        if (!std::isfinite(r)) {
            result(i) = 0.0;
            continue;
        }
        const double third_order = plus ? (2.0*r + 1.0)/3.0 : (r + 2.0)/3.0;
        result(i) = std::max(0.0, std::min(std::min(beta*r, beta), third_order));
    }
}

void prepare_mixture_face(const Grid& grid, MixtureFaceBundle out) {
    // Every (i,k) slot is written by one physical-cell worker. Pre-size once,
    // outside the OpenMP region; no worker resizes shared Armadillo storage.
    out.conserved.set_size(grid.n_state);
    out.flux.set_size(grid.n_state);
    out.spectral_radius.set_size(grid.n_state);
    out.dp_deint_rho.set_size(grid.ns);
}

void build_mixture_face_cell(const Grid& grid, const MixtureVec& mixture,
                             const Vec& phi_face, MixtureFaceBundle out,
                             arma::uword i) {
    const auto sz = arma::size(grid.ns, num_of_eq);
    const auto msz = mixture_size(grid);
    const double log_rho = mixture(arma::sub2ind(msz, i, MIX_LOG_RHO));
    const double velocity = mixture(arma::sub2ind(msz, i, MIX_V));
    const double thermal = mixture(arma::sub2ind(msz, i, MIX_THERMAL));
    const MixtureFaceState face = grid.pressure_reconstruct
        ? equilibrium_mixture_face_state_from_log_pressure(
              grid.eos_gamma_table, log_rho, velocity, thermal,
              static_cast<double>(phi_face(i)), grid.eos_trace_fraction_floor,
              grid.eos_gamma_debug_clamp)
        : equilibrium_mixture_face_state_from_logs(
              grid.eos_gamma_table, log_rho, velocity, thermal,
              static_cast<double>(phi_face(i)), grid.eos_trace_fraction_floor,
              grid.eos_gamma_debug_clamp);
    const ProjectedMixture& u = face.conserved;
    const std::array<double, 7> flux = equilibrium_mixture_flux(face);
    const double values[7] = {u.rho_i, u.rho_n, u.momentum_i, u.momentum_n,
                              u.energy_i, u.energy_n, u.energy_e};
    const float a = static_cast<float>(std::abs(velocity) + face.sound_speed);
    out.dp_deint_rho(i) = static_cast<float>(face.dp_deint_rho);
    for (arma::uword k = 0; k < num_of_eq; ++k) {
        out.conserved(arma::sub2ind(sz, i, k)) = static_cast<float>(values[k]);
        out.flux(arma::sub2ind(sz, i, k)) = static_cast<float>(flux[k]);
        out.spectral_radius(arma::sub2ind(sz, i, k)) = a;
    }
}

void build_mixture_face_batch2(
    const Grid& grid,
    const MixtureVec& mixture0, const Vec& phi0, MixtureFaceBundle out0,
    const MixtureVec& mixture1, const Vec& phi1, MixtureFaceBundle out1) {
    prepare_mixture_face(grid, out0);
    prepare_mixture_face(grid, out1);
    parallel_for_cells(grid.ns, [&](std::size_t raw_i) {
        const arma::uword i = static_cast<arma::uword>(raw_i);
        build_mixture_face_cell(grid, mixture0, phi0, out0, i);
        build_mixture_face_cell(grid, mixture1, phi1, out1, i);
    });
}

void build_mixture_face_batch4(
    const Grid& grid,
    const MixtureVec& mixture0, const Vec& phi0, MixtureFaceBundle out0,
    const MixtureVec& mixture1, const Vec& phi1, MixtureFaceBundle out1,
    const MixtureVec& mixture2, const Vec& phi2, MixtureFaceBundle out2,
    const MixtureVec& mixture3, const Vec& phi3, MixtureFaceBundle out3) {
    prepare_mixture_face(grid, out0);
    prepare_mixture_face(grid, out1);
    prepare_mixture_face(grid, out2);
    prepare_mixture_face(grid, out3);
    parallel_for_cells(grid.ns, [&](std::size_t raw_i) {
        const arma::uword i = static_cast<arma::uword>(raw_i);
        build_mixture_face_cell(grid, mixture0, phi0, out0, i);
        build_mixture_face_cell(grid, mixture1, phi1, out1, i);
        build_mixture_face_cell(grid, mixture2, phi2, out2, i);
        build_mixture_face_cell(grid, mixture3, phi3, out3, i);
    });
}

// Roe-local corrector flux for the three authoritative equilibrium-manifold
// totals (rho, rho*v, E_i+E_n). This is deliberately a local characteristic
// average, not a claim of exact general-EOS Roe Property U. The seven carrier
// rows are only a storage representation: after the explicit update the active
// integrator projects the three totals back to the Saha equilibrium manifold.
// Any invalid characteristic state falls back face-locally to the validated
// Rusanov flux.
void build_roe_local_mixture_flux(const Grid& grid,
                                  const MixtureFaceBundle& left,
                                  const MixtureFaceBundle& right,
                                  const Vec& a_face, Vec& output) {
    output.set_size(grid.n_state);
    const auto sz = arma::size(grid.ns, num_of_eq);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        auto P = [&](const Vec& v, arma::uword row) -> double {
            return static_cast<double>(v(arma::sub2ind(sz, i, row)));
        };
        auto rusanov_face = [&]() {
            const double a = P(a_face, cons::RHO_I);
            for (arma::uword row = 0; row < num_of_eq; ++row) {
                output(arma::sub2ind(sz, i, row)) = static_cast<float>(
                    0.5 * (P(left.flux, row) + P(right.flux, row)
                           - a * (P(right.conserved, row) - P(left.conserved, row))));
            }
        };

        const double rho_l = P(left.conserved, cons::RHO_I)
                           + P(left.conserved, cons::RHO_N);
        const double rho_r = P(right.conserved, cons::RHO_I)
                           + P(right.conserved, cons::RHO_N);
        const double mom_l = P(left.conserved, cons::MOM_I)
                           + P(left.conserved, cons::MOM_N);
        const double mom_r = P(right.conserved, cons::MOM_I)
                           + P(right.conserved, cons::MOM_N);
        const double energy_l = P(left.conserved, cons::E_I)
                              + P(left.conserved, cons::E_N);
        const double energy_r = P(right.conserved, cons::E_I)
                              + P(right.conserved, cons::E_N);
        if (!(rho_l > 0.0) || !(rho_r > 0.0)
            || !std::isfinite(energy_l) || !std::isfinite(energy_r)) {
            rusanov_face();
            continue;
        }

        const double v_l = mom_l / rho_l;
        const double v_r = mom_r / rho_r;
        const double p_l = P(left.flux, cons::MOM_I) + P(left.flux, cons::MOM_N)
                         - rho_l * v_l * v_l;
        const double p_r = P(right.flux, cons::MOM_I) + P(right.flux, cons::MOM_N)
                         - rho_r * v_r * v_r;
        const double h_l = (energy_l + p_l) / rho_l;
        const double h_r = (energy_r + p_r) / rho_r;
        const double c_l = P(left.spectral_radius, cons::RHO_I) - std::abs(v_l);
        const double c_r = P(right.spectral_radius, cons::RHO_I) - std::abs(v_r);
        const double b_l = static_cast<double>(left.dp_deint_rho(i));
        const double b_r = static_cast<double>(right.dp_deint_rho(i));
        if (!(p_l > 0.0) || !(p_r > 0.0) || !(c_l > 0.0) || !(c_r > 0.0)
            || !(b_l > 0.0) || !(b_r > 0.0)
            || !std::isfinite(h_l) || !std::isfinite(h_r)) {
            rusanov_face();
            continue;
        }

        const double sr_l = std::sqrt(rho_l);
        const double sr_r = std::sqrt(rho_r);
        const double denom = sr_l + sr_r;
        const double velocity = (sr_l * v_l + sr_r * v_r) / denom;
        const double enthalpy = (sr_l * h_l + sr_r * h_r) / denom;
        const double b = (sr_l * b_l + sr_r * b_r) / denom;
        const double c2 = (sr_l * c_l * c_l + sr_r * c_r * c_r) / denom;
        if (!(b > 0.0) || !(c2 > 0.0) || !std::isfinite(enthalpy)) {
            rusanov_face();
            continue;
        }
        const double c = std::sqrt(c2);

        const double d_rho = rho_r - rho_l;
        const double d_mom = mom_r - mom_l;
        const double d_energy = energy_r - energy_l;
        const double theta = b / c2 * (d_energy - velocity * d_mom
            + (velocity * velocity - enthalpy) * d_rho);
        const double acoustic_momentum = (d_mom - velocity * d_rho) / c;
        const double alpha_minus = 0.5 * (d_rho + theta - acoustic_momentum);
        const double alpha_zero = -theta;
        const double alpha_plus = 0.5 * (d_rho + theta + acoustic_momentum);
        const double lambda_minus = std::abs(velocity - c);
        const double lambda_zero = std::abs(velocity);
        const double lambda_plus = std::abs(velocity + c);

        const double d0 = lambda_minus * alpha_minus
                        + lambda_zero * alpha_zero
                        + lambda_plus * alpha_plus;
        const double d1 = lambda_minus * alpha_minus * (velocity - c)
                        + lambda_zero * alpha_zero * velocity
                        + lambda_plus * alpha_plus * (velocity + c);
        const double d2 = lambda_minus * alpha_minus * (enthalpy - velocity * c)
                        + lambda_zero * alpha_zero * (enthalpy - c2 / b)
                        + lambda_plus * alpha_plus * (enthalpy + velocity * c);
        if (!std::isfinite(d0) || !std::isfinite(d1) || !std::isfinite(d2)) {
            rusanov_face();
            continue;
        }

        const double x_l = P(left.conserved, cons::RHO_I) / rho_l;
        const double x_r = P(right.conserved, cons::RHO_I) / rho_r;
        const double x = std::max(0.0, std::min(1.0,
            (sr_l * x_l + sr_r * x_r) / denom));
        const double dissipation[7] = {
            x * d0, (1.0 - x) * d0,
            x * d1, (1.0 - x) * d1,
            x * d2, (1.0 - x) * d2,
            0.0
        };
        for (arma::uword row = 0; row < num_of_eq; ++row) {
            output(arma::sub2ind(sz, i, row)) = static_cast<float>(
                0.5 * (P(left.flux, row) + P(right.flux, row) - dissipation[row]));
        }
    }
}

void mixture_source(const Grid& grid, const Vec& state,
                    const DecodedMixtureField& decoded, Vec& source) {
    decoded.require_matches(grid, state);
    source.zeros(grid.n_state);
    const auto sz = arma::size(grid.ns, num_of_eq);
    const Vec gb = -(grid.phi_g_iph - grid.phi_g_imh) / grid.ds_i;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        auto at = [&](arma::uword k) -> double {
            return static_cast<double>(state(arma::sub2ind(sz, i, k)));
        };
        const MixtureThermo& thermo = decoded.cells[i];
        source(arma::sub2ind(sz, i, cons::MOM_I)) = static_cast<float>(
            thermo.p_i * grid.B_i(i) * grid.dinvB_ds_i(i)
            + at(cons::RHO_I) * gb(i));
        source(arma::sub2ind(sz, i, cons::MOM_N)) = static_cast<float>(
            thermo.p_n * grid.B_i(i) * grid.dinvB_ds_i(i)
            + at(cons::RHO_N) * gb(i));
    }
}

// Diagnostic-only copy-out of the production face reconstruction and Rusanov
// TOTAL-mass flux at the i+1/2 faces (chromosphere.hpp::GammaFaceFluxCapture).
// Reads the same arrays the continuity row is differenced from; writes nothing
// the solver consumes.
void capture_gamma_face_flux(const Grid& grid, const MixtureVec& w,
                            const MixtureVec& wl, const MixtureVec& wr,
                            const MixtureFaceBundle& fl,
                            const MixtureFaceBundle& fr,
                            const Vec& a_face, const Vec& flux,
                            const MixtureVec& r, const MixtureVec& r_ip1,
                            const MixtureVec& lp_r, const MixtureVec& lm_rip1) {
    GammaFaceFluxCapture& c = grid.face_flux_capture;
    c.resize(grid.ns);
    const auto sz = arma::size(grid.ns, num_of_eq);
    const auto msz = mixture_size(grid);
    const bool has_eq = grid.eq_wb && !grid.eq_residual.is_empty();
    for (arma::uword i = 0; i < grid.ns; ++i) {
        auto W = [&](const MixtureVec& v, arma::uword slot) {
            return v(arma::sub2ind(msz, i, slot));
        };
        auto P = [&](const Vec& v, arma::uword row) {
            return static_cast<double>(v(arma::sub2ind(sz, i, row)));
        };
        c.rho_cell[i] = std::exp(W(w, MIX_LOG_RHO));
        c.v_cell[i]   = W(w, MIX_V);
        c.rho_L[i]    = std::exp(W(wl, MIX_LOG_RHO));
        c.rho_R[i]    = std::exp(W(wr, MIX_LOG_RHO));
        c.v_L[i]      = W(wl, MIX_V);
        c.v_R[i]      = W(wr, MIX_V);
        // Both (T, p_total) pairs, whichever slot 2 actually carried. The missing
        // member of the pair comes from the SAME authoritative closure the solver
        // used for this face, so the mismatch diagnostic is directly comparable
        // between the two reconstructions. Capture steps only.
        auto thermal_pair = [&](const MixtureVec& v, double rho,
                                double& temperature, double& pressure) {
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
        // Every spectral_radius row of a one-sided bundle carries |V|+c_s of that
        // state (build_mixture_face), so the sound speed is recoverable exactly.
        c.cs_L[i]     = P(fl.spectral_radius, cons::RHO_I) - std::abs(c.v_L[i]);
        c.cs_R[i]     = P(fr.spectral_radius, cons::RHO_I) - std::abs(c.v_R[i]);
        c.a_face[i]   = P(a_face, cons::RHO_I);
        const double mass_flux_L = P(fl.flux, cons::RHO_I) + P(fl.flux, cons::RHO_N);
        const double mass_flux_R = P(fr.flux, cons::RHO_I) + P(fr.flux, cons::RHO_N);
        const double rho_cons_L  = P(fl.conserved, cons::RHO_I) + P(fl.conserved, cons::RHO_N);
        const double rho_cons_R  = P(fr.conserved, cons::RHO_I) + P(fr.conserved, cons::RHO_N);
        c.f_central[i] = 0.5 * (mass_flux_L + mass_flux_R);
        c.f_total[i]   = P(flux, cons::RHO_I) + P(flux, cons::RHO_N);
        c.f_diff[i]    = grid.roe_characteristic_flux
            ? c.f_total[i] - c.f_central[i]
            : -0.5 * c.a_face[i] * (rho_cons_R - rho_cons_L);
        c.eq_residual_mass[i] = has_eq
            ? P(grid.eq_residual, cons::RHO_I) + P(grid.eq_residual, cons::RHO_N)
            : 0.0;
        c.r_rho[i]         = W(r,        MIX_LOG_RHO);
        c.phi_plus_rho[i]  = W(lp_r,     MIX_LOG_RHO);
        c.r_ip1_rho[i]     = W(r_ip1,    MIX_LOG_RHO);
        c.phi_minus_rho[i] = W(lm_rip1,  MIX_LOG_RHO);
        c.r_v[i]           = W(r,        MIX_V);
        c.phi_plus_v[i]    = W(lp_r,     MIX_V);
        // Slot-2 limiter inputs/outputs: log(T) by default, log(p) in pressure mode.
        c.r_T[i]           = W(r,        MIX_THERMAL);
        c.phi_plus_T[i]    = W(lp_r,     MIX_THERMAL);
    }
    c.valid = true;
}

Vec rhs_explicit_mixture(const Grid& grid, const Vec& xn_state,
                         const DecodedMixtureField& decoded) {
    ProfileScope timer(ProfileRegion::Rhs);
    decoded.require_matches(grid, xn_state);

    auto& m=grid.gamma_rhs_scratch.mixture;
    MixtureVec& w=m[0]; MixtureVec& w_ip1=m[1]; MixtureVec& w_im1=m[2];
    MixtureVec& w_ip2=m[3]; MixtureVec& w_im2=m[4];
    MixtureVec& dw_iph=m[5]; MixtureVec& dw_imh=m[6];
    MixtureVec& r=m[7]; MixtureVec& r_ip1=m[8]; MixtureVec& r_im1=m[9];
    MixtureVec& wr_iph=m[13]; MixtureVec& wl_iph=m[14];
    MixtureVec& wr_imh=m[15]; MixtureVec& wl_imh=m[16];
    MixtureVec& wt=m[17]; MixtureVec& wt_ip1=m[18]; MixtureVec& wt_im1=m[19];
    MixtureVec& lp_r=m[20]; MixtureVec& lm_r=m[21];
    MixtureVec& lm_rip1=m[22]; MixtureVec& lp_rim1=m[23];
    mixture_stencil_view(grid,decoded,0,w);
    mixture_stencil_view(grid,decoded,1,w_ip1);
    mixture_stencil_view(grid,decoded,-1,w_im1);
    mixture_stencil_view(grid,decoded,2,w_ip2);
    mixture_stencil_view(grid,decoded,-2,w_im2);

    dw_iph=w_ip1-w;
    dw_imh=w-w_im1;
    r=dw_imh/dw_iph;
    r_ip1=dw_iph/(w_ip2-w_ip1);
    r_im1=(w_im1-w_im2)/dw_imh;
    // Non-uniform reconstruction weights and slope-ratio metrics depend only on
    // the static mesh, so they are built once per mesh generation rather than
    // rebuilt (with six packed temporaries) on every timestep.
    GammaRhsScratch& scratch = grid.gamma_rhs_scratch;
    const MixtureVec& W1 = scratch.W1;
    const MixtureVec& W3 = scratch.W3;
    const MixtureVec& W4 = scratch.W4;
    if (!grid.uniform_mesh) {
        if (!scratch.weights_valid
            || scratch.weights_generation != grid.metrics_generation()
            || scratch.W1.n_elem != grid.ns*MIX_ROWS) {
            const Vec w1_i = 0.5f * grid.ds_i / grid.ds_iph_i;
            bcast_mixture(grid,w1_i,scratch.W1);
            bcast_mixture(grid,ip1(grid,w1_i,SLICE),scratch.W3);
            bcast_mixture(grid,im1(grid,w1_i,SLICE),scratch.W4);
            bcast_mixture(grid,grid.ds_iph_i/grid.ds_imh_i,scratch.metric_r);
            bcast_mixture(grid,ip1(grid,grid.ds_iph_i,SLICE)/grid.ds_iph_i,
                          scratch.metric_r_ip1);
            bcast_mixture(grid,grid.ds_imh_i/im1(grid,grid.ds_imh_i,SLICE),
                          scratch.metric_r_im1);
            scratch.weights_generation = grid.metrics_generation();
            scratch.weights_valid = true;
        }
        r%=scratch.metric_r;
        r_ip1%=scratch.metric_r_ip1;
        r_im1%=scratch.metric_r_im1;
    }
    r.elem(arma::find_nonfinite(r)).zeros();
    r_ip1.elem(arma::find_nonfinite(r_ip1)).zeros();
    r_im1.elem(arma::find_nonfinite(r_im1)).zeros();
    const double beta = grid.limiter_beta;
    auto lim_plus=[&](const MixtureVec& q,MixtureVec& out) {
        if (grid.mc3_limiter) mixture_mc3(q,beta,true,out);
        else mixture_minmod(q,out);
    };
    auto lim_minus=[&](const MixtureVec& q,MixtureVec& out) {
        if (grid.mc3_limiter) mixture_mc3(q,beta,false,out);
        else mixture_minmod(q,out);
    };
    lim_plus(r,lp_r); lim_minus(r,lm_r);
    lim_minus(r_ip1,lm_rip1); lim_plus(r_im1,lp_rim1);
    // First (predictor) pass. The predictor update below differences ONLY
    // fl_iph.flux − fr_imh.flux, so the i+½ RIGHT and i−½ LEFT reconstructions
    // (wr_iph, wl_imh) and their face bundles (fr_iph, fl_imh) are dead work here
    // — both are recomputed unconditionally by the corrector pass further down.
    // Building 6 face bundles per step instead of 8 is byte-identical.
    if (grid.uniform_mesh) {
        wl_iph=w+0.5*lp_r%(w_ip1-w);
        wr_imh=w-0.5*lm_r%(w_ip1-w);
    } else {
        wl_iph=w+W1%lp_r%(w_ip1-w);
        wr_imh=w-W1%lm_r%(w_ip1-w);
    }

    auto& p=grid.gamma_rhs_scratch.packed;
    MixtureFaceBundle fr_iph{p[0],p[1],p[2],p[3]};
    MixtureFaceBundle fl_iph{p[4],p[5],p[6],p[7]};
    MixtureFaceBundle fr_imh{p[8],p[9],p[10],p[11]};
    MixtureFaceBundle fl_imh{p[12],p[13],p[14],p[15]};
    build_mixture_face_batch2(grid,
        wl_iph, grid.phi_g_iph, fl_iph,
        wr_imh, grid.phi_g_imh, fr_imh);

    // Internal MUSCL predictor: update conservatively, then independently invert
    // the caloric EOS. This path intentionally never calls legacy cons2prim.
    Vec& predicted=p[16];
    predicted=xn_state-grid.dt_state/grid.ds_state%(fl_iph.flux-fr_imh.flux);
    DecodedMixtureField& decoded_predicted=grid.gamma_rhs_scratch.predicted;
    decode_predicted_caloric_primitives_into(
        grid,predicted,decoded,decoded_predicted);
    mixture_stencil_view(grid,decoded_predicted,0,wt);
    mixture_stencil_view(grid,decoded_predicted,1,wt_ip1);
    mixture_stencil_view(grid,decoded_predicted,-1,wt_im1);

    if (grid.uniform_mesh) {
        wr_iph=0.5*(w_ip1+wt_ip1)-0.5*lm_rip1%(w_ip2-w_ip1);
        wl_iph=0.5*(w+wt)+0.5*lp_r%(w_ip1-w);
        wr_imh=0.5*(w+wt)-0.5*lm_r%(w_ip1-w);
        wl_imh=0.5*(w_im1+wt_im1)+0.5*lp_rim1%(w-w_im1);
    } else {
        wr_iph=0.5*(w_ip1+wt_ip1)-W3%lm_rip1%(w_ip2-w_ip1);
        wl_iph=0.5*(w+wt)+W1%lp_r%(w_ip1-w);
        wr_imh=0.5*(w+wt)-W1%lm_r%(w_ip1-w);
        wl_imh=0.5*(w_im1+wt_im1)+W4%lp_rim1%(w-w_im1);
    }

    build_mixture_face_batch4(grid,
        wr_iph, grid.phi_g_iph, fr_iph,
        wl_iph, grid.phi_g_iph, fl_iph,
        wr_imh, grid.phi_g_imh, fr_imh,
        wl_imh, grid.phi_g_imh, fl_imh);
    Vec& a_imh=p[17]; Vec& a_iph=p[18];
    Vec& flux_imh=p[19]; Vec& flux_iph=p[20]; Vec& rhs=p[21];
    a_imh=arma::max(fl_imh.spectral_radius,fr_imh.spectral_radius);
    a_iph=arma::max(fl_iph.spectral_radius,fr_iph.spectral_radius);
    if (grid.roe_characteristic_flux) {
        build_roe_local_mixture_flux(grid,fl_imh,fr_imh,a_imh,flux_imh);
        build_roe_local_mixture_flux(grid,fl_iph,fr_iph,a_iph,flux_iph);
    } else {
        flux_imh=0.5f*(fl_imh.flux+fr_imh.flux
            -a_imh%(fr_imh.conserved-fl_imh.conserved));
        flux_iph=0.5f*(fl_iph.flux+fr_iph.flux
            -a_iph%(fr_iph.conserved-fl_iph.conserved));
    }
    Vec& source=p[22];
    mixture_source(grid,xn_state,decoded,source);
    rhs=-grid.B_state%(flux_iph/grid.B_state_iph
        -flux_imh/grid.B_state_imh)/grid.ds_state
        +source;
    if (grid.eq_wb && !grid.eq_residual.is_empty()) rhs -= grid.eq_residual;
    if (grid.capture_face_flux)
        capture_gamma_face_flux(grid,w,wl_iph,wr_iph,fl_iph,fr_iph,
                                a_iph,flux_iph,r,r_ip1,lp_r,lm_rip1);
    return Vec(rhs);
}

} // namespace

// ============================================================================
// Explicit RHS: TVD-MUSCL reconstruction + predictor-corrector +
// Rusanov flux differencing, minus the explicit source.
// (1.5D field-aligned; z-direction not considered. CoMFi Re_MUSCL analogue.)
// ============================================================================

Vec rhs_explicit_state(const Grid& grid, const Vec& xn_state) {
    if (!grid.eos_gamma_table.empty()) {
        const DecodedMixtureField decoded = decode_mixture_field(grid, xn_state);
        return rhs_explicit_mixture(grid, xn_state, decoded);
    }
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

Vec rhs_explicit_state(const Grid& grid, const Vec& xn_state,
                       const DecodedMixtureField& decoded) {
    if (grid.eos_gamma_table.empty()) return rhs_explicit_state(grid, xn_state);
    return rhs_explicit_mixture(grid, xn_state, decoded);
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
