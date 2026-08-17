/*!
 * @file single_fluid/integrator.cpp
 * @brief Release solver: implicit physical conduction and the release timestep.
 * @ingroup release_solver
 *
 *     U^n  --MUSCL-Hancock/Godunov hydro-->  U*  --implicit physical conduction-->  U^{n+1}
 *
 * Those are the only two stages, composed by first-order Lie splitting. There is
 * no equilibrium projection between them (U = (rho, rho u, E) is already the
 * authoritative state, and Saha equilibrium enters only through the EOS closure
 * that decodes it), and there are no volumetric energy stages: TRAC, nonthermal
 * beam heating, volumetric coronal heating, radiative cooling and
 * artificial/numerical conduction are NOT part of the release and live only in
 * the historical two-fluid research solver (src/two_fluid/).
 *
 * The conduction stage is a nonlinear backward-Euler solve on the mixture energy
 * row alone: Newton on the cell temperatures with a tridiagonal Jacobian, the
 * effective heat capacity and the conductivity taken from one fused Saha
 * evaluation per cell per pass. Mass and momentum are untouched, and the accepted
 * energy is the flux-updated internal energy plus the unchanged kinetic and
 * gravitational parts, so the stage is conservative by construction.
 */

#include "single_fluid/mixture.hpp"
#include "profiling.hpp"
#include "parallel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace chromosphere {

void MixtureConductionScratch::resize(std::size_t n) {
    rho.resize(n); e_old.resize(n); temperature.resize(n); target.resize(n);
    conductivity.resize(n); capacity.resize(n); n_e.resize(n); n_hi.resize(n);
    e_at_T.resize(n);
    g_left.resize(n); g_right.resize(n);
    a.resize(n); b.resize(n); c.resize(n); rhs.resize(n); delta.resize(n);
}

namespace {

inline arma::SizeMat state_size(const Grid& grid) {
    return arma::size(grid.ns, num_of_mixture_eq);
}

double mixture_kappa_e(double n_e, double n_hi, double temperature) {
    return 9.2048e-12*n_e*std::pow(temperature, 2.5)
         / (n_e + 2.836e-11*n_hi*temperature*temperature);
}

double mixture_kappa_n(double n_e, double n_hi, double temperature) {
    return 0.0342006*n_hi*temperature
         / (1.20613*n_e*std::sqrt(2.0*temperature)
            + 1.70573*n_hi*std::sqrt(temperature));
}

// The release conductivity is PHYSICAL only: Spitzer electron conduction plus
// neutral-hydrogen conduction. No TRAC broadening and no mesh-scaled artificial
// diffusivity enter the release operator anywhere.
double mixture_conductivity(double n_e, double n_hi, double temperature) {
    return mixture_kappa_e(n_e, n_hi, temperature)
         + mixture_kappa_n(n_e, n_hi, temperature);
}

void thomas_solve_double_inplace(
    const std::vector<double>& a, std::vector<double>& b,
    const std::vector<double>& c, std::vector<double>& d,
    std::vector<double>& x) {
    const std::size_t n = b.size();
    for (std::size_t i = 1; i < n; ++i) {
        const double m = a[i]/b[i-1];
        b[i] -= m*c[i-1];
        d[i] -= m*d[i-1];
    }
    x[n-1] = d[n-1]/b[n-1];
    for (std::size_t i = n-1; i > 0; --i)
        x[i-1] = (d[i-1]-c[i-1]*x[i])/b[i-1];
}

CaloricMixtureThermo decode_ghost(const Grid& grid, const Vec& ghost,
                                  double phi_face) {
    return decode_equilibrium_caloric_mixture(
        grid.eos_gamma_table, ghost(mix::RHO), ghost(mix::MOM),
        ghost(mix::ENERGY), phi_face, std::numeric_limits<double>::quiet_NaN(),
        grid.eos_gamma_debug_clamp);
}

// The outer state as the CONDUCTION rows see it. When the hydro and conduction
// temperatures are decoupled, the imposed temperature lives at the PHYSICAL
// boundary face, not at a fictitious ghost centre. Reconstruct the material
// state at that face from the externally imposed hydro back-pressure and the
// fixed face temperature, then evaluate kappa there. Without the override this
// remains the decoded hydro ghost state.
struct OuterConductionWall { double T, n_e, n_HI; };

OuterConductionWall outer_conduction_wall(const Grid& grid,
                                          const CaloricMixtureThermo& ghost) {
    if (!grid.outer_conduction_temperature_override)
        return {ghost.T, ghost.n_e, ghost.n_HI};
    const double T = static_cast<double>(grid.outer_conduction_temperature);
    const double rho = equilibrium_density_from_pressure(ghost.p, T);
    const double n_h = rho/eos_constants::m_h;
    const double x = saha_ionization_fraction_n_h(n_h, T);
    return {T, x*n_h, (1.0-x)*n_h};
}

} // namespace

// ============================================================================
// Energy replacement: the accepted conduction temperature written back into the
// energy row, keeping rho and rho u untouched.
// ============================================================================

namespace {

Vec mixture_set_internal_energy(
    const Grid& grid, const Vec& state,
    const std::vector<double>& target_internal_energy) {
    if (target_internal_energy.size() != grid.ns)
        throw std::invalid_argument("mixture internal-energy target has wrong size");
    Vec updated = state;
    const auto sz = state_size(grid);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const double rho = static_cast<double>(state(arma::sub2ind(sz, i, mix::RHO)));
        const double momentum =
            static_cast<double>(state(arma::sub2ind(sz, i, mix::MOM)));
        const double phi = mixture_cell_phi(grid, i);
        const double target = target_internal_energy[i];
        if (!std::isfinite(target))
            throw std::domain_error(
                "non-finite mixture internal-energy target at cell "
                + std::to_string(i) + ", time=" + std::to_string(grid.sim_time)
                + ", rho=" + std::to_string(rho));
        const double e_min = equilibrium_internal_energy(
            rho, grid.eos_gamma_table.min_temperature());
        const double e_max = equilibrium_internal_energy(
            rho, grid.eos_gamma_table.max_temperature());
        const double tolerance = 64.0*std::numeric_limits<double>::epsilon()
                               * std::max({std::abs(target), e_min, e_max});
        if (target < e_min-tolerance || target > e_max+tolerance)
            throw std::out_of_range(
                "mixture energy target outside the EOS caloric domain at cell "
                + std::to_string(i) + ", time=" + std::to_string(grid.sim_time)
                + ", rho=" + std::to_string(rho)
                + ", target_e=" + std::to_string(target)
                + ", valid_e=[" + std::to_string(e_min) + ","
                + std::to_string(e_max) + "]");
        // Round the packed float back into the caloric domain: independently
        // rounding E can otherwise place e_int = E - KE - rho phi a few ulp
        // outside the table edge and make the next decode throw.
        const double carried = 0.5*momentum*momentum/rho + rho*phi;
        float energy = static_cast<float>(target + carried);
        while (static_cast<double>(energy) - carried < e_min)
            energy = std::nextafter(energy, std::numeric_limits<float>::infinity());
        while (static_cast<double>(energy) - carried > e_max)
            energy = std::nextafter(energy, -std::numeric_limits<float>::infinity());
        if (!std::isfinite(energy))
            throw std::runtime_error(
                "failed to pack the mixture energy inside the EOS domain");
        updated(arma::sub2ind(sz, i, mix::ENERGY)) = energy;
    }
    return updated;
}

} // namespace

// ============================================================================
// Nonlinear backward-Euler physical conduction
// ============================================================================
//
// Solved for T with the analytic effective heat capacity C_V^eff = de_int/dT and
// only d(kappa)/dT lagged (a quasi-Newton scheme):
//
//   C_V^eff dT + dt [ g_L (dT_i - dT_{i-1}) - g_R (dT_{i+1} - dT_i) ]
//       = -( e_int(rho,T) - e_old - dt div q ),
//   g_{L,R} = (B_i / ds_i) (kappa_face / B_face) / ds_face.
//
// Mass and momentum never enter, so the stage conserves them exactly; the
// accepted energy is E = e_target + 1/2 rho u^2 + rho phi.
Vec mixture_apply_conduction(const Grid& grid, const Vec& state, double dt) {
    if (!(dt > 0.0)) return state;
    ProfileScope timer(ProfileRegion::Conduction);
    const arma::uword ns = grid.ns;
    MixtureConductionScratch& s = grid.mixture_conduction_scratch;
    s.resize(ns);
    const auto sz = state_size(grid);

    ParallelFailure failure;
    std::array<double, kMaximumParallelThreads> residual_max_slots{};
    std::array<double, kMaximumParallelThreads> step_max_slots{};
    int team_size = 1;
    bool stop = false;
    bool converged = false;
    bool caloric_seeded = true;
    bool previous_step_was_small = false;
    int converged_after_updates = 0;
    CaloricMixtureThermo inner{};
    CaloricMixtureThermo outer{};
    OuterConductionWall wall{};
    double k_inner = 0.0;
    double k_outer = 0.0;
    double diag_kr_top = 0.0;
    const int conduction_threads = std::min(
        kMaximumConductionThreads, parallel_max_threads());

    auto face_k = [&](double kh, double kt, double dh, double dtw) {
        if (grid.uniform_mesh) return 0.5*(kh+kt);
        kh = std::max(kh, 1.0e-30);
        kt = std::max(kt, 1.0e-30);
        return (dh+dtw)/(dh/kh+dtw/kt);
    };

    // One Grid owns one mutable conduction scratch set and is intentionally not
    // concurrently reentrant. All physical-cell writes below are index-disjoint.
#pragma omp parallel if(conduction_threads > 1) num_threads(conduction_threads) shared(stop,converged,caloric_seeded,previous_step_was_small,converged_after_updates,inner,outer,wall,k_inner,k_outer,diag_kr_top,team_size,failure,residual_max_slots,step_max_slots)
    {
        const int tid = parallel_thread_index();

#pragma omp for schedule(static)
        for (long long raw_i = 0; raw_i < static_cast<long long>(ns); ++raw_i) {
            const arma::uword i = static_cast<arma::uword>(raw_i);
            try {
                auto at = [&](arma::uword row) {
                    return static_cast<double>(state(arma::sub2ind(sz,i,row)));
                };
                const double phi = mixture_cell_phi(grid, i);
                CaloricState seeded;
                const CaloricMixtureThermo th = decode_equilibrium_caloric_mixture(
                    grid.eos_gamma_table, at(mix::RHO), at(mix::MOM),
                    at(mix::ENERGY), phi, grid.eos_temperature_hint(i),
                    grid.eos_gamma_debug_clamp, &seeded);
                grid.store_eos_temperature_hint(i, th.T);
                s.rho[i] = th.rho;
                s.e_old[i] = th.internal_energy;
                s.temperature[i] = th.T;
                s.n_e[i] = seeded.n_e;
                s.n_hi[i] = seeded.n_hi;
                s.e_at_T[i] = seeded.internal_energy;
                s.capacity[i] = seeded.heat_capacity;
            } catch (...) {
                failure.capture(i, std::current_exception());
            }
        }

#pragma omp single
        {
            team_size = parallel_team_size();
            if (failure.failed()) {
                stop = true;
            } else {
                try {
                    inner = decode_ghost(grid, grid.mix_inner_boundary0,
                                         grid.phi_g_imh(0));
                    outer = decode_ghost(grid, grid.mix_outer_boundary0,
                                         grid.phi_g_iph(ns-1));
                    wall = outer_conduction_wall(grid, outer);
                    k_inner = mixture_conductivity(inner.n_e, inner.n_HI, inner.T);
                    k_outer = mixture_conductivity(wall.n_e, wall.n_HI, wall.T);
                } catch (...) {
                    failure.capture(ns, std::current_exception());
                    stop = true;
                }
            }
        }

        if (!stop) {
            for (int iteration = 0; iteration <= 40; ++iteration) {
#pragma omp for schedule(static)
                for (long long raw_i = 0; raw_i < static_cast<long long>(ns); ++raw_i) {
                    const arma::uword i = static_cast<arma::uword>(raw_i);
                    try {
                        if (!caloric_seeded) {
                            const CaloricState cs = equilibrium_caloric_state(
                                s.rho[i], s.temperature[i]);
                            s.n_e[i] = cs.n_e;
                            s.n_hi[i] = cs.n_hi;
                            s.e_at_T[i] = cs.internal_energy;
                            s.capacity[i] = cs.heat_capacity;
                        }
                        s.conductivity[i] = mixture_conductivity(
                            s.n_e[i], s.n_hi[i], s.temperature[i]);
                    } catch (...) {
                        failure.capture(i, std::current_exception());
                    }
                }

#pragma omp single
                {
                    caloric_seeded = false;
                    if (failure.failed()) stop = true;
                }
                if (stop) break;

#pragma omp for schedule(static)
                for (long long raw_i = 0; raw_i < static_cast<long long>(ns); ++raw_i) {
                    const arma::uword i = static_cast<arma::uword>(raw_i);
                    try {
                        const double kl = i == 0
                            ? face_k(s.conductivity[i],k_inner,grid.ds_i(i),grid.ds_i(i))
                            : face_k(s.conductivity[i],s.conductivity[i-1],
                                     grid.ds_i(i),grid.ds_i(i-1));
                        // The imposed outer temperature is a PHYSICAL-FACE datum:
                        // the wall conductivity is used directly and the flux
                        // distance is the half top cell, not a ghost-centre span.
                        const bool physical_outer_face =
                            i+1 == ns && grid.outer_conduction_temperature_override;
                        const double kr = i+1 == ns
                            ? (physical_outer_face
                                ? k_outer
                                : face_k(s.conductivity[i],k_outer,grid.ds_i(i),grid.ds_i(i)))
                            : face_k(s.conductivity[i],s.conductivity[i+1],
                                     grid.ds_i(i),grid.ds_i(i+1));
                        const double ds_right = physical_outer_face
                            ? 0.5*static_cast<double>(grid.ds_i(i))
                            : static_cast<double>(grid.ds_iph_i(i));
                        if (i+1 == ns) diag_kr_top = kr;
                        s.g_left[i] = grid.B_i(i)/grid.ds_i(i)
                                  * kl/grid.B_imh(i)/grid.ds_imh_i(i);
                        s.g_right[i] = grid.B_i(i)/grid.ds_i(i)
                                   * kr/grid.B_iph(i)/ds_right;
                    } catch (...) {
                        failure.capture(i, std::current_exception());
                    }
                }

#pragma omp single
                {
                    if (failure.failed()) stop = true;
                }
                if (stop) break;

                double local_residual_max = 0.0;
#pragma omp for schedule(static)
                for (long long raw_i = 0; raw_i < static_cast<long long>(ns); ++raw_i) {
                    const arma::uword i = static_cast<arma::uword>(raw_i);
                    try {
                        const double tl = i==0 ? inner.T : s.temperature[i-1];
                        const double tr = i+1==ns ? wall.T : s.temperature[i+1];
                        double divergence = s.g_right[i]*(tr-s.temperature[i])
                                          -s.g_left[i]*(s.temperature[i]-tl);
                        if (i == 0 && grid.inner_conduction_neumann) {
                            divergence += s.g_left[i]*(s.temperature[i]-tl);
                            s.g_left[i] = 0.0;
                        }
                        if (i+1 == ns && grid.impose_outer_heat_flux) {
                            divergence = -s.g_left[i]*(s.temperature[i]-tl)
                                       + grid.outer_heat_flux/grid.ds_i(i);
                            s.g_right[i] = 0.0;
                        }
                        s.target[i] = s.e_old[i]+dt*divergence;
                        const double residual = s.e_at_T[i]-s.target[i];
                        s.a[i] = -dt*s.g_left[i];
                        s.b[i] = s.capacity[i]+dt*(s.g_left[i]+s.g_right[i]);
                        s.c[i] = -dt*s.g_right[i];
                        s.rhs[i] = -residual;
                        local_residual_max = std::max(local_residual_max,
                            std::abs(residual)/std::max(s.e_old[i],1.0e-30));
                    } catch (...) {
                        failure.capture(i, std::current_exception());
                    }
                }
                residual_max_slots[static_cast<std::size_t>(tid)] = local_residual_max;
#pragma omp barrier

#pragma omp single
                {
                    if (failure.failed()) {
                        stop = true;
                    } else {
                        s.a[0]=0.0;
                        s.c[ns-1]=0.0;
                        double max_scaled_residual = 0.0;
                        for (int slot = 0; slot < team_size; ++slot)
                            max_scaled_residual = std::max(
                                max_scaled_residual,
                                residual_max_slots[static_cast<std::size_t>(slot)]);
                        if (max_scaled_residual < 2.0e-11) {
                            converged = true;
                            converged_after_updates = iteration;
                            stop = true;
                        } else if (previous_step_was_small) {
                            failure.capture(ns, std::make_exception_ptr(std::runtime_error(
                                "mixture conduction stagnated before residual convergence: "
                                "scaled_residual="+std::to_string(max_scaled_residual))));
                            stop = true;
                        } else if (iteration == 40) {
                            stop = true;
                        } else {
                            try {
                                thomas_solve_double_inplace(s.a,s.b,s.c,s.rhs,s.delta);
                            } catch (...) {
                                failure.capture(ns, std::current_exception());
                                stop = true;
                            }
                        }
                    }
                }
                if (stop) break;

                double local_step_max = 0.0;
#pragma omp for schedule(static)
                for (long long raw_i = 0; raw_i < static_cast<long long>(ns); ++raw_i) {
                    const arma::uword i = static_cast<arma::uword>(raw_i);
                    try {
                        double candidate = s.temperature[i]+s.delta[i];
                        if (!std::isfinite(candidate))
                            throw std::runtime_error(
                                "mixture conduction Newton produced a non-finite T");
                        const double tmin = grid.eos_gamma_table.min_temperature();
                        const double tmax = grid.eos_gamma_table.max_temperature();
                        if (candidate<=tmin)
                            candidate=0.5*(s.temperature[i]+tmin);
                        if (candidate>=tmax)
                            candidate=0.5*(s.temperature[i]+tmax);
                        local_step_max = std::max(local_step_max,
                            std::abs(candidate-s.temperature[i])/
                            std::max(s.temperature[i],1.0));
                        s.temperature[i]=candidate;
                    } catch (...) {
                        failure.capture(i, std::current_exception());
                    }
                }
                step_max_slots[static_cast<std::size_t>(tid)] = local_step_max;
#pragma omp barrier

#pragma omp single
                {
                    if (failure.failed()) {
                        stop = true;
                    } else {
                        double max_relative_step = 0.0;
                        for (int slot = 0; slot < team_size; ++slot)
                            max_relative_step = std::max(
                                max_relative_step,
                                step_max_slots[static_cast<std::size_t>(slot)]);
                        previous_step_was_small = max_relative_step < 2.0e-11;
                    }
                }
                if (stop) break;
            }
        }
    }

    failure.rethrow_lowest();
    if (!converged)
        throw std::runtime_error("mixture conduction Newton did not converge");
    profile_note_conduction_iterations(
        static_cast<std::uint64_t>(converged_after_updates));

    if (grid.capture_outer_conduction) {
        OuterConductionCapture& oc = grid.outer_conduction_capture;
        oc.T_top  = s.temperature[ns-1];
        oc.T_wall = wall.T;
        oc.kappa_face = diag_kr_top;
        oc.ds_face    = grid.outer_conduction_temperature_override
            ? 0.5*static_cast<double>(grid.ds_i(ns-1))
            : static_cast<double>(grid.ds_iph_i(ns-1));
        oc.area_ratio = grid.B_i(ns-1)/grid.B_iph(ns-1);
        oc.q_face = oc.kappa_face*(oc.T_wall-oc.T_top)/oc.ds_face;
        oc.imposed_neumann = grid.impose_outer_heat_flux;
        oc.valid = true;
    }

    return mixture_set_internal_energy(grid, state, s.target);
}

double mixture_conduction_residual_max(const Grid& grid, const Vec& before,
                                       const Vec& after, double dt) {
    if (grid.eos_gamma_table.empty())
        throw std::logic_error("mixture conduction residual requires a Gamma1 table");
    const arma::uword ns = grid.ns;
    const auto sz = state_size(grid);
    std::vector<double> conductivity(ns), e_old(ns), e_now(ns), temperature(ns);
    for (arma::uword i = 0; i < ns; ++i) {
        const double phi = mixture_cell_phi(grid, i);
        auto decode = [&](const Vec& st) {
            return decode_equilibrium_caloric_mixture(
                grid.eos_gamma_table,
                static_cast<double>(st(arma::sub2ind(sz, i, mix::RHO))),
                static_cast<double>(st(arma::sub2ind(sz, i, mix::MOM))),
                static_cast<double>(st(arma::sub2ind(sz, i, mix::ENERGY))), phi);
        };
        const CaloricMixtureThermo old_th = decode(before);
        const CaloricMixtureThermo new_th = decode(after);
        e_old[i] = old_th.internal_energy;
        e_now[i] = new_th.internal_energy;
        temperature[i] = new_th.T;
        conductivity[i] = mixture_conductivity(new_th.n_e, new_th.n_HI,
                                               new_th.T);
    }
    const CaloricMixtureThermo inner = decode_ghost(
        grid, grid.mix_inner_boundary0, grid.phi_g_imh(0));
    const CaloricMixtureThermo outer = decode_ghost(
        grid, grid.mix_outer_boundary0, grid.phi_g_iph(ns-1));
    const double k_inner = mixture_conductivity(inner.n_e, inner.n_HI, inner.T);
    const OuterConductionWall wall = outer_conduction_wall(grid, outer);
    const double k_outer = mixture_conductivity(wall.n_e, wall.n_HI, wall.T);
    auto face_k = [&](double kh, double kt, double dh, double dtw) {
        if (grid.uniform_mesh) return 0.5*(kh+kt);
        kh=std::max(kh,1.0e-30); kt=std::max(kt,1.0e-30);
        return (dh+dtw)/(dh/kh+dtw/kt);
    };
    double maximum = 0.0;
    for (arma::uword i = 0; i < ns; ++i) {
        const double kl = i == 0
            ? face_k(conductivity[i],k_inner,grid.ds_i(i),grid.ds_i(i))
            : face_k(conductivity[i],conductivity[i-1],grid.ds_i(i),grid.ds_i(i-1));
        const bool physical_outer_face =
            i+1 == ns && grid.outer_conduction_temperature_override;
        const double kr = i+1 == ns
            ? (physical_outer_face
                ? k_outer
                : face_k(conductivity[i],k_outer,grid.ds_i(i),grid.ds_i(i)))
            : face_k(conductivity[i],conductivity[i+1],grid.ds_i(i),grid.ds_i(i+1));
        double gl = grid.B_i(i)/grid.ds_i(i)*kl/grid.B_imh(i)/grid.ds_imh_i(i);
        const double ds_right = physical_outer_face
            ? 0.5*static_cast<double>(grid.ds_i(i))
            : static_cast<double>(grid.ds_iph_i(i));
        double gr = grid.B_i(i)/grid.ds_i(i)*kr/grid.B_iph(i)/ds_right;
        const double tl = i ? temperature[i-1] : inner.T;
        const double tr = i+1 < ns ? temperature[i+1] : wall.T;
        double divergence = gr*(tr-temperature[i])-gl*(temperature[i]-tl);
        if (i == 0 && grid.inner_conduction_neumann)
            divergence += gl*(temperature[i]-tl);
        if (i+1 == ns && grid.impose_outer_heat_flux)
            divergence = -gl*(temperature[i]-tl)+grid.outer_heat_flux/grid.ds_i(i);
        const double residual = e_now[i] - e_old[i] - dt*divergence;
        maximum = std::max(maximum, std::abs(residual)/std::max(e_old[i],1.0e-30));
    }
    return maximum;
}

// ============================================================================
// Release timestep
// ============================================================================

Vec mixture_advance(Grid& grid, const Vec& state, const Vec& dt_i,
                    const MixtureField& decoded) {
    if (grid.eos_gamma_table.empty())
        throw std::logic_error("the mixture solver requires a Gamma1 table");
    decoded.require_matches(grid, state);
    const double dt = dt_i(0);   // uniform by mixture_timestep construction

    Vec next = state;
    {
        const Vec rhs = mixture_rhs_explicit(grid, state, decoded, dt);
        const auto sz = state_size(grid);
        for (arma::uword i = 0; i < grid.ns; ++i)
            for (arma::uword row = 0; row < num_of_mixture_eq; ++row) {
                const arma::uword k = arma::sub2ind(sz, i, row);
                next(k) = state(k) + dt_i(i)*rhs(k);
            }
    }
    if (grid.enable_conduction)
        next = mixture_apply_conduction(grid, next, dt);
    return next;
}

} // namespace chromosphere
