/*!
 * @file single_fluid/exact_rs.cpp
 * @brief Port of SWMF `share/Library/src/ModExactRS.f90`.
 * @ingroup release_solver
 *
 * See single_fluid/exact_rs.hpp. Every branch below mirrors a branch of the
 * Fortran; the comments name the Fortran routine each block came from so the two
 * can be diffed by eye.
 */

#include "single_fluid/exact_rs.hpp"

#include <algorithm>
#include <cmath>

namespace chromosphere {
namespace swmf_rs {

namespace {

inline bool admissible(const ExactRsState& s) {
    return std::isfinite(s.rho) && s.rho > 0.0
        && std::isfinite(s.p)   && s.p   > 0.0
        && std::isfinite(s.u);
}

/// SWMF `pressure_function`: the Toro pressure function f_K(p) of one side, its
/// derivative, and that side's total perturbation speed relative to the state.
/// `p` [Pa] is the current star-pressure iterate; `rho_k`, `p_k`, `c_k` are the
/// side's density, pressure and sound speed.
void pressure_function(double p, double rho_k, double p_k, double c_k,
                       double gamma, double& f, double& fd, double& wk) {
    if (p <= p_k) {
        // Rarefaction.
        const double p_ratio = p/p_k;
        const double p_ratio_pow = std::pow(p_ratio, (gamma - 1.0)/(2.0*gamma));
        f  = c_k*(p_ratio_pow - 1.0)*2.0/(gamma - 1.0);
        fd = p_ratio_pow/(p_ratio*rho_k*c_k);
        wk = c_k;
    } else {
        // Shock.
        const double bk = p_k*(gamma - 1.0)/(gamma + 1.0) + p;
        const double qrt = std::sqrt(2.0/((gamma + 1.0)*rho_k*bk));
        f  = (p - p_k)*qrt;
        fd = (1.0 - 0.5*(p - p_k)/bk)*qrt;
        wk = 0.50*(gamma + 1.0)*bk*qrt;
    }
}

/// SWMF `guess_p` (P. Voinovich, 1992): the acoustic/linearised star pressure,
/// floored at a tiny fraction of the impedance-weighted mean so the first Newton
/// evaluation is always at a positive pressure.
double guess_p(const ExactRsState& left, const ExactRsState& right,
               double c_left, double c_right) {
    constexpr double tiny_fraction_of = 1.0e-8;
    const double impedance_l = left.rho*c_left;
    const double impedance_r = right.rho*c_right;
    const double impedance_total_inv = 1.0/(impedance_l + impedance_r);
    const double p_avr =
        impedance_total_inv*(left.p*impedance_r + right.p*impedance_l);
    return std::max(tiny_fraction_of*p_avr,
                    p_avr + (left.u - right.u)*impedance_total_inv
                            *impedance_l*impedance_r);
}

/// SWMF `simple_wave`, the routine nested inside `exact_rs_sample`. It works in
/// the frame in which the contact is at rest and the sampled side is on the
/// right, so the caller mirrors the left state into it and un-mirrors the result.
/// `s_rel` [m s^-1] is the similarity coordinate in that frame.
void simple_wave(double s_rel, double rho_k, double un_k, double p_k, double c_k,
                 double gamma, double p_star,
                 double& rho, double& un, double& p) {
    if (p_star <= p_k) {
        // Rarefaction.
        if (s_rel < c_k*std::pow(p_star/p_k, (gamma - 1.0)/(2.0*gamma))) {
            // Between the contact and the rarefaction tail: the star state.
            un = 0.0;
            p  = p_star;
        } else {
            // Inside the fan. Solve
            //   c_k - g7 un_k = c - g7 un   (the Riemann invariant)
            //   s_rel = c + un              (the characteristic through the point)
            const double g7 = (gamma - 1.0)/2.0;
            const double c = (c_k - g7*(un_k - s_rel))/(1.0 + g7);
            un = s_rel - c;
            p  = p_k*std::pow(c/c_k, 2.0*gamma/(gamma - 1.0));
        }
        rho = rho_k*std::pow(p/p_k, 1.0/gamma);
    } else {
        // Shock: the sampled point is the star state behind it.
        rho = rho_k*((gamma + 1.0)*p_star + p_k*(gamma - 1.0))
                   /((gamma - 1.0)*p_star + p_k*(gamma + 1.0));
        un  = 0.0;
        p   = p_star;
    }
}

} // namespace

ExactRsSolution exact_rs_pu_star(const ExactRsState& left,
                                 const ExactRsState& right,
                                 double gamma_left, double gamma_right) {
    ExactRsSolution out{};
    out.status = ExactRsStatus::Ok;
    out.p_star = 0.0;
    out.u_star = 0.0;
    out.wl = 0.0;
    out.wr = 0.0;
    out.iterations = 0;
    if (!admissible(left) || !admissible(right)
        || !(gamma_left > 1.0) || !(gamma_right > 1.0)) {
        out.status = ExactRsStatus::BadInput;
        return out;
    }

    const double c_left  = std::sqrt(gamma_left*left.p/left.rho);
    const double c_right = std::sqrt(gamma_right*right.p/right.rho);
    out.c_left = c_left;
    out.c_right = c_right;

    // Safe signal-speed bounds, used verbatim by SWMF whenever it bails out.
    auto bail = [&](ExactRsStatus status) {
        out.status = status;
        out.wl = std::min(std::min(0.0, left.u - c_left), right.u - c_right);
        out.wr = std::max(std::max(0.0, left.u + c_left), right.u + c_right);
        return out;
    };

    // Speed of expansion into vacuum. If the states separate faster than this a
    // vacuum opens and there is no star region.
    const double u_expansion_l = 2.0*c_left/(gamma_left - 1.0);
    const double u_expansion_r = 2.0*c_right/(gamma_right - 1.0);
    const double u_vacuum = u_expansion_l + u_expansion_r;
    const double u_diff = right.u - left.u;
    if (u_diff > kVacuumSafetyFactor*u_vacuum) return bail(ExactRsStatus::Vacuum);

    double p_old = guess_p(left, right, c_left, c_right);
    double f_l = 0.0, fd_l = 0.0, w_l = 0.0;
    double f_r = 0.0, fd_r = 0.0, w_r = 0.0;
    double change = 2.0*kTolerance;
    double p_star = p_old;
    int iteration = 0;
    while (change > kTolerance && iteration < kMaxIterations) {
        pressure_function(p_old, left.rho,  left.p,  c_left,  gamma_left,
                          f_l, fd_l, w_l);
        pressure_function(p_old, right.rho, right.p, c_right, gamma_right,
                          f_r, fd_r, w_r);
        p_star = p_old - (f_l + f_r + u_diff)/(fd_l + fd_r);
        if (!(p_star > 0.0)) return bail(ExactRsStatus::NegativePressure);
        change = 2.0*std::abs((p_star - p_old)/(p_star + p_old));
        p_old = p_star;
        ++iteration;
    }
    out.iterations = iteration;
    if (!std::isfinite(p_star)) return bail(ExactRsStatus::NegativePressure);
    if (change > kTolerance) {
        // Report the budget overrun but keep the last iterate: SWMF accepts it,
        // and the caller decides whether an unconverged star state is usable.
        out.status = ExactRsStatus::NotConverged;
    }

    out.p_star = p_star;
    out.u_star = 0.50*(left.u + right.u + f_r - f_l);
    // SWMF's end-of-routine fixup: turn the per-side perturbation speeds into
    // absolute leftmost/rightmost signal speeds.
    out.wl = left.u  - w_l;
    out.wr = right.u + w_r;
    if (!std::isfinite(out.u_star) || !std::isfinite(out.wl)
        || !std::isfinite(out.wr))
        return bail(ExactRsStatus::NegativePressure);
    return out;
}

ExactRsState exact_rs_sample(double s, const ExactRsSolution& solution,
                             const ExactRsState& left, const ExactRsState& right,
                             double gamma_left, double gamma_right) {
    ExactRsState out{};
    if (s > solution.wr) return right;   // outside the rightmost wave
    if (s < solution.wl) return left;    // outside the leftmost wave
    if (s <= solution.u_star) {
        // Left of the contact. Mirror into the simple-wave frame and back.
        simple_wave(solution.u_star - s, left.rho, solution.u_star - left.u,
                    left.p, solution.c_left, gamma_left, solution.p_star,
                    out.rho, out.u, out.p);
        out.u = solution.u_star - out.u;
    } else {
        simple_wave(s - solution.u_star, right.rho, right.u - solution.u_star,
                    right.p, solution.c_right, gamma_right, solution.p_star,
                    out.rho, out.u, out.p);
        out.u = out.u + solution.u_star;
    }
    return out;
}

} // namespace swmf_rs
} // namespace chromosphere
