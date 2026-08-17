/*!
 * @file single_fluid/exact_rs.hpp
 * @brief A direct port of the SWMF exact ideal-gas Riemann solver
 *        (`share/Library/src/ModExactRS.f90`), used by the RELEASE `swmf_godunov`
 *        numerical flux.
 * @ingroup release_solver
 *
 * This is a line-for-line port of the two routines the SWMF Godunov scheme calls,
 * `exact_rs_pu_star` and `exact_rs_sample`, together with their two helpers
 * `pressure_function` and `guess_p`. The algorithm is Toro's exact two-rarefaction/
 * two-shock pressure iteration (Toro, *Riemann Solvers and Numerical Methods for
 * Fluid Dynamics*, 2nd ed., ch. 4), in the SWMF revision that allows
 * `gamma_left != gamma_right` and replaces Toro's guess with P. Voinovich's
 * impedance-weighted one.
 *
 * It solves the IDEAL-GAS Riemann problem only: `p = (gamma - 1) e_int`. The
 * chromosphere release EOS is not an ideal gas, and this header makes no attempt
 * to hide that — the non-ideal content is carried outside the solve, by the
 * energy-offset treatment of `single_fluid/mixture.cpp`, exactly as SWMF's
 * `util/CRASH/src/test_godunov.f90` does for the CRASH EOS.
 *
 * This solver IS on the release path: `single_fluid/mixture.cpp` calls it at every
 * face, at the frozen-composition index `gamma = 5/3`, with the ionization energy
 * carried in the passive offset. The mixture Roe solver it replaced is retained as
 * the equilibrium-linearization reference flux (`ISO_RIEMANN=roe-local`).
 *
 * ## Deviations from the Fortran
 *
 * The Fortran module keeps `PStar`, `UnStar`, `WL`, `WR`, `CL`, `CR` and the six
 * input primitives in OpenMP-threadprivate module variables, and `sample` reads
 * them by host association. This port returns them in an ExactRsSolution instead,
 * so the solver is reentrant and every face is independent. The arithmetic,
 * branch structure, iteration budget and tolerances are unchanged.
 *
 * Where the Fortran calls `CON_stop` on a vacuum state, or silently retries with
 * `PStar = 0.1 pOld` when the caller passed no `UseAnotherRS` flag, this port
 * always behaves as the SWMF caller that DOES pass `UseAnotherRS`: it reports the
 * failure and leaves the fallback decision to the caller.
 */
#pragma once

namespace chromosphere {
namespace swmf_rs {

/** @addtogroup release_solver
 *  @{
 */

/// Outcome of one exact Riemann solve. Anything other than `Ok` means the
/// returned star state must not be used and the caller has to fall back.
enum class ExactRsStatus {
    Ok,               ///< The iteration converged to a positive star pressure.
    BadInput,         ///< A side state had non-positive or non-finite rho, p or u.
    Vacuum,           ///< u_R - u_L exceeded the two-rarefaction escape speed: a
                      ///< vacuum would open between the waves and there is no
                      ///< star region to solve for.
    NegativePressure, ///< A Newton step produced p* <= 0.
    NotConverged      ///< The 10-iteration budget ran out before the 1e-3
                      ///< relative-change tolerance was met.
};

/// One side of a Riemann problem, in primitive variables.
struct ExactRsState {
    double rho;  ///< mass density [kg m^-3]
    double u;    ///< face-normal velocity [m s^-1]
    double p;    ///< pressure [Pa]
};

/// Star-region solution plus everything `exact_rs_sample` needs to evaluate the
/// wave fan. Only valid when `status == ExactRsStatus::Ok`, except for `wl`/`wr`,
/// which are also filled with the safe signal-speed bounds on a `Vacuum` or
/// `NegativePressure` failure so a caller can still size its fallback dissipation.
struct ExactRsSolution {
    ExactRsStatus status;  ///< Why the solve stopped; see ExactRsStatus.
    double p_star;         ///< pressure in the star region [Pa]
    double u_star;         ///< velocity of the contact discontinuity [m s^-1]
    /// Speed of the leftmost signal [m s^-1]: `u_L - c_L` for a left rarefaction
    /// (its head), or the left shock speed. This is SWMF's `WL` after the
    /// `WL = UnL - WL` fixup at the end of `exact_rs_pu_star`.
    double wl;
    /// Speed of the rightmost signal [m s^-1], the mirror image of `wl`.
    double wr;
    double c_left;         ///< left-state sound speed sqrt(gamma p / rho) [m s^-1]
    double c_right;        ///< right-state sound speed sqrt(gamma p / rho) [m s^-1]
    int iterations;        ///< Newton iterations actually taken [-] (1..10).
};

/// Relative-change tolerance of the star-pressure iteration [-]. SWMF's `TolP`.
constexpr double kTolerance = 0.0010;
/// Hard cap on Newton iterations [-]. SWMF's `nIterMax`.
constexpr int kMaxIterations = 10;
/// Fraction of the two-rarefaction escape speed above which the solve is declared
/// a vacuum [-]. SWMF's `cSafetyFactor`.
constexpr double kVacuumSafetyFactor = 0.999;

/// Solve for the pressure and velocity in the star region of the ideal-gas
/// Riemann problem `(left | right)`, allowing a different ratio of specific heats
/// on each side. Port of SWMF `exact_rs_pu_star` with `UseAnotherRS` present.
///
/// @param left     left state (rho [kg m^-3], u [m s^-1], p [Pa])
/// @param right    right state, same units
/// @param gamma_left  ratio of specific heats of the left state [-]
/// @param gamma_right ratio of specific heats of the right state [-]
/// @return the star state and wave speeds; check `status` before using them.
ExactRsSolution exact_rs_pu_star(const ExactRsState& left,
                                 const ExactRsState& right,
                                 double gamma_left, double gamma_right);

/// Sample the self-similar Riemann solution along the ray `x/t = s`. Port of SWMF
/// `exact_rs_sample`, including its nested `simple_wave`. The Godunov flux needs
/// `s = 0`, but the parameter is kept general so the shock-tube regression test can
/// reconstruct the whole wave pattern and compare it against an independent
/// solution.
///
/// @param s        similarity coordinate x/t [m s^-1]
/// @param solution a converged `exact_rs_pu_star` result
/// @param left     the same left state that solution was computed from
/// @param right    the same right state
/// @param gamma_left  ratio of specific heats of the left state [-]
/// @param gamma_right ratio of specific heats of the right state [-]
/// @return the sampled primitive state (rho [kg m^-3], u [m s^-1], p [Pa])
ExactRsState exact_rs_sample(double s, const ExactRsSolution& solution,
                             const ExactRsState& left, const ExactRsState& right,
                             double gamma_left, double gamma_right);

/** @} */

} // namespace swmf_rs
} // namespace chromosphere
