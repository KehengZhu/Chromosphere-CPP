/*!
 * Chromosphere model header file.
 * This header contains code from CoMFi which has an MIT License
 */

#pragma once

#include <armadillo> // Vector library
#include <cmath>

using namespace std;
using namespace arma;

/// Armadillo's fortran style column vector
typedef arma::Col<float> Vec;

/*!
 * Chromosphere model namespace.
 */
namespace chromosphere {
/* By Keheng: the data structure xn is shown in
    II.2.8 Data Structure of Qusai's PhD dissertation.
    The order is a little different, though.
    https://github.com/qalshidi/phd-dissertation
*/
// SOLUTION INDEX
const uword num_of_eq = 6; // Number of equations
// Conserved Variable Indices
const uword CNI = 0;        // Ion density
const uword CNN = 1;        // Neutral density
const uword CNV = 2;        // Ion velocity
const uword CNU = 3;        // Neutral velocity
const uword CEI = 4;        // Ion energy
const uword CEN = 5;        // Neutral energy
// Primitive Variable Indices
const uword PNI = 0;        // Ion density
const uword PNN = 1;        // Neutral density
const uword PV  = 2;        // Ion velocity
const uword PU  = 3;        // Neutral velosity
const uword PPI = 4;        // Ion pressure
const uword PPN = 5;        // Neutral pressure

// Finite difference methods
const uword FORWARD = 0;
const uword BACKWARD = 1;
const uword CENTRAL = 2;

// for the permutation functions
const uword SLICE = 1;
const uword CUBE = num_of_eq;

// Small enough
const float SMALL_ENOUGH = 1.0e-6;

// These are all global variables and may cause problems when attempts are made to parallize the code.
// If we want to parallelize the code with OpenMP, we better encapsulate these variable into a structure.

extern uword ns;       // Number of grid points in the x-direction
extern uword num_of_elem; // Total number of elements (ns * num_of_eq)
// constants
extern float gammamono; // Ratio of specific heats (adiabatic index for a monatomic ideal gas)
extern float alpha_p;   // Divergence error propagation parameter (alpha_p)
extern float pi;        // The value of π (pi)

// Physical constants
extern float m_i;       // Mass of proton in kilograms (also m0, normalization constant for mass)
extern float m_n;       //  mass of neutrals
extern float m_e;       //  mass of electron
extern float g;         // Gravitational acceleration in normalized units, default is Sun value of 0.27395 km/s^2
extern float mu_0;      // Permeability of free space (magnetic constant)
extern float k_b;       // Boltzmann constant in SI units
extern float e_;        // Elementary charge in Coulombs

// Domain normalization constants
extern Vec   ds_i;      // The length of i-th cell. Index range: 0~N-1
extern float CFL;       // Courant-Friedrichs-Lewy condition

// stream-aligned quantities. length: ns
extern Vec cellX_i, cellY_i, cellZ_i; // Cell coordinates XYZ.
extern Vec B_imh;                  // B at the left face of i-th cell.
extern Vec B_iph;                  // B at the right face of i-the cell.
extern Vec B_i;                    // B at the i-th cell.
extern Vec dinvB_ds_i;             // \frac{\partial}{\partial s}(\frac{1}{B}) at i-th cell.
extern Vec gPotential_imh;
extern Vec gPotential_iph;
extern Vec outer_boundary0_i;
extern Vec outer_boundary1_i;
extern Vec inner_boundary0_i;
extern Vec inner_boundary1_i;

const bool USE_NEUMANN_BC = false;

// Broadcasted variables to avoid redundant computations. Length: ns*num_of_eq
extern Vec B_iimh, B_iiph, B_ii, dinvB_ds_ii, dt_ii, ds_ii;

// Initial conditions
Vec shock_tube_ic();

/*!
 * By Keheng:
 * Advance methods: Euler and RK4
 */

// Get local time steps
Vec cal_dt_i(const Vec& xn);
Vec advance_Euler_ii(const Vec& xn, const Vec& dt);
Vec advance_RK4(const Vec& xn, const Vec& dt);

// refer to (10) in the paper.
/* ----------------------
arma operators:
%	 element-wise multiplication of two objects (Schur product)
------------------------- */
inline float nu_nn(const float& nn, const float& T) {
    const float sigma_nn = 7.73e-19; // m-2
    const float m_nn = 1.007825 * arma::datum::m_u;
    const float nn_coeff =
        sigma_nn * std::sqrt(16.0 * arma::datum::k / (arma::datum::pi * m_nn));
    return (nn_coeff * nn * std::sqrt(std::abs(T))); // ion-neutral collision rate
}

inline Vec nu_nn(const Vec& nn, const Vec& T) {
    const float sigma_nn = 7.73e-19; // m-2
    const float m_nn = 1.007825 * arma::datum::m_u;
    const float nn_coeff =
        sigma_nn * std::sqrt(16.0 * arma::datum::k / (arma::datum::pi * m_nn));
    return (nn_coeff * nn % arma::sqrt(arma::abs(T))); // ion-neutral collision rate
}

inline float nu_in(const float& nn, const float& T) {
    const float sigma_in = 1.16e-18; // m-2
    const float m_in = 1.007276466879 * 1.007825 / (1.007276466879 + 1.007825) * arma::datum::m_u;
    const float in_coeff =
        sigma_in * std::sqrt(8.0 * arma::datum::k / (arma::datum::pi * m_in));
    return (in_coeff * nn * std::sqrt(std::abs(T))); // ion-neutral collision rate
}

// These quantities are confusing: 1.007276466879 and 1.007825
// They are not the mass of proton and neutron
// the unit is normalized unit because it is divided by t_0^{-1}
inline Vec nu_in(const Vec& nn, const Vec& Ti, const Vec& Tn) {
    const float bohr_r = 53e-12;
    return (2.0*bohr_r)*(2.0*bohr_r)*nn%arma::sqrt(arma::abs(8.0*pi*k_b*(Ti+Tn)/m_i));
}


inline Vec nu_en(const Vec& nn, const Vec& T) {
    // Return electron-neutral collision frequency in SI units.
    return 1.95e-16*nn%arma::sqrt(arma::abs(T));
}

/*
\nu_{e i}(v)=\frac{n_i e_i^2 e^2}{4 \pi m_e^2 \epsilon_0^2 v^3} \ln \Lambda
*/
inline Vec nu_ei(const Vec& ni, const Vec& Te) {
    // first calculate by SI units and then normalize
    const float qe = 1.60217663e-19; // C
    const float eps_0 = 8.854187817e-12; // F/m
    const float me = 9.10938356e-31; // kg
    const Vec v_th = arma::sqrt(arma::abs(arma::datum::k * Te / me)); // m/s
    const float coef = ((double)qe*qe*qe*qe)/(4.0*pi*(double)me*me*eps_0*eps_0)*20.0;
    return coef*ni/(v_th%v_th%v_th);
}

inline Vec kappa_e(const Vec& ne, const Vec& nn, const Vec& Te) {
    // Calculate electron heat conductivity in SI units.
    return (9.2048e-12*ne%arma::pow(Te, 2.5)/(ne + 3.5609e-12*nn%Te%Te));
}

inline Vec kappa_n(const Vec& ni, const Vec& nn, const Vec& Ti, const Vec& Tn) {
    // Calculate neutral heat conductivity in SI units.
    return (0.0342006*nn%Tn)/(1.20613*ni%arma::sqrt(Tn+Ti) + 1.70573*nn%arma::sqrt(Tn));
}

/*!
 * Permutation horizontally right 1: i+1
 */
Vec ip1(const Vec& xn, const uword nk = CUBE);

/*!
 * Permutation horizontally left 1: i-1
 */
Vec im1(const Vec& xn, const uword nk = CUBE);

Vec ip2(const Vec& xn, const uword nk = CUBE);
Vec im2(const Vec& xn, const uword nk = CUBE);

/*!
 * Derivative in the x-direction
 * type = FORWARD (0): forward difference
 * type = BACKWARD (1): backward difference
 * type = CENTRAL (2): central difference
 */
Vec deriv_s(const Vec& xn, const uword type = FORWARD);

/* Norm of a vector*/
Vec vec_abs(const Vec& Ax, const Vec& Ay, const Vec& Az);
Vec vec_abs2(const Vec& Ax, const Vec& Ay, const Vec& Az);

Vec cons2prim(const Vec& cons);
Vec prim2cons(const Vec& prim);

/*!
 * Get scalar values for vector of size ns
 */
Vec get_scalar(const Vec& xn_ii, const uword& index);

/*!
 * Put scalars to unknown of larger vector
 */
Vec scalar_to(const Vec& xn_i, const uword& index);

/*!
 * Flux limiter. Currently using ospre
 */
Vec flux_lim(const Vec& r);

/*!
 * Calculate explicit terms in right hand equation.
 */
Vec rhs_explicit_ii(const Vec& xn_ii, const float& dt_i);

/*!
 * Calculate implicit terms in right hand equation.
 */
Vec rhs_implicit_ii(const Vec& xn_ii);

/*!
 * Calculate the fast speed eigenvalue in the x-direction
 */
Vec find_spectral_radius_ii(const Vec& xn_ii);

/*!
 * Calculate flux in the x-direction
 */
Vec cal_F_ii(const Vec& xn_ii);

/*!
 * Get max speed for courant condition calculation
 */
Vec get_max_v_i(const Vec& xn_ii);

void print_xn(const Vec& xn);

} // namespace chromosphere

//--------------------------------------------------------------------
// Standalone driver: Model C7 chromosphere setup
namespace chromosphere {

/*!
 * Initialize physical constants, allocate stream-aligned arrays.
 * Sets ns, CFL, gammamono, masses, k_b, e_, mu_0, pi, g, gPotential, B.
 */
void chromo_init(uword ns_in, float CFL_in);

/*!
 * Build Model C7 initial condition (Table 1 of writeup),
 * populating the conserved-variable state vector and boundary buffers.
 * Returns the initial xn (length ns*num_of_eq).
 */
Vec model_c7_ic();

/*!
 * Refresh outer boundary at each step: ion/neutral densities and temperatures
 * are pinned to Model C7; outer velocity is halved each call (outflow damping).
 */
void model_c7_update_bc(const Vec& xn);

} // namespace chromosphere
