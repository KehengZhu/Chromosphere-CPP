#include "chromosphere.hpp"

#include <iostream>

namespace chromosphere {

// ============================================================================
// Time-step calculation
// ============================================================================

Vec cal_max_v_i(const Grid& grid, const Vec& xn_state) {
    Vec spc_r_state = cal_spectral_radius_state(grid, xn_state);
    return get_scalar(grid, spc_r_state, 0);
}

Vec cal_dt_i(const Grid& grid, const Vec& xn_state) {
    Vec dt_i = grid.CFL * grid.ds_i / cal_max_v_i(grid, xn_state);
    float dt_min_i = arma::min(dt_i);
    return dt_min_i * arma::ones<Vec>(grid.ns);
}

// ============================================================================
// Time integrators
// ============================================================================
//
// Both integrators mutate grid.dt_state because rhs_explicit_state reads it
// during the MUSCL prediction step. The pattern is: broadcast dt_i across all
// num_of_eq slots, then call the RHS.
//

static void broadcast_dt(Grid& grid, const Vec& dt_i) {
    grid.dt_state.zeros(grid.n_state);
    for (arma::uword i = 0; i < num_of_eq; ++i) {
        grid.dt_state += scalar_to(grid, dt_i, i);
    }
}

Vec advance_Euler_state(Grid& grid, const Vec& xn_state, const Vec& dt_i) {
    const float tol = 1e-6f;
    broadcast_dt(grid, dt_i);

    Vec RE_n_state = rhs_explicit_state(grid, xn_state);
    Vec RI_np1_kp1_state(arma::size(xn_state), arma::fill::zeros);
    Vec RI_np1_k_state  (arma::size(xn_state), arma::fill::zeros);
    Vec xnp1_state;
    Vec error, norm;
    arma::uword it_cnt = 0;

    // Picard iteration on R_I (writeup §3.7) diverges whenever
    // dt * |dR_I/du| > 1 — i.e., when ion-neutral drag (α·dt ~ 60 at Model C7
    // conditions) or heat conduction makes R_I stiff. A true implicit solve
    // (Newton on the R_I Jacobian) is required to enable this branch.
    // Until then R_I stays zero, so the step is pure explicit Euler.
    do {
        RI_np1_k_state = RI_np1_kp1_state;
        xnp1_state = xn_state + grid.dt_state % (RE_n_state + RI_np1_k_state);
        // RI_np1_kp1_state = rhs_implicit_state(grid, xnp1_state);
        RI_np1_kp1_state.zeros(arma::size(xn_state));
        error = (RI_np1_kp1_state - RI_np1_k_state) % (RI_np1_kp1_state - RI_np1_k_state);
        norm  = RI_np1_k_state % RI_np1_k_state;
        ++it_cnt;
    } while (arma::sum(error) / arma::sum(norm) > tol && it_cnt < 25);

    return xnp1_state;
}

Vec advance_RK4(Grid& grid, const Vec& xn_state, const Vec& dt_i) {
    broadcast_dt(grid, dt_i);
    const Vec k1 = grid.dt_state % rhs_explicit_state(grid, xn_state);
    const Vec k2 = grid.dt_state % rhs_explicit_state(grid, xn_state + 0.5 * k1);
    const Vec k3 = grid.dt_state % rhs_explicit_state(grid, xn_state + 0.5 * k2);
    const Vec k4 = grid.dt_state % rhs_explicit_state(grid, xn_state + k3);
    return xn_state + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
}

// ============================================================================
// Debug
// ============================================================================

void print_xn(const Grid& grid, const Vec& xn) {
    for (arma::uword i = 0; i < grid.ns; ++i) {
        for (arma::uword j = 0; j < num_of_eq; ++j) {
            std::cout << xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, j)) << " ";
        }
        std::cout << std::endl;
    }
}

} // namespace chromosphere
