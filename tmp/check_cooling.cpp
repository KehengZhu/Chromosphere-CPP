// Standalone sanity check: evaluate the CL2012 cooling on a C7 column and
// print the height-integrated cooling rate. Compare to Morosin 2022:
//   <Q>_integrated ~ 22-28 kW/m^2 (height integral of Q_rad along z).
//
// Build: see tmp/check_cooling.sh

#include "chromosphere.hpp"
#include "physics.hpp"
#include "scenarios/model_c7.hpp"

#include <cstdio>

using namespace chromosphere;

int main() {
    Grid grid;
    grid.init(100, 0.25f);
    Vec xn = model_c7_ic(grid);

    Vec prim  = cons2prim(grid, xn);
    Vec rho_i = get_scalar(grid, prim, prim::RHO_I);
    Vec rho_n = get_scalar(grid, prim, prim::RHO_N);
    Vec p_i   = get_scalar(grid, prim, prim::P_I);
    Vec n_i   = rho_i / grid.m_i;
    Vec n_n   = rho_n / grid.m_n;
    Vec T_e   = p_i / (2.0f * n_i * grid.k_b);

    Vec Q = radiative_loss_thick(grid, n_i, n_n, T_e);

    // Height-integrated cooling [W/m^2]
    float Q_int = 0.0f;
    for (arma::uword i = 0; i < grid.ns; ++i) Q_int += Q[i] * grid.ds_i[i];

    std::printf("    i    T_e[K]    n_i[/m3]    n_n[/m3]    Q_rad[W/m3]\n");
    for (arma::uword i = 0; i < grid.ns; i += 10) {
        std::printf("  %3llu  %.3e  %.3e  %.3e  %.3e\n",
                    (unsigned long long)i, T_e[i], n_i[i], n_n[i], Q[i]);
    }
    std::printf("\n  height-integrated Q_rad = %.3e W/m^2 (%.2f kW/m^2)\n",
                Q_int, Q_int * 1.0e-3f);
    std::printf("  Morosin 2022 reference:  22-28 kW/m^2\n");

    return 0;
}
