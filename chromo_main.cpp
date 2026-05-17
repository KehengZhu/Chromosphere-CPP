#include "chromosphere.hpp"
#include "scenarios/model_c7.hpp"

#include <armadillo>
#include <fstream>
#include <iostream>

int main() {
    using namespace chromosphere;

    const arma::uword n_cells = 100;
    const float       cfl     = 0.25f;

    Grid grid;
    grid.init(n_cells, cfl);
    Vec xn = model_c7_ic(grid);

    const float c_s_target = 2.0e4f; // m/s, ion sound speed scale (writeup §4)
    const float total_time = 10.0f * arma::sum(grid.ds_i) / c_s_target;

    std::ofstream logf("output.log");
    logf << "Total time = " << total_time << std::endl;
    std::cout << "Total time = " << total_time << std::endl;

    float time = 0.0f;
    int   step = 0;
    while (time < total_time && step < 10000) {
        Vec dt = cal_dt_i(grid, xn);
        const float dt_avg = arma::mean(dt);

        if (step % 100 == 0) {
            std::cout << "step = " << step
                      << "  dt = " << dt_avg
                      << "  time = " << time << std::endl;
        }

        model_c7_update_bc(grid, xn);
        xn = advance_Euler_state(grid, xn, dt);
        time += dt_avg;
        ++step;
    }

    // Output: line 1 = cumulative cell heights (km, offset by 1.003e3 to match Fortran);
    // remaining lines = num_of_eq conserved variables per cell.
    std::ofstream fout("output.txt");
    float cum = 0.0f;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        cum += grid.ds_i(i);
        fout << "  " << (cum + 1.003e3f);
    }
    fout << '\n';
    for (arma::uword i = 0; i < grid.ns; ++i) {
        for (arma::uword k = 0; k < num_of_eq; ++k) {
            fout << "  " << xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, k));
        }
        fout << '\n';
    }
    std::cout << "END! step=" << step << " time=" << time << std::endl;
    return 0;
}
