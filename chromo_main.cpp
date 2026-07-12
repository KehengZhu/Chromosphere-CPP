#include "chromosphere.hpp"
#include "scenarios/scenario.hpp"

#include <armadillo>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    using namespace chromosphere;

    // Usage: chromo_main [output_path] [mode] [ionization] [scenario] [data_path] [time_mult] [cooling]
    //   output_path:  defaults to "output.txt"
    //   mode:         "full" (semi-implicit, default) or "explicit" (R_I ≡ 0)
    //   ionization:   "ionization" (default — Stage E enabled) or "no-ionization"
    //   scenario:     "model_column" (default — the unified chromosphere→corona column;
    //                 accepts "model_isentropic" / "model_gentle" as aliases) |
    //                 "model_c7" | "model_flare" | "pfss_field_line" | "analytic_canopy"
    //   data_path:    required for tabulated scenarios (e.g. pfss_field_line);
    //                 ignored otherwise (pass "-" or "" for scenarios that don't use it)
    //   time_mult:    multiplier on the default total_time (default 1.0). Step
    //                 cap scales accordingly so a longer run is not truncated.
    //   cooling:      "cooling" (Stage R: CL2012 optically-thick + optically-thin
    //                 radiative cooling, enabled by DEFAULT) or "no-cooling"
    const std::string out_path      = (argc > 1) ? argv[1] : "outputs/output.txt";
    const std::string mode          = (argc > 2) ? argv[2] : "full";
    const std::string ioniz_arg     = (argc > 3) ? argv[3] : "ionization";
    const std::string scenario_name = (argc > 4) ? argv[4] : "model_column";
    const std::string data_path     = (argc > 5) ? argv[5] : "";
    const float       time_mult     = (argc > 6) ? std::stof(argv[6]) : 1.0f;
    const std::string cool_arg      = (argc > 7) ? argv[7] : "cooling";
    const bool explicit_only        = (mode == "explicit");
    const bool ionization_on        = (ioniz_arg != "no-ionization");
    const bool cooling_on           = (cool_arg != "no-cooling");

    float cfl = 0.25f;
    if (const char* e = std::getenv("CHROMO_CFL")) { try { cfl = std::stof(e); } catch (...) {} }

    Scenario sc;
    try {
        sc = make_scenario(scenario_name, data_path);
    } catch (const std::exception& e) {
        std::cerr << "scenario error: " << e.what() << std::endl;
        return 1;
    }

    Grid grid;
    grid.init(sc.peek_ns(), cfl);
    grid.enable_ionization = ionization_on;
    grid.enable_radiative_cooling = cooling_on;
    // "neutrals off" experiment: SINGLE_FLUID=1 slaves neutrals to the ion fluid
    // (single-fluid limit). Default (unset/0) is the full two-fluid model.
    if (const char* e = std::getenv("SINGLE_FLUID")) {
        try { grid.single_fluid = (std::stof(e) != 0.0f); } catch (...) {}
    }
    // Separate electron temperature T_e ≠ T_i (docs/electron_temperature_plan.md).
    // ENABLE_TE=1 switches on the three-temperature model; default (unset/0) is
    // the single-temperature baseline that reproduces the pre-T_e physics.
    if (const char* e = std::getenv("ENABLE_TE")) {
        try { grid.enable_Te = (std::stof(e) != 0.0f); } catch (...) {}
    }
    Vec xn = sc.ic(grid);

    const float c_s_target = 2.0e4f; // m/s, ion sound speed scale (writeup §4)
    const float total_time = time_mult * 10.0f * arma::sum(grid.ds_i) / c_s_target;
    const int   step_cap   = static_cast<int>(std::max(10000.0f, 10000.0f * time_mult));

    std::ofstream logf("outputs/output.log", std::ios::app);
    logf << "[" << out_path << " mode=" << mode << "] Total time = " << total_time << std::endl;
    std::cout << "[" << mode << "] Total time = " << total_time << std::endl;

    // Multi-snapshot output. Format:
    //   line 1: "ns num_of_eq"
    //   line 2: ns cumulative cell heights in km (offset by grid.out_base_km, the
    //           domain base — 1003 km for C7-based scenarios, 0 km for the
    //           photosphere-anchored model_isentropic C7 IC; ds_i is in m, hence
    //           the 1e-3 conversion)
    //   then, repeated: a "# t = T step = S" marker followed by ns lines of
    //   num_of_eq space-separated conserved-variable values.
    std::ofstream fout(out_path);
    fout << grid.ns << " " << num_of_eq << '\n';
    {
        float cum_km = 0.0f;
        for (arma::uword i = 0; i < grid.ns; ++i) {
            cum_km += grid.ds_i(i) * 1.0e-3f;
            fout << "  " << (cum_km + grid.out_base_km);
        }
        fout << '\n';
    }

    auto write_frame = [&](float t_now, int step_now) {
        fout << "# t = " << t_now << " step = " << step_now << '\n';
        for (arma::uword i = 0; i < grid.ns; ++i) {
            for (arma::uword k = 0; k < num_of_eq; ++k) {
                fout << "  " << xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, k));
            }
            fout << '\n';
        }
    };

    // Scale snapshot cadence with run length so that frame count stays
    // bounded (~500 frames per run) regardless of time_mult. Keeps render
    // time roughly constant when sweeping time_mult — animation length at
    // 30 fps is ~17 s independent of physical simulation duration.
    const int frame_stride = std::max(10,
        static_cast<int>(10.0f * std::max(1.0f, time_mult / 5.0f)));
    float     time         = 0.0f;
    int       step         = 0;
    write_frame(time, step);

    while (time < total_time && step < step_cap) {
        grid.sim_time = time;   // expose current time to time-dependent terms (beam window)
        Vec dt = cal_dt_i(grid, xn);
        const float dt_avg = arma::mean(dt);

        if (step % 100 == 0) {
            std::cout << "[" << mode << "] step = " << step
                      << "  dt = " << dt_avg
                      << "  time = " << time
                      << "  T_c = " << grid.trac_cutoff_T << std::endl;
        }

        sc.update_bc(grid, xn);
        xn = explicit_only ? advance_Euler_explicit_state(grid, xn, dt)
                           : advance_Euler_state(grid, xn, dt);
        time += dt_avg;
        ++step;

        if (step % frame_stride == 0) write_frame(time, step);
    }
    if (step % frame_stride != 0) write_frame(time, step);

    std::cout << "[" << mode << "] END! step=" << step << " time=" << time << std::endl;
    return 0;
}
