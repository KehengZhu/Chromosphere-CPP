/*! Whole-grid CRASH zAv vs the actual Chromosphere2026 C++ Saha function. */
#include "eos.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

using chromosphere::saha_ionization_fraction_n_h;

namespace {

double analytic_gamma1(double n_h, double temperature) {
    const double x = saha_ionization_fraction_n_h(n_h,temperature);
    const double theta = chromosphere::eos_constants::chi_h
                       / (chromosphere::eos_constants::k_b*temperature);
    // Differentiate x^2/(1-x)=Phi(T)/n_H analytically in log coordinates.
    const double dx_dln_t = x*(1.0-x)*(1.5+theta)/(2.0-x);
    const double dx_dln_rho = -x*(1.0-x)/(2.0-x);
    const double chi_rho = 1.0+dx_dln_rho/(1.0+x);
    const double chi_t = 1.0+dx_dln_t/(1.0+x);
    const double rho = n_h*chromosphere::eos_constants::m_h;
    const double pressure = (1.0+x)*n_h*chromosphere::eos_constants::k_b
                          * temperature;
    const double cv = chromosphere::equilibrium_heat_capacity(rho,temperature);
    return chi_rho+pressure/(temperature*cv)*chi_t*chi_t;
}

struct WorstGamma1 {
    double error = 0.0;
    double temperature = 0.0;
    double n_h = 0.0;
    double analytic = 0.0;
    double crash = 0.0;
};

void update_worst(WorstGamma1& worst, double error, double temperature,
                  double n_h, double analytic, double crash) {
    if (error <= worst.error) return;
    worst = {error,temperature,n_h,analytic,crash};
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: eos_validate_saha RAW_CRASH_TABLE RUNTIME_GAMMA1_TABLE\n";
        return 2;
    }
    const std::string path = argv[1];
    chromosphere::EosGammaTable gamma_table;
    try {
        gamma_table = chromosphere::EosGammaTable::load(argv[2]);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
    std::ifstream input(path.c_str());
    if (!input) {
        std::cerr << "cannot open " << path << '\n';
        return 2;
    }

    std::size_t rows = 0, line_number = 0;
    double max_abs = 0.0, max_transition_rel = 0.0, max_tiny_log = 0.0;
    double max_gamma1_node_error = 0.0;
    WorstGamma1 worst_gamma1_abs, worst_gamma1_rel;
    std::string line;
    while (std::getline(input, line)) {
        ++line_number;
        const std::size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#') continue;
        std::istringstream parser(line);
        double value[11];
        for (double& item : value) {
            if (!(parser >> item)) {
                std::cerr << path << ':' << line_number << ": expected 11 numeric columns\n";
                return 2;
            }
        }
        std::string extra;
        if (parser >> extra) {
            std::cerr << path << ':' << line_number << ": extra column\n";
            return 2;
        }
        const double temperature = value[0];
        const double n_h = value[1];
        const double x_crash = value[3];
        const double gamma1_crash = value[7];
        if (!std::isfinite(temperature) || !std::isfinite(n_h)
            || !std::isfinite(x_crash) || !std::isfinite(gamma1_crash)
            || temperature <= 0.0 || n_h <= 0.0
            || x_crash < 0.0 || x_crash > 1.0) {
            std::cerr << path << ':' << line_number
                      << ": invalid T, n_H, zAv, or GammaS\n";
            return 2;
        }
        const double x_cpp = saha_ionization_fraction_n_h(n_h, temperature);
        max_abs = std::max(max_abs, std::abs(x_crash - x_cpp));
        if (x_cpp >= 1e-10 && x_cpp <= 1.0 - 1e-10)
            max_transition_rel = std::max(max_transition_rel,
                std::abs(x_crash - x_cpp) / x_cpp);
        if (x_cpp < 1e-8 && x_crash > 0.0)
            max_tiny_log = std::max(max_tiny_log,
                std::abs(std::log(x_crash) - std::log(x_cpp)));
        max_gamma1_node_error = std::max(max_gamma1_node_error,
            std::abs(gamma_table.gamma1(temperature,n_h)-gamma1_crash));
        const double gamma1_cpp = analytic_gamma1(n_h,temperature);
        const double gamma1_abs = std::abs(gamma1_cpp-gamma1_crash);
        const double gamma1_rel = gamma1_abs/std::abs(gamma1_crash);
        update_worst(worst_gamma1_abs,gamma1_abs,temperature,n_h,
                     gamma1_cpp,gamma1_crash);
        update_worst(worst_gamma1_rel,gamma1_rel,temperature,n_h,
                     gamma1_cpp,gamma1_crash);
        ++rows;
    }
    if (rows == 0) {
        std::cerr << path << ": no data rows\n";
        return 2;
    }

    std::cout << "C++ Saha vs CRASH zAv over " << rows << " rows:\n"
              << "  max_abs=" << max_abs
              << " transition_rel=" << max_transition_rel
              << " tiny_log=" << max_tiny_log << '\n'
              << "  runtime Gamma1 node max_abs=" << max_gamma1_node_error << '\n'
              << "  analytic Gamma1 vs CRASH GammaS max_abs="
              << worst_gamma1_abs.error << " at T=" << worst_gamma1_abs.temperature
              << " n_H=" << worst_gamma1_abs.n_h << " analytic="
              << worst_gamma1_abs.analytic << " CRASH=" << worst_gamma1_abs.crash << '\n'
              << "  analytic Gamma1 vs CRASH GammaS max_rel="
              << worst_gamma1_rel.error << " at T=" << worst_gamma1_rel.temperature
              << " n_H=" << worst_gamma1_rel.n_h << " analytic="
              << worst_gamma1_rel.analytic << " CRASH=" << worst_gamma1_rel.crash << '\n';
    if (max_abs > 1e-6 || max_transition_rel > 2e-3 || max_tiny_log > 2e-3) {
        std::cerr << "FAIL: C++ Saha / CRASH zAv tolerance exceeded\n";
        return 1;
    }
    if (max_gamma1_node_error > 1e-12) {
        std::cerr << "FAIL: runtime Gamma1 loader differs from raw CRASH nodes\n";
        return 1;
    }
    if (worst_gamma1_abs.error > 1e-3 || worst_gamma1_rel.error > 1e-3) {
        std::cerr << "FAIL: analytic Saha Gamma1 / CRASH GammaS tolerance exceeded\n";
        return 1;
    }
    std::cout << "  PASS: actual C++ Saha closure matches the CRASH grid\n";
    std::cout << "  PASS: actual C++ Gamma1 loader reproduces every table node\n";
    std::cout << "  PASS: independent analytic Saha Gamma1 matches CRASH GammaS\n";
    return 0;
}
