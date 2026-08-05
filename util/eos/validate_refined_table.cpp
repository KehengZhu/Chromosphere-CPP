/*! Stage-9 off-node validation against a direct-CRASH 2x-refined grid. */
#include "eos.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct RawRow {
    double temperature, n_h, energy, gamma_energy, gamma1, cv;
};

std::vector<RawRow> read_raw(const std::string& path) {
    std::ifstream input(path.c_str());
    if (!input) throw std::runtime_error("cannot open " + path);
    std::vector<RawRow> rows;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input,line)) {
        ++line_number;
        const std::size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#') continue;
        std::istringstream parser(line);
        double value[11];
        for (double& item : value)
            if (!(parser >> item))
                throw std::runtime_error(path+":"+std::to_string(line_number)
                                         +": expected 11 columns");
        rows.push_back({value[0],value[1],value[5],value[6],value[7],value[10]});
    }
    return rows;
}

double relative_error(double actual, double expected) {
    return std::abs(actual-expected)/std::max(std::abs(expected),1.0e-300);
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: eos_validate_refined_table COARSE_RUNTIME REFINED_RUNTIME REFINED_RAW\n";
        return 2;
    }
    try {
        const chromosphere::EosGammaTable coarse =
            chromosphere::EosGammaTable::load(argv[1]);
        const chromosphere::EosGammaTable refined =
            chromosphere::EosGammaTable::load(argv[2]);
        const std::vector<RawRow> raw = read_raw(argv[3]);
        if (refined.temperature_size() != 2*coarse.temperature_size()-1
            || refined.density_size() != 2*coarse.density_size()-1
            || raw.size() != refined.temperature_size()*refined.density_size())
            throw std::runtime_error("refined grid is not exactly 2x the production intervals");

        double coarse_direct = 0.0, coarse_refined = 0.0, refined_direct = 0.0;
        double gamma_energy_error = 0.0, cv_error = 0.0, inversion_error = 0.0;
        std::size_t samples = 0;
        for (std::size_t j = 1; j+1 < refined.density_size(); j += 2) {
          for (std::size_t i = 1; i+1 < refined.temperature_size(); i += 2) {
            ++samples;
            const RawRow& point = raw[j*refined.temperature_size()+i];
            const double gamma_coarse = coarse.gamma1(point.temperature,point.n_h);
            const double gamma_refined = refined.gamma1(point.temperature,point.n_h);
            coarse_direct = std::max(coarse_direct,std::abs(gamma_coarse-point.gamma1));
            coarse_refined = std::max(coarse_refined,std::abs(gamma_coarse-gamma_refined));
            refined_direct = std::max(refined_direct,std::abs(gamma_refined-point.gamma1));

            const double rho = point.n_h*chromosphere::eos_constants::m_h;
            const chromosphere::GammaState analytic =
                chromosphere::gamma_state(coarse,rho,point.temperature);
            gamma_energy_error = std::max(gamma_energy_error,
                std::abs(analytic.gamma_energy-point.gamma_energy));
            const double cv_norm = chromosphere::equilibrium_heat_capacity(
                rho,point.temperature)/(point.n_h*chromosphere::eos_constants::k_b);
            cv_error = std::max(cv_error,relative_error(cv_norm,point.cv));
            const double recovered = chromosphere::temperature_from_rho_eint(
                coarse,rho,point.energy,point.temperature);
            inversion_error = std::max(inversion_error,
                                       relative_error(recovered,point.temperature));
          }
        }

        std::cout << "Stage-9 2x-grid/direct-CRASH validation over " << samples
                  << " exhaustive production-cell midpoints:\n"
                  << "  Gamma1 coarse-vs-direct max_abs=" << coarse_direct << '\n'
                  << "  Gamma1 coarse-vs-refined max_abs=" << coarse_refined << '\n'
                  << "  Gamma1 refined-node-vs-direct max_abs=" << refined_direct << '\n'
                  << "  gamma_E analytic-vs-direct max_abs=" << gamma_energy_error << '\n'
                  << "  Cv_eff analytic-vs-direct max_rel=" << cv_error << '\n'
                  << "  T(rho,e_CRASH) max_rel=" << inversion_error << '\n';
        if (coarse_direct > 6.0e-3 || coarse_refined > 6.0e-3
            || refined_direct > 1.0e-12 || gamma_energy_error > 3.0e-5
            || cv_error > 2.0e-3 || inversion_error > 3.0e-6)
            throw std::runtime_error("Stage-9 refined/off-node tolerance exceeded");
        std::cout << "  PASS: runtime between-node errors remain within Stage-9 gates\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
