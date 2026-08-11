/*! Smoke-test the tracked production Gamma1 table independently of CRASH. */
#include "eos.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require_close(double actual, double expected, double tolerance,
                   const char* description) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance)
        throw std::runtime_error(std::string(description) + " mismatch");
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: eos_validate_production_table TABLE\n";
        return 2;
    }
    try {
        const chromosphere::EosGammaTable table =
            chromosphere::EosGammaTable::load(argv[1]);
        if (table.temperature_size() != 501 || table.density_size() != 57)
            throw std::runtime_error("production dimensions mismatch");
        require_close(table.min_temperature(), 3.2e3, 2e-10, "minimum T");
        require_close(table.max_temperature(), 1.0e8, 2e-6, "maximum T");
        require_close(table.min_n_h(), 1.0e12, 2e-2, "minimum n_H");
        require_close(table.max_n_h(), 1.0e26, 2e12, "maximum n_H");

        // Stage-9 limiting states: both the analytic caloric index and the
        // independently tabulated CRASH sound-speed index recover a monatomic
        // gas when hydrogen is either neutral or fully ionized.
        const double five_thirds = 5.0/3.0;
        const chromosphere::GammaState neutral = chromosphere::gamma_state(
            table, table.max_n_h()*chromosphere::eos_constants::m_h,
            table.min_temperature());
        const chromosphere::GammaState ionized = chromosphere::gamma_state(
            table, table.min_n_h()*chromosphere::eos_constants::m_h,
            table.max_temperature());
        if (!(neutral.x_eq < 1.0e-8) || !(ionized.x_eq > 1.0-1.0e-12))
            throw std::runtime_error("neutral/ionized limiting composition mismatch");
        require_close(neutral.gamma_energy, five_thirds, 2e-8,
                      "neutral-limit gamma_energy");
        require_close(neutral.gamma_sound, five_thirds, 2e-12,
                      "neutral-limit Gamma1");
        // Even at 1e8 K the chosen energy zero retains the already-paid
        // ionization potential, so gamma_E approaches 5/3 as O(chi/kT).
        require_close(ionized.gamma_energy, five_thirds, 4e-4,
                      "ionized-limit gamma_energy");
        require_close(ionized.gamma_sound, five_thirds, 2e-12,
                      "ionized-limit Gamma1");

        // Bilinear interpolation must meet continuously on every internal grid
        // line. Probe each line from both adjacent cells in log coordinates.
        double max_continuity_jump = 0.0;
        for (std::size_t j = 0; j < table.density_size(); ++j) {
            const double n_h = table.n_h_at(j);
            for (std::size_t i = 1; i+1 < table.temperature_size(); ++i) {
                const double t = table.temperature_at(i);
                const double eps = 1.0e-10*std::min(
                    std::log(t/table.temperature_at(i-1)),
                    std::log(table.temperature_at(i+1)/t));
                const double left = table.gamma1(t*std::exp(-eps), n_h);
                const double right = table.gamma1(t*std::exp(eps), n_h);
                max_continuity_jump = std::max(max_continuity_jump,
                                               std::abs(left-right));
            }
        }
        for (std::size_t i = 0; i < table.temperature_size(); ++i) {
            const double t = table.temperature_at(i);
            for (std::size_t j = 1; j+1 < table.density_size(); ++j) {
                const double n_h = table.n_h_at(j);
                const double eps = 1.0e-10*std::min(
                    std::log(n_h/table.n_h_at(j-1)),
                    std::log(table.n_h_at(j+1)/n_h));
                const double below = table.gamma1(t, n_h*std::exp(-eps));
                const double above = table.gamma1(t, n_h*std::exp(eps));
                max_continuity_jump = std::max(max_continuity_jump,
                                               std::abs(below-above));
            }
        }
        if (max_continuity_jump > 2.0e-10)
            throw std::runtime_error("Gamma1 interpolation is discontinuous at a grid line");

        const double t0 = 4.197356881658640;
        const double t1 = 4.206346581702001;
        const double n0 = 20.0;
        const double n1 = 20.25;
        require_close(table.gamma1(std::pow(10.0, t0), std::pow(10.0, n0)),
                      1.656361479772092, 2e-13, "representative node");
        require_close(table.gamma1(std::pow(10.0, 0.5*(t0+t1)),
                                   std::pow(10.0, 0.5*(n0+n1))),
                      1.6541815739025656, 2e-13, "representative interpolation");

        const double densities[] = {1.0e12, 1.0e16, 1.0e20, 1.0e26};
        const double temperatures[] = {3.2e3, 5.0e3, 1.0e4, 1.0e5, 1.0e7, 1.0e8};
        for (double n_h : densities) {
            const double rho = n_h * chromosphere::eos_constants::m_h;
            for (double temperature : temperatures) {
                const double energy =
                    chromosphere::equilibrium_internal_energy(rho, temperature);
                const double recovered = chromosphere::temperature_from_rho_eint(
                    table, rho, energy);
                require_close(recovered / temperature, 1.0, 2e-11,
                              "production caloric round trip");
                const chromosphere::GammaState gamma =
                    chromosphere::gamma_state(table, rho, temperature);
                if (!(gamma.x_eq >= 0.0 && gamma.x_eq <= 1.0)
                    || !(gamma.gamma_energy > 1.0)
                    || !(gamma.gamma_sound > 1.0))
                    throw std::runtime_error("production GammaState bounds mismatch");
            }
        }
        const double n_h = 1.0e20;
        const double rho = n_h * chromosphere::eos_constants::m_h;
        const double temperature = 1.0e4;
        // Conserved-state round trip of the release representation: pack
        // (rho, rho u, E) from (rho, u, T) and decode it back.
        const double velocity = 1.5e4;
        const double phi = 3.0e6;
        const double e_int =
            chromosphere::equilibrium_internal_energy(rho, temperature);
        const double energy = chromosphere::mixture_total_energy(
            rho, rho*velocity, e_int, phi);
        const chromosphere::MixtureThermo decoded =
            chromosphere::decode_equilibrium_mixture(
                table, rho, rho*velocity, energy, phi);
        require_close(decoded.rho/rho, 1.0, 2e-14,
                      "production mixture density round trip");
        require_close(decoded.T/temperature, 1.0, 2e-12,
                      "production mixture temperature round trip");
        require_close(decoded.internal_energy/e_int, 1.0, 2e-14,
                      "production mixture internal-energy round trip");
        require_close(decoded.n_e/(decoded.x*n_h), 1.0, 2e-14,
                      "production derived electron density");
        require_close((decoded.n_e+decoded.n_HI)/n_h, 1.0, 2e-14,
                      "production derived carrier closure");
        std::cout << "PASS: tracked Gamma1 table dimensions, limits, all-grid-line continuity,\n"
                  << "      node/interpolation, 24 caloric round trips, and the mixture state round trip"
                  << " (max continuity jump=" << max_continuity_jump << ")\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
