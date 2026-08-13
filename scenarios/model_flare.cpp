/*!
 * @file scenarios/model_flare.cpp
 * @brief Explosive-evaporation flare scenario: the C7 column plus a transient
 *        nonthermal electron beam and an open supersonic outer boundary.
 * @ingroup scenarios
 */
#include "model_flare.hpp"
#include "model_c7.hpp"

#include <cstdlib>
#include <string>

namespace chromosphere {

Vec model_flare_ic(Grid& grid) {
    // Force the full physical model: the explosive/gentle threshold is the
    // competition between beam heating and radiative loss, so Stage R cooling
    // must be active; ionization tracks the evaporated material going ionized.
    grid.enable_ionization        = true;
    grid.enable_radiative_cooling = true;
    // The full Stage E network (S_i + S_CR ionization, α_r + α_c recombination)
    // is the default (see chromosphere.hpp). The flare specifically *requires*
    // three-body recombination α_c = κ_c n_e: explosive evaporation shock-
    // compresses the condensation to densities where α_c (volumetric rate ∝
    // n_e³) is the dominant recombination sink; without it the dense, rapidly-
    // ionizing condensation has no adequate sink — the ionization fraction
    // drifts, the ion/neutral pressure partition corrupts, and the hydro goes
    // NaN a few s after the beam fires. model_c7_ic inverts P_phot over exactly
    // the active network, so C7 stays a fixed point regardless of the flags.

    // Build the C7 chromosphere + lower-TR column (also sets the inner discrete-
    // HSE base, the ambient RTV q(T), numerical diffusivity, TRAC). model_c7_ic
    // reads enable_radiative_cooling to gate q(T)/TRAC, hence the order above.
    Vec xn = model_c7_ic(grid);

    // --- Flare beam-heating driver (docs/explosive_evaporation_plan.md) ------
    // Beam energy flux F_e [W/m²]. Default 5×10⁷ W m⁻² = 5×10¹⁰ erg cm⁻² s⁻¹ —
    // Fisher et al. (1985) middle explosive run (well above F_crit ≈ 7×10⁶ W m⁻²,
    // peaks at ~500 km/s). Overridable via FLARE_BEAM_FLUX so the gentle control
    // (1×10⁶ W m⁻² = 1×10⁹ erg cm⁻² s⁻¹, below threshold) reuses this build.
    auto env_f = [](const char* key, float fallback) -> float {
        if (const char* e = std::getenv(key)) {
            try { return std::stof(std::string(e)); } catch (...) {}
        }
        return fallback;
    };

    float beam_flux = env_f("FLARE_BEAM_FLUX", 5.0e7f);

    grid.enable_beam_heating    = true;
    grid.beam_flux              = beam_flux;
    grid.beam_t_on              = env_f("FLARE_T_ON",   2.0f);   // let IC settle ~1 acoustic time
    grid.beam_duration          = env_f("FLARE_DUR",   10.0f);   // flat-top heating (Fisher: a few s)
    grid.beam_ramp              = env_f("FLARE_RAMP",   1.0f);    // cosine ramps
    grid.beam_h_lo_km           = env_f("FLARE_H_LO", 1600.0f);  // deposition window [km]
    grid.beam_h_hi_km           = env_f("FLARE_H_HI", 2000.0f);  // (upper chromosphere, below TR)

    // Open the outer face so the supersonic evaporation upflow can leave the
    // domain (the C7 Mach-0.05 cap would choke it). Env-toggle for diagnosis.
    grid.outer_free_outflow     = (env_f("FLARE_FREE_OUTFLOW", 1.0f) != 0.0f);

    // Optional extra numerical diffusivity to stabilize the evaporation shock on
    // the coarse grid (multiplier on the C7 grid-scaled value set by model_c7_ic).
    grid.numerical_diffusivity_per_length *= env_f("FLARE_DIFF_MULT", 1.0f);

    return xn;
}

void model_flare_update_bc(Grid& grid, const Vec& xn) {
    // Same closures as C7; model_c7_update_bc honors grid.outer_free_outflow.
    model_c7_update_bc(grid, xn);
}

} // namespace chromosphere
