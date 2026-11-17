/*!
 * @file scenarios/scenario.cpp
 * @brief Scenario name dispatch.
 * @ingroup scenarios
 */
#include "scenario.hpp"
#include "model_column.hpp"
#ifdef CHROMO_ENABLE_LEGACY
#include "legacy/scenarios/legacy_scenarios.hpp"
#endif

#include <stdexcept>
#include <cstdlib>
#include <string>

namespace {
// Set an environment variable only if the user has not already set it — used to
// apply a scenario "preset" (default env knobs) that the user can still override.
void set_env_default(const char* key, const char* value) {
    setenv(key, value, /*overwrite=*/0);
}

// Grid resolution for the column scenarios: a resolved corona (ISO_CORONA, a
// two-fluid-only option the release rejects) needs more cells than the
// chromosphere-only default. ISO_NS overrides, and the model_column release
// preset pins it to 500.
arma::uword column_peek_ns() {
    auto set = [](const char* k) {
        const char* e = std::getenv(k);
        return e && std::string(e) != "0" && !std::string(e).empty();
    };
    arma::uword ns = set("ISO_CORONA") ? 1000 : 600;
    if (const char* e = std::getenv("ISO_NS")) { try { ns = std::stoul(e); } catch (...) {} }
    return ns;
}

// The frozen scenarios under legacy/scenarios/: the two-fluid ones and the
// two-temperature model_column_2t. They are compiled only with
// -DCHROMO_ENABLE_LEGACY=ON; without it, asking for one fails immediately.
bool is_legacy_scenario(const std::string& name) {
    return name == "model_c7" || name == "model_flare" || name == "model_gentle"
        || name == "analytic_canopy" || name == "pfss_field_line"
        || name == "model_column_2t";
}
} // namespace

namespace chromosphere {

Scenario make_scenario(const std::string& name, const std::string& data_path) {
    Scenario sc;
    sc.name = name;

    if (is_legacy_scenario(name)) {
#ifdef CHROMO_ENABLE_LEGACY
        return make_legacy_scenario(name, data_path);
#else
        (void)data_path;
        throw std::runtime_error(
            "make_scenario: '" + name + "' is a legacy scenario, and this "
            "build does not include the legacy code. Reconfigure with "
            "-DCHROMO_ENABLE_LEGACY=ON to build it (see legacy/README.md).");
#endif
    }
    // The canonical field-aligned chromosphere→corona column (the DEFAULT
    // scenario). It always runs the single-fluid release solver; its historical
    // two-fluid sibling `model_gentle` is a separate legacy scenario.
    if (name == "model_column") {
        // Canonical release model: the validated Saha/Gamma1 N=500/R4
        // physical-conduction-only single-fluid column. The Gamma1 table is
        // supplied UNCONDITIONALLY (only the path is overridable), so
        // `model_column` names exactly one solver. The remaining entries are
        // override-preserving defaults for the release model/grid, so
        // selecting model_column means the release configuration without a
        // launcher-specific collection of ISO_* assignments. No two-fluid,
        // fixed-gamma, finite-rate-ionization, cooling, TRAC or artificial-
        // conduction knob appears here: the release solver has none of them,
        // and model_column_ic rejects them outright if they are set.
        if (std::getenv("ISO_GAMMA"))
            throw std::runtime_error(
                "make_scenario: model_column is the single-fluid release "
                "scenario and derives every thermodynamic index from the "
                "Saha/Gamma1 closure; ISO_GAMMA (the historical fixed-gamma "
                "two-fluid index) has no meaning here. Use the legacy "
                "model_gentle scenario for the historical fixed-gamma column.");
        set_env_default("GAMMA_TABLE", "data/eos/gamma1_hydrogen_v1.dat");
        set_env_default("ISO_H_BASE",                     "1600");
        set_env_default("ISO_DH",                         "553");
        set_env_default("ISO_NS",                         "500");
        set_env_default("ISO_HEAT_FLUX",                  "1");
        set_env_default("ISO_T_TOP",                      "22000");
        set_env_default("ISO_REFINE_PROFILE",             "outer");
        set_env_default("ISO_REFINE_FACTOR",              "4");
        set_env_default("ISO_REFINE_S_LO_KM",             "500");
        set_env_default("ISO_REFINE_TRANSITION_KM",       "20");
        sc.peek_ns   = []() { return column_peek_ns(); };
        sc.ic        = [](Grid& g) { return model_column_ic(g); };
        sc.update_bc = [](Grid& g, const Vec& xn) { model_column_update_bc(g, xn); };
        return sc;
    }
    throw std::runtime_error("make_scenario: unknown scenario name '" + name + "'");
}

} // namespace chromosphere
