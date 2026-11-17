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

// Coarse-equivalent cell count of the release column: CHROMO_MESH_NS, default
// 600 (the model_column preset pins it to 500).
arma::uword column_peek_ns() {
    arma::uword ns = 600;
    if (const char* e = std::getenv("CHROMO_MESH_NS")) { try { ns = std::stoul(e); } catch (...) {} }
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
        // launcher-specific collection of assignments. No two-fluid,
        // fixed-gamma, finite-rate-ionization, cooling, TRAC or artificial-
        // conduction knob appears here: the release solver has none of them,
        // and model_column_ic rejects every retired ISO_* name outright.
        set_env_default("CHROMO_EOS_TABLE", "data/eos/gamma1_hydrogen_v1.dat");
        set_env_default("CHROMO_DOMAIN_BASE_KM",          "1600");
        set_env_default("CHROMO_DOMAIN_SPAN_KM",          "553");
        set_env_default("CHROMO_MESH_NS",                 "500");
        set_env_default("CHROMO_CONDUCTION",              "1");
        set_env_default("CHROMO_BC_T_TOP_K",              "22000");
        set_env_default("CHROMO_MESH_REFINE_PROFILE",     "outer");
        set_env_default("CHROMO_MESH_REFINE_FACTOR",      "4");
        set_env_default("CHROMO_MESH_REFINE_START_KM",    "500");
        set_env_default("CHROMO_MESH_REFINE_GRADE_KM",    "20");
        sc.peek_ns   = []() { return column_peek_ns(); };
        sc.ic        = [](Grid& g) { return model_column_ic(g); };
        sc.update_bc = [](Grid& g, const Vec& xn) { model_column_update_bc(g, xn); };
        return sc;
    }
    throw std::runtime_error("make_scenario: unknown scenario name '" + name + "'");
}

} // namespace chromosphere
