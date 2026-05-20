#include "scenario.hpp"
#include "analytic_canopy.hpp"
#include "model_c7.hpp"
#include "pfss_field_line.hpp"

#include <stdexcept>

namespace chromosphere {

Scenario make_scenario(const std::string& name, const std::string& data_path) {
    Scenario sc;
    sc.name = name;

    if (name == "model_c7") {
        sc.peek_ns   = []() -> arma::uword { return 100; };
        sc.ic        = [](Grid& g) { return model_c7_ic(g); };
        sc.update_bc = [](Grid& g, const Vec& xn) { model_c7_update_bc(g, xn); };
        return sc;
    }
    if (name == "analytic_canopy") {
        sc.peek_ns   = []() { return analytic_canopy_peek_ns(); };
        sc.ic        = [](Grid& g) { return analytic_canopy_ic(g); };
        sc.update_bc = [](Grid& g, const Vec& xn) { analytic_canopy_update_bc(g, xn); };
        return sc;
    }
    if (name == "pfss_field_line") {
        if (data_path.empty())
            throw std::runtime_error("make_scenario: pfss_field_line requires a data_path "
                                     "(6th positional argument to chromo_main)");
        sc.peek_ns   = [data_path]() { return pfss_peek_ns(data_path); };
        sc.ic        = [data_path](Grid& g) { return pfss_ic(g, data_path); };
        sc.update_bc = [](Grid& g, const Vec& xn) { pfss_update_bc(g, xn); };
        return sc;
    }

    throw std::runtime_error("make_scenario: unknown scenario name '" + name + "'");
}

} // namespace chromosphere
