/*!
 * Parser for the tabulated scenario data-file format produced by
 * util/extract_field_line.py. Used by any scenario that initialises the grid
 * from an external file (pfss_field_line for now; Bifrost traces later).
 *
 * Format (ASCII, line-oriented):
 *
 *   # Comments are lines starting with '#'. Blank lines are skipped.
 *   [META]
 *     key=value     (one per line, free-form)
 *   [CELLS]
 *     i  ds_m  B_imh_T  B_iph_T  phi_g_imh  phi_g_iph  ne_im3  nn_im3  T_K
 *     (ns rows, indexed 0..ns-1 in column 1; whitespace-separated)
 *   [GHOSTS]
 *     tag  B_T  phi_g  ne  nn  T  T_e_factor
 *     (exactly 4 rows; tag in {outer_0, outer_1, inner_0, inner_1})
 *
 * Section order is fixed: [META] first (must declare ns=), then [CELLS],
 * then [GHOSTS]. peek_ns() only parses [META] for the ns= value.
 */

#pragma once

#include <string>

#include "chromosphere.hpp"

namespace chromosphere {

struct ScenarioDataFile {
    arma::uword ns          = 0;
    float       g_si        = 274.0f;
    float       B_outer_T   = 0.0f;
    float       phi_g_offset_Jpkg = 0.0f;
    // Field-line topology. "closed" ⇒ the outer face is the loop apex (reflecting
    // symmetry BC); anything else (default "open") ⇒ coronal outflow at the top.
    std::string topology    = "open";

    // CELLS, all length ns.
    Vec ds_m;
    Vec B_imh_T, B_iph_T;
    Vec phi_g_imh, phi_g_iph;
    Vec ne_im3, nn_im3, T_K;

    // GHOSTS, length 4 in order {outer_0, outer_1, inner_0, inner_1}.
    Vec ghost_B_T;
    Vec ghost_phi_g_Jpkg;
    Vec ghost_ne_im3, ghost_nn_im3, ghost_T_K;
    Vec ghost_T_e_factor;
};

/// Read just the ns= entry from the [META] section. Throws on missing.
arma::uword peek_ns_from_file(const std::string& path);

/// Read and parse the entire file. Throws std::runtime_error on any
/// structural problem (missing section, wrong row count, parse failure).
ScenarioDataFile parse_scenario_data_file(const std::string& path);

} // namespace chromosphere
