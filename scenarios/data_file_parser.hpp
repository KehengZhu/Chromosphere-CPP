/*!
 * @file scenarios/data_file_parser.hpp
 * @brief Parser for the tabulated [META]/[CELLS]/[GHOSTS] scenario data-file format.
 * @ingroup scenarios
 *
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

/** @addtogroup scenarios
 *  @{
 */

/// One fully parsed scenario data file: the [META] scalars, the ns per-cell
/// [CELLS] columns, and the four [GHOSTS] rows, in SI units throughout.
///
/// This is a faithful transcription of the file, not a Grid: the consumer
/// (pfss_ic) copies the geometry columns into the Grid and builds the conserved
/// state from (n_e, n_n, T). Vector members are sized only after a successful
/// parse_scenario_data_file().
struct ScenarioDataFile {
    /// Number of interior cells declared by `ns=` in [META]; the length of every
    /// [CELLS] column. Must match Grid::ns.
    arma::uword ns          = 0;
    /// Gravitational acceleration the generator assumed when it integrated the
    /// φ_g columns [m/s^2] (solar surface value 274). Recorded for provenance;
    /// the per-face φ_g columns are what the solver actually uses.
    float       g_si        = 274.0f;
    /// Field strength at the outer face of the last cell [T], i.e. B_iph(ns-1).
    /// Recorded for provenance; B_imh_T / B_iph_T are authoritative.
    float       B_outer_T   = 0.0f;
    /// Additive gauge offset already applied to the φ_g columns [J/kg]. The
    /// generators write 0, i.e. φ_g = 0 at the inner face of the trace.
    float       phi_g_offset_Jpkg = 0.0f;
    /// Field-line topology. "closed" ⇒ the outer face is the loop apex (reflecting
    /// symmetry BC); anything else (default "open") ⇒ coronal outflow at the top.
    std::string topology    = "open";
    /// Loop half-length L [m] (footpoint→apex arc), written by make_loop_dat.py.
    /// Used to calibrate the ambient coronal heating E_H0(L, s_H, T_max). 0 if the
    /// file omits it (older traces) — the heating overlay then falls back to the
    /// arc-length sum.
    float       loop_half_length_m = 0.0f;

    /// CELLS, all length ns.
    /// Cell width along the field line [m].
    Vec ds_m;
    /// Magnetic field magnitude at the inner (i−½) cell face [T].
    Vec B_imh_T;
    /// Magnetic field magnitude at the outer (i+½) cell face [T].
    Vec B_iph_T;
    /// Gravitational potential at the inner (i−½) cell face [J/kg], increasing
    /// outward along the trace.
    Vec phi_g_imh;
    /// Gravitational potential at the outer (i+½) cell face [J/kg].
    Vec phi_g_iph;
    /// Electron number density [m^-3]; also the ion (proton) number density
    /// under quasi-neutrality, so ρ_i = ne_im3 · m_i.
    Vec ne_im3;
    /// Neutral-hydrogen number density [m^-3], so ρ_n = nn_im3 · m_n.
    Vec nn_im3;
    /// Cell temperature [K], shared by ions, neutrals and electrons in the IC.
    Vec T_K;

    /// GHOSTS, length 4 in order {outer_0, outer_1, inner_0, inner_1}.
    /// Ghost magnetic field magnitude [T].
    Vec ghost_B_T;
    /// Ghost gravitational potential [J/kg] as written in the file. pfss_ic
    /// overrides it with the adjacent interior face potential (phi_g_iph(ns-1)
    /// at the outer end, phi_g_imh(0) at the inner end) to match the φ_g gauge
    /// the RHS assumes.
    Vec ghost_phi_g_Jpkg;
    /// Ghost electron/ion number density [m^-3].
    Vec ghost_ne_im3;
    /// Ghost neutral-hydrogen number density [m^-3].
    Vec ghost_nn_im3;
    /// Ghost temperature [K], before the T_e_factor scaling.
    Vec ghost_T_K;
    /// Dimensionless multiplier applied to ghost_T_K when the ghost energies are
    /// built (T_eff = ghost_T_K · ghost_T_e_factor). 1 leaves the tabulated
    /// temperature alone; model_c7 uses 2 at the outer ghost to mimic the hotter
    /// coronal electrons.
    Vec ghost_T_e_factor;
};

/// Read just the ns= entry from the [META] section. Throws on missing.
arma::uword peek_ns_from_file(const std::string& path);

/// Read and parse the entire file. Throws std::runtime_error on any
/// structural problem (missing section, wrong row count, parse failure).
ScenarioDataFile parse_scenario_data_file(const std::string& path);

/** @} */

} // namespace chromosphere
