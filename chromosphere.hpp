/*!
 * @file chromosphere.hpp
 * @brief Shared solver infrastructure: the Grid, the two state-row index sets,
 *        and the reusable scratch / diagnostic-capture structures.
 * @ingroup shared
 *
 * Chromosphere model — SHARED solver infrastructure only.
 *
 * The project carries TWO solvers over one shared field-aligned mesh, each in
 * its own source directory:
 *
 *   * the RELEASE solver, `src/single_fluid/` (`single_fluid/mixture.hpp`): a
 *     single-fluid equilibrium-mixture model whose conserved state is the
 *     three-component vector U = (rho, rho u, E). This is what the canonical
 *     `model_column` production configuration advances;
 *
 *   * the HISTORICAL two-fluid / multi-temperature / finite-rate-ionization
 *     research solver, `src/two_fluid/` (`two_fluid/two_fluid.hpp`), whose
 *     conserved state is the seven-row ion/neutral/electron carrier vector.
 *     It is NOT the release path.
 *
 * This header declares only what BOTH sides share: the Grid (mesh geometry,
 * prescribed magnetic field, gravity, physical constants, runtime toggles),
 * the two solvers' state-row index sets, and the reusable scratch/diagnostic
 * structures the Grid owns. Neither solver's entry points are declared here,
 * and neither solver includes the other's header. The dependency direction is
 *
 *     chromosphere.hpp + eos.hpp  (shared)
 *            ^                ^
 *      single_fluid        two_fluid
 *
 * The two solvers do NOT share a conserved-state width, a packing convention,
 * or a timestep path.
 *
 * All solver state lives in a Grid struct that the caller owns; functions
 * take `const Grid&` (or `Grid&` for those that mutate scratch buffers).
 * No global mutable state — multiple Grids can coexist.
 *
 * This header contains code derived from CoMFi (MIT License).
 */

#pragma once

#include <armadillo>
#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>
#include "eos.hpp"

namespace chromosphere {

/** @addtogroup grid
 *  @{
 */

/// Armadillo column vector — the workhorse type for cell-wise data, including
/// the conserved state and the static mesh metric arrays.
///
/// **The release storage type is `float`.** `CHROMO_STATE_FLOAT64` promotes it to
/// `double` and exists for ONE purpose: a controlled experiment isolating float32
/// representation round-off, which is a known and documented limitation of this
/// solver (the per-step mass-density increment falls below half a float32 ULP of
/// rho through most of the lower chromosphere, so the density there cannot
/// evolve). It is **not** a supported configuration, it is not exercised by the
/// release validation, and it must not be enabled for production or for any
/// quoted release number. Build it into a separate tree:
///
///     cmake -S . -B build_omp_f64 -DCMAKE_BUILD_TYPE=Release \
///           -DCHROMO_ENABLE_OPENMP=ON -DCHROMO_STATE_FLOAT64=ON
#ifdef CHROMO_STATE_FLOAT64
typedef arma::Col<double> Vec;
#else
typedef arma::Col<float> Vec;
#endif

// ============================================================================
// Release state indices — the single-fluid equilibrium-mixture solver
// ============================================================================

/// Conserved rows of the release state U = (rho, rho u, E).
///
///   RHO     total mass density                       [kg m^-3]
///   MOM     field-aligned momentum density rho u     [kg m^-2 s^-1]
///   ENERGY  total energy density
///           E = e_int(rho,T) + 1/2 rho u^2 + rho phi_g   [J m^-3]
///
/// e_int is the equilibrium internal energy INCLUDING the hydrogen ionization
/// energy, so the ionization fraction x, the electron density n_e and the
/// neutral density n_HI are derived thermodynamic quantities of (rho, e_int),
/// not independent conserved variables. See single_fluid/mixture.hpp.
namespace mix {
    const arma::uword RHO    = 0;  ///< total mass density rho                  [kg m^-3]
    const arma::uword MOM    = 1;  ///< field-aligned momentum density rho u    [kg m^-2 s^-1]
    const arma::uword ENERGY = 2;  ///< total energy density E                  [J m^-3]
}
/// Number of conserved rows in the release state; the packed state length is
/// `ns * num_of_mixture_eq`.
const arma::uword num_of_mixture_eq = 3;

// ============================================================================
// Historical two-fluid state indices — NON-RELEASE research solver
// ============================================================================

/// Number of conserved rows in the historical two-fluid state; the packed state
/// length is `ns * num_of_eq`. Deliberately different from num_of_mixture_eq —
/// the two solvers do not share a state width.
const arma::uword num_of_eq = 7;

/// Conserved-variable indices (writeup eq 61).
///
/// Three-temperature extension (docs/electron_temperature_plan.md): a 7th
/// variable E_E carries the electron internal energy. To keep the change a
/// clean, energy-conserving generalization, E_I retains its original meaning —
/// the TOTAL charged-fluid energy (protons + electrons + bulk KE + gravity) —
/// so the conservative MUSCL/Rusanov flux, gravity/area source, and sound speed
/// are byte-for-byte unchanged. The proton temperature is then a *derived*
/// quantity, T_i = (p_total − p_e)/(n_i k_B), and the electron temperature is
/// T_e = p_e/(n_i k_B) with p_e = ⅔ E_E. This "total-energy + electron-energy"
/// split is the standard two-temperature MHD formulation (AWSoM/BATSRUS;
/// Sokolov 2021) and guarantees ENABLE_TE=0 reproduces the single-temperature
/// baseline exactly (E_E never feeds back into E_I / momentum / ρ).
namespace cons {
    const arma::uword RHO_I = 0;  ///< ion mass density       ρ_i
    const arma::uword RHO_N = 1;  ///< neutral mass density   ρ_n
    const arma::uword MOM_I = 2;  ///< ion momentum density   ρ_i V
    const arma::uword MOM_N = 3;  ///< neutral momentum dens. ρ_n U
    const arma::uword E_I   = 4;  ///< TOTAL charged energy   e_i (protons+electrons+KE+φ)
    const arma::uword E_N   = 5;  ///< neutral total energy   e_n
    const arma::uword E_E   = 6;  ///< electron internal energy ε_e = (3/2) p_e
}

/// Primitive-variable indices. P_I is the TOTAL charged pressure (p_proton+p_e);
/// P_E is the electron partial pressure. Proton pressure = P_I − P_E (derived).
namespace prim {
    const arma::uword RHO_I = 0;  ///< ion mass density ρ_i                [kg m^-3]
    const arma::uword RHO_N = 1;  ///< neutral mass density ρ_n            [kg m^-3]
    const arma::uword V     = 2;  ///< ion velocity V                      [m s^-1]
    const arma::uword U     = 3;  ///< neutral velocity U                  [m s^-1]
    const arma::uword P_I   = 4;  ///< total charged pressure p_proton + p_e  [Pa]
    const arma::uword P_N   = 5;  ///< neutral pressure p_n                [Pa]
    const arma::uword P_E   = 6;  ///< electron partial pressure p_e       [Pa]
}

/// Per-Grid reusable storage for the release nonlinear conduction solve. This is
/// scratch only, never a thermodynamic cache; keeping it on the owning Grid
/// preserves reentrancy and avoids function-static state.
struct MixtureConductionScratch {
    std::vector<double> rho;          ///< total mass density of each cell        [kg m^-3]
    std::vector<double> e_old;        ///< internal energy density entering the solve [J m^-3]
    std::vector<double> temperature;  ///< current Newton iterate for the cell temperature [K]
    std::vector<double> target;       ///< backward-Euler right-hand side per cell [J m^-3]

    /// Per-cell quantities of the CURRENT Newton iterate. `e_at_T` is
    /// e_int(rho, temperature) of that iterate; the ionization fraction, the heat
    /// capacity, the conductivity inputs and this residual energy all come from one
    /// fused Saha evaluation per cell per pass (they used to be three independent
    /// Saha solves).
    ///   - `conductivity`  physical kappa_e + kappa_n            [W m^-1 K^-1]
    ///   - `capacity`      effective heat capacity d e_int/dT    [J m^-3 K^-1]
    ///   - `n_e` / `n_hi`  electron / neutral-hydrogen density   [m^-3]
    ///   - `e_at_T`        e_int at the iterate temperature      [J m^-3]
    std::vector<double> conductivity;  ///< physical kappa_e + kappa_n [W m^-1 K^-1]
    std::vector<double> capacity;      ///< effective heat capacity d e_int/dT [J m^-3 K^-1]
    std::vector<double> n_e;           ///< electron density [m^-3]
    std::vector<double> n_hi;          ///< neutral-hydrogen density [m^-3]
    std::vector<double> e_at_T;        ///< e_int at the iterate temperature [J m^-3]

    /// Face conductance, already carrying the flux-tube area factor and the face
    /// spacing, so the flux to a neighbour is the conductance times the temperature
    /// difference [W m^-3 K^-1].
    ///@{
    std::vector<double> g_left;   ///< conductance to the i-1 neighbour
    std::vector<double> g_right;  ///< conductance to the i+1 neighbour
    ///@}

    /// Tridiagonal Newton system solved for the temperature correction each pass.
    ///@{
    std::vector<double> a;      ///< sub-diagonal
    std::vector<double> b;      ///< diagonal
    std::vector<double> c;      ///< super-diagonal
    std::vector<double> rhs;    ///< residual of the backward-Euler equation [J m^-3]
    std::vector<double> delta;  ///< solved temperature correction [K]
    ///@}

    /// Size every buffer above for an `n`-cell solve.
    void resize(std::size_t n);
};

struct Grid;

/// Equilibrium thermodynamics decoded from one immutable packed release state
/// U = (rho, rho u, E). Physical cells and two explicitly decoded ghost layers
/// share the same state-carried gravitational-potential convention.
struct MixtureField {
    /// Identity of the state this field was decoded from. require_matches() checks
    /// all four, so a field can never be silently used against a different Grid, a
    /// different state object, a reallocated buffer, or a resized state.
    ///@{
    const Grid* source_grid = nullptr;     ///< Grid the decode was performed on
    const Vec* source_state = nullptr;     ///< state object that was decoded
    const Vec::elem_type* source_memory = nullptr;  ///< its backing buffer, as `Vec::memptr()`
    arma::uword source_elements = 0;       ///< its element count
    ///@}

    /// Cache generation of the decoded state. The driver bumps it once per step, so
    /// a predicted or post-source state gets its own generation and cannot be
    /// mistaken for the state at the start of the step.
    std::uint64_t state_generation = 0;

    /// Copy of the four ghost buffers (two inner, two outer, three rows each) as
    /// they stood at decode time. The ghosts are not part of the state vector, so
    /// require_matches() compares this signature to detect a boundary refresh that
    /// would invalidate the decoded ghost layers.
    std::array<float, 4*num_of_mixture_eq> boundary_signature{};

    /// Decoded equilibrium thermodynamics of every physical cell, index 0..ns-1.
    std::vector<MixtureThermo> cells;
    /// (ns+4) x 3 reconstruction variables (log rho, u, log p | log T).
    arma::Col<double> extended_primitive;
    /// Decoded temperature of every extended cell (ns+4, ghosts included), in the
    /// same indexing as extended_primitive. It is NOT a reconstruction variable:
    /// it exists only so the pressure-based face builder can seed its (rho,p)->T
    /// Newton iteration with the temperature of the cell each one-sided face
    /// state was extrapolated from.
    arma::Col<double> extended_temperature;

    /// Throw unless this field really is the decode of `state` on `grid`, ghosts
    /// included. Called at the top of every stage that consumes a decoded field, so
    /// a stale cache is a hard error rather than a silently wrong result.
    void require_matches(const Grid& grid, const Vec& state) const;
};

/// One side of every interface, as the release flux needs it. All three-row
/// arrays use the (ns, num_of_mixture_eq) packing of the conserved state.
struct MixtureFaceArrays {
    arma::Col<double> conserved;  ///< ns*3, (rho, rho u, E)
    arma::Col<double> flux;  ///< ns*3, (rho u, rho u^2 + p, (E+p) u)
    arma::Col<double> velocity;  ///< ns
    arma::Col<double> pressure;  ///< ns
    arma::Col<double> sound_speed;  ///< ns
    /// ns, the general-EOS pressure response at fixed total density, i.e. the
    /// partial derivative of p with respect to e_int at constant rho. Distinct from
    /// Gamma1: Gamma1 fixes the isentropic sound speed, this fixes the energy
    /// component of the Euler characteristic vectors.
    arma::Col<double> dp_deint;

    /// Size every array above for an `ns`-cell mesh.
    void resize(arma::uword ns);
};

/// Running tally of what the release `swmf_godunov` flux did, so a fallback
/// is diagnosable instead of silent. Accumulated over the whole run and reported
/// once at the end by the driver; also readable from a test. Counting is plain
/// (non-atomic) because the numerical-flux loop of mixture_rhs_explicit is serial.
struct GodunovFluxStats {
    std::uint64_t faces = 0;              ///< Faces the Godunov flux was asked for.
    std::uint64_t exact = 0;              ///< Faces the exact solve supplied.
    std::uint64_t fallback_bad_input = 0; ///< Fell back: a side state was not admissible.
    std::uint64_t fallback_vacuum = 0;    ///< Fell back: the vacuum guard tripped.
    std::uint64_t fallback_negative_p = 0;///< Fell back: a Newton step gave p* <= 0.
    std::uint64_t fallback_no_converge = 0;///< Fell back: the 10-iteration budget ran out.
    std::uint64_t fallback_bad_sample = 0;///< Fell back: the x/t=0 sample was not admissible.
    /// Largest number of Newton iterations any single face needed [-].
    std::uint64_t max_iterations = 0;
    /// Sum of Newton iterations over all exact solves, for a mean-cost estimate [-].
    std::uint64_t total_iterations = 0;

    /// Total number of faces that fell back to Rusanov.
    std::uint64_t fallbacks() const {
        return fallback_bad_input + fallback_vacuum + fallback_negative_p
             + fallback_no_converge + fallback_bad_sample;
    }
    /// Reset every counter to zero.
    void clear() { *this = GodunovFluxStats(); }
};

/// Reusable storage for the release MUSCL predictor/reconstruction. The numeric
/// arrays are deliberately untyped scratch; thermodynamic validity remains in
/// the state-bound MixtureField objects.
struct MixtureRhsScratch {
    /// Decode of the MUSCL-Hancock predictor state, carrying its own cache
    /// generation so it can never be confused with the state at the step start.
    MixtureField predicted;

    /// Untyped per-cell working arrays of the reconstruction (limited slopes,
    /// slope ratios, one-sided face primitives and their exponentiated values).
    std::array<arma::Col<double>, 24> work;

    /// The four one-sided face states of each cell, in the order
    /// {right-of-i+1/2, left-of-i+1/2, right-of-i-1/2, left-of-i-1/2}.
    std::array<MixtureFaceArrays, 4> faces;

    arma::Col<double> flux_iph;  ///< numerical flux at the upper face i+1/2, ns*3
    arma::Col<double> flux_imh;  ///< numerical flux at the lower face i-1/2, ns*3
    arma::Col<double> source;    ///< geometric pressure-area + gravity source, ns*3

    Vec predicted_state;  ///< packed conserved state after the Hancock half step
    Vec rhs;              ///< assembled explicit right-hand side returned to the caller

    /// Non-uniform MUSCL reconstruction weights and slope-ratio metrics. They are
    /// functions of the STATIC mesh alone (ds_i via ds_iph_i/ds_imh_i), so they are
    /// built once per mesh generation instead of once per timestep. On a uniform
    /// mesh they reduce to the constants of the classical uniform-grid formulas.
    ///@{
    arma::Col<double> W1, W3, W4;  ///< half-width weights of the owning cell
    arma::Col<double> metric_r;      ///< center-to-center gradient-ratio correction at cell i
    arma::Col<double> metric_r_ip1;  ///< the same correction evaluated at cell i+1
    arma::Col<double> metric_r_im1;  ///< the same correction evaluated at cell i-1
    ///@}

    /// Grid::metrics_generation() the weights above were built from; a mismatch
    /// means the mesh changed and they must be rebuilt.
    std::uint64_t weights_generation = 0;
    bool weights_valid = false;  ///< false until the weights have been built once
};

/// READ-ONLY capture of the release reconstruction and numerical mass flux for
/// ONE explicit RHS evaluation. Filled by mixture_rhs_explicit only when
/// Grid::capture_face_flux is true, and never read back by any solver stage — so
/// enabling it cannot change a numerical result. Index i refers to the UPPER face
/// i+1/2 of cell i, which is exactly the face the continuity row differences
/// (flux_iph); the cell-centred slots are the reconstruction INPUTS w of the same
/// evaluation, so a capture record is self-contained.
struct MixtureFaceFluxCapture {
    /** @name Cell-centred reconstruction inputs (length ns)
     *  The decoded state the reconstruction actually saw. */
    ///@{
    std::vector<double> rho_cell;  ///< total mass density   [kg m^-3]
    std::vector<double> v_cell;    ///< velocity             [m s^-1]
    std::vector<double> T_cell;    ///< temperature          [K]
    ///@}

    /** @name One-sided reconstructed face states at i+1/2
     *  `L` is extrapolated from cell i, `R` from cell i+1. */
    ///@{
    std::vector<double> rho_L, rho_R;  ///< face mass density   [kg m^-3]
    std::vector<double> v_L, v_R;      ///< face velocity       [m s^-1]
    std::vector<double> T_L, T_R;      ///< face temperature    [K]
    std::vector<double> cs_L, cs_R;    ///< face sound speed    [m s^-1]
    ///@}

    /** @name Face pressures [Pa]
     *  Authoritative total pressure of the same corrector face states, rebuilt
     *  through the production face builder on capture steps only. `p_L` versus
     *  `p_R` is the reconstruction-induced pressure mismatch that the
     *  `(log rho, V, log p)` reconstruction targets; `p_cell` is the decoded
     *  cell-centre value. */
    ///@{
    std::vector<double> p_cell;  ///< decoded cell-centre total pressure
    std::vector<double> p_L;     ///< face pressure of the state from cell i
    std::vector<double> p_R;     ///< face pressure of the state from cell i+1
    ///@}

    /** @name Numerical mass flux and its split [kg m^-2 s^-1]
     *  `a_face` is the EQUILIBRIUM acoustic spectral radius
     *  max(|V_L|+c_L, |V_R|+c_R) with c = sqrt(Gamma1 p/rho). It is the Rusanov
     *  coefficient when Rusanov is selected and the per-face fallback coefficient
     *  of both other solvers; under the release Godunov flux and under Roe it is
     *  otherwise diagnostic only. It is NOT the frozen signal speed the release
     *  Riemann solve and the CFL use.
     *  The total flux is decomposed into its central and dissipative parts. */
    ///@{
    std::vector<double> a_face;     ///< spectral radius        [m s^-1]
    std::vector<double> f_central;  ///< central (averaged-flux) part
    std::vector<double> f_diff;     ///< dissipative part
    std::vector<double> f_total;    ///< the flux actually differenced
    ///@}

    /** @name Limiter inputs and outputs actually used for this face
     *  Per reconstruction slot: `r` / `phi_plus` build the L state from cell i,
     *  `r_ip1` / `phi_minus` build the R state from cell i+1. Dimensionless. */
    ///@{
    std::vector<double> r_rho;          ///< log-density slope ratio at cell i
    std::vector<double> phi_plus_rho;   ///< limited log-density slope building the L state
    std::vector<double> r_ip1_rho;      ///< log-density slope ratio at cell i+1
    std::vector<double> phi_minus_rho;  ///< limited log-density slope building the R state
    std::vector<double> r_v;            ///< velocity slope ratio at cell i
    std::vector<double> phi_plus_v;     ///< limited velocity slope building the L state
    std::vector<double> r_T;            ///< thermal-slot slope ratio at cell i
    std::vector<double> phi_plus_T;     ///< limited thermal-slot slope building the L state
    ///@}

    /// True once a capture has been filled; false after resize() and whenever the
    /// capture was not armed for the step just taken.
    bool valid = false;

    /// Zero every array for an `ns`-cell mesh and clear `valid`.
    void resize(arma::uword ns) {
        std::vector<double>* all[] = {
            &rho_cell,&v_cell,&T_cell,
            &rho_L,&rho_R,&v_L,&v_R,&T_L,&T_R,&cs_L,&cs_R,
            &p_cell,&p_L,&p_R,
            &a_face,&f_central,&f_diff,&f_total,
            &r_rho,&phi_plus_rho,&r_ip1_rho,&phi_minus_rho,
            &r_v,&phi_plus_v,&r_T,&phi_plus_T};
        for (std::vector<double>* v : all) v->assign(ns, 0.0);
        valid = false;
    }
};

/// READ-ONLY capture of the OUTER-face (i = ns-1, upper face) quantities of the
/// FINAL CONVERGED Newton iteration of mixture_apply_conduction — the ones
/// that actually built g_right[ns-1] for the accepted solve. Filled only when
/// Grid::capture_outer_conduction is true and never read back by any solver
/// stage, so enabling it cannot change a numerical result.
///
/// The release outer-face conductive flux INTO the top cell is
///     q_face = kappa_face * (T_wall - T_top)/ds_face,
/// with kappa_face the PHYSICAL conductivity (kappa_e + kappa_n) — the release
/// conduction operator carries no TRAC broadening and no artificial diffusivity.
/// When the scenario imposes an external conductive reservoir temperature, that
/// datum sits at the PHYSICAL outer face, so
///     T_wall  = the imposed face temperature (e.g. 22 000 K),
///     kappa_face = kappa evaluated at that face state, and
///     ds_face = 1/2 ds_i(ns-1), the top-cell CENTRE-to-FACE distance.
/// Without the override the datum is the hydro ghost CENTRE and
/// ds_face = ds_iph_i(ns-1). area_ratio = B_i/B_iph is the flux-tube factor that
/// converts the face flux into the cell's volumetric divergence.
/// Positive q_face = heating the top cell.
struct OuterConductionCapture {
    double T_top  = 0.0;  ///< temperature at the top-cell centre [K]
    double T_wall = 0.0;  ///< imposed thermal datum at the outer face [K]
    double kappa_face = 0.0;  ///< physical face conductivity kappa_e + kappa_n [W m^-1 K^-1]
    double ds_face = 0.0;  ///< distance from the top-cell centre to the thermal datum [m]
    double area_ratio = 0.0;  ///< flux-tube factor B_i(ns-1)/B_iph(ns-1) [dimensionless]
    double q_face = 0.0;  ///< conductive flux into the top cell, positive = heating [W m^-2]
    bool   imposed_neumann = false;  ///< true when impose_outer_heat_flux replaced the Dirichlet datum
    bool   valid = false;  ///< true once a converged solve has filled this record
};

// ============================================================================
// Grid: owns all solver state.
// ============================================================================

/// The object that owns all solver state: mesh geometry, the prescribed magnetic
/// field, gravity, the physical constants, every runtime toggle, the ghost
/// buffers of both solvers, and the reusable scratch workspaces.
///
/// The caller owns the Grid. Solver functions take `const Grid&`, or `Grid&` when
/// they write a scratch buffer, so several Grids can coexist and every test can
/// construct its own — there is no global mutable state anywhere in the project.
///
/// Usage is always the same three steps: init() allocates for a cell count, the
/// scenario fills the geometry (and may resize() when it generates a locally
/// refined mesh), and broadcast() derives the packed and static-mesh caches from
/// it. Anything that edits `ds_i`, `B_i`, `B_imh`, `B_iph` or `dinvB_ds_i`
/// afterwards must call broadcast() again.
///
/// Which solver a Grid belongs to is decided by data, not by a flag: a non-empty
/// `eos_gamma_table` selects the release single-fluid solver, and each solver's
/// entry points reject the other's Grid. Fields below are therefore marked
/// RELEASE, LEGACY TWO-FLUID ONLY, or BOTH; a field belonging to the other solver
/// is ignored, not merely unused.
///
/// @see @ref grid, @ref architecture
struct Grid {
    /** @name Sizes and resolution */
    ///@{
    arma::uword ns              = 0;  ///< number of physical cells
    arma::uword n_mixture_state = 0;  ///< ns * num_of_mixture_eq  (release)
    arma::uword n_state         = 0;  ///< ns * num_of_eq          (legacy two-fluid)
    /// Courant number of the acoustic timestep limit. The conservative hard-coded
    /// default is 0.25; CHROMO_CFL=0.50 is the validated production value for the
    /// reduced model_column release configuration.
    float       CFL             = 0.25f;

    ///@}

    /** @name Physical constants (SI) */
    ///@{
    /// Fixed adiabatic index of the LEGACY two-fluid solver only. The release
    /// solver derives every thermodynamic index from the Saha/Gamma1 closure and
    /// never reads gamma_mono.
    float gamma_mono = 5.0f / 3.0f;
    /// CRASH Gamma1 table. A non-empty table selects the RELEASE single-fluid
    /// equilibrium-mixture solver (mixture.hpp); an empty one selects the legacy
    /// fixed-gamma two-fluid solver.
    EosGammaTable eos_gamma_table;
    /// Allow table lookups outside the tabulated density range by clamping to the
    /// nearest edge instead of throwing. A DEBUG bypass for offline validation
    /// tools; production runs leave it false so an out-of-domain state is a hard
    /// error rather than a silently extrapolated one.
    bool eos_gamma_debug_clamp = false;
    /// Per-cell temperature hint for the EOS inversion — the last T decoded in this
    /// cell by ANY stage. Purely an accelerator: it seeds Newton, never changes the
    /// converged root or the residual tolerance, and a stale or absent hint only
    /// costs the ordinary safeguarded solve. Measured on the h1600 ns=1000 column:
    /// a seeded call converges in 1.13 iterations with ZERO bisections, an unseeded
    /// one takes 8.68 iterations and 4.05 bisections. Mutable so the const-Grid
    /// decode paths (flux, conduction) can refresh it.
    mutable std::vector<double> eos_T_hint;

    /// Last decoded temperature in cell i, or NaN when none is recorded.
    double eos_temperature_hint(arma::uword i) const {
        return i < eos_T_hint.size() ? eos_T_hint[i]
                                     : std::numeric_limits<double>::quiet_NaN();
    }
    /// Record cell i's decoded temperature for the next inversion in that cell.
    void store_eos_temperature_hint(arma::uword i, double temperature) const {
        if (eos_T_hint.size() != ns)
            throw std::logic_error("EOS temperature hints were not pre-sized");
        if (i >= ns)
            throw std::out_of_range("EOS temperature hint cell is out of range");
        eos_T_hint[i] = temperature;
    }
    /// Adiabatic-index factors derived from gamma_mono so the equation of state
    /// and every energy↔pressure / heat-capacity conversion tracks a single γ.
    /// The thermal energy density of an ideal gas is ε = p/(γ−1), so the pressure
    /// is p = (γ−1)[E − ½ρv² − ρφ] and the heat capacity per volume is C = nk/(γ−1).
    /// At the γ=5/3 default these reduce to the historical literals: gm1()=2/3,
    /// inv_gm1()=3/2, half_gm1()=1/3. NOTE: these are the *thermodynamic* index
    /// only — Spitzer κ∝T^{5/2} and the ion–neutral collisional energy-exchange
    /// coefficient (3 k_B α/(m_i+m_n)) are kinetic-theory factors, intentionally
    /// NOT tied to gamma_mono.
    inline float gm1()      const { return gamma_mono - 1.0f; }           // γ−1
    inline float inv_gm1()  const { return 1.0f / (gamma_mono - 1.0f); }  // 1/(γ−1)
    inline float half_gm1() const { return 0.5f * (gamma_mono - 1.0f); }  // (γ−1)/2
    float m_i        = static_cast<float>(eos_constants::m_h);
    float m_n        = static_cast<float>(eos_constants::m_h);
    float m_e        = static_cast<float>(eos_constants::m_e);
    float g          = 0.27395e3f;
    float mu_0       = 4.0f * static_cast<float>(arma::datum::pi) * 1.0e-7f;
    float k_b        = static_cast<float>(eos_constants::k_b);
    float q_e        = 1.602176634e-19f;
    float chi_H_J    = static_cast<float>(eos_constants::chi_h);  ///< 13.6 eV

    ///@}

    /** @name Runtime toggles */
    ///@{
    /// LEGACY TWO-FLUID ONLY. Enable Stage E (hydrogen ionization / recombination)
    /// in advance_Euler_state. The release solver has no finite-rate ionization: its
    /// ionization state is the Saha equilibrium value of the EOS closure.
    bool enable_ionization = false;

    /// Sub-options for the Stage E ionization network (only active when
    /// enable_ionization = true). Both default ON: the full network is
    ///   ionization:    P_phot (photoionization) + S_CR (multilevel collisional)
    ///                  + S_i (Voronov ground-state electron-impact)
    ///   recombination: α_r (case-B radiative) + α_c = κ_c n_e (three-body)
    /// No single channel dominates across the parameter range this code spans
    /// (cool weakly-ionized chromosphere ~6 kK through the dense flare
    /// condensation to the ~10⁷ K, fully-ionized evaporated plasma), so the
    /// complete network is the correct default; each channel is the leading
    /// term in some regime (see visualization/c7_ioniz_recomb_rates.png):
    ///
    /// enable_direct_collisional_ionization: S_i (Voronov 1997 ground-state
    ///   electron-impact rate). Subdominant to the multilevel S_CR by 2–3 orders
    ///   in the cool chromosphere, but it is the *physically robust* ionization
    ///   channel above ~10⁵ K (transition region / flare), where the cool-plasma
    ///   S_CR fit is an extrapolation (S_i overtakes S_CR near 1.1×10⁵ K).
    ///
    /// enable_threebody_recombination: α_c = κ_c n_e (Hinnov & Hirschberg 1962 /
    ///   Stevefelt 1975). Its volumetric rate ∝ n_e³ makes it the dominant
    ///   recombination sink wherever the gas is compressed — the dense base and,
    ///   critically, the shock-compressed flare condensation, where dropping it
    ///   destabilizes the run. Enabling it makes the Stage E solve cubic rather
    ///   than quadratic. Negligible only in the tenuous upper chromosphere.
    bool enable_direct_collisional_ionization = true;
    bool enable_threebody_recombination       = true;

    /// TWO-FLUID ONLY. Enable optically-thick chromospheric radiative cooling
    /// (CL2012 recipe, writeup §2.4) as Stage R in the operator-split two-fluid
    /// integrator. The single-fluid release timestep has no radiative stage and
    /// ignores this flag entirely. Default off so the existing test suite remains
    /// a clean regression baseline; flip to true to activate.
    bool enable_radiative_cooling = false;

    /// LEGACY TWO-FLUID ONLY. Single-fluid limit ("neutrals off"). When true, the ion-neutral drag
    /// (Stage B) and T_i↔T_n equilibration (Stage C) rigidly slave the neutral
    /// fluid to the ions each step: V→U→V_cm (mass-weighted bulk velocity) and
    /// T_i→T_n→T_cm (heat-capacity-weighted bulk temperature). This is the exact
    /// infinite-coupling limit of the two-fluid equations — neutrals still carry
    /// mass and follow the ionization network, but have no independent dynamics,
    /// i.e. a single-fluid partially-ionized treatment. Default false = full
    /// two-fluid ("neutrals on"). Toggled at runtime by SINGLE_FLUID=1.
    bool single_fluid = false;

    /// LEGACY TWO-FLUID ONLY. Separate electron temperature T_e ≠ T_i
    /// (docs/electron_temperature_plan.md). The release solver is a common-temperature
    /// equilibrium mixture and never evolves a separate electron energy.
    /// The 7th conserved variable E_E (electron internal energy) is ALWAYS carried;
    /// this flag only decides whether the electrons evolve independently:
    ///   * enable_Te = false (default): the single-temperature baseline. Every
    ///     stage runs the original code path on the combined charged pool
    ///     (P_I = total charged pressure, "T_i" = T_charged = p_total/(2 n_i k_B)),
    ///     and at the end of each step E_E is slaved to ½ the charged thermal
    ///     energy so that T_e ≡ T_i is reported. E_E never feeds back, so this
    ///     reproduces the pre-T_e physics to round-off (cf. SINGLE_FLUID).
    ///   * enable_Te = true: the three-temperature model (T_e, T_i, T_n). The
    ///     beam deposits into electrons, Spitzer κ_e conducts on T_e, the
    ///     optically-thin radiation and the χ_H ionization cost are electron-pool
    ///     terms, and a 3-way point-implicit ν_ei/ν_en/ν_in relaxation (Stage C)
    ///     equilibrates the three temperatures. Decoupling (T_e ≠ T_i) appears at
    ///     the beam-driven onset and the tenuous loop-top where ν_ei is slow
    ///     (Bradshaw 2006; Manchester 2012). Toggled at runtime by ENABLE_TE=1.
    bool enable_Te = false;

    /// LEGACY TWO-FLUID ONLY. Well-balanced gravity in the explicit MUSCL
    /// reconstruction. The release reconstruction decodes each cell with its own
    /// phi_g, so it has no shifted-potential bias to correct.
    /// The TVD-MUSCL
    /// step limits PRIMITIVE variables, obtaining the neighbour pressures by
    /// cons2prim of the spatially-shifted conserved states (ip1/im1/ip2/im2 in
    /// rhs_explicit_state). But cons2prim subtracts the LOCAL-index φ_g, while the
    /// shifted cell's total energy carries that cell's ρφ_g — so each neighbour
    /// pressure is biased by ⅔ρ·Δφ_g = ⅔ρg·Δs, which turns the hydrostatic slope
    /// dp/ds = −ρg into −⅓ρg (the reconstructed pressure gradient is only ⅓ of the
    /// true one). A hydrostatic atmosphere then feels a spurious, resolution-
    /// independent (γ−1)g ≈ ⅔g downward force and cannot stay static. (Stage-D
    /// conduction already shifts φ_g correctly — rhs_implicit_state — so only the
    /// explicit reconstruction is affected.) When true, rhs_explicit_state corrects
    /// the shifted-state pressures with the shifted φ_g, so a hydrostatic profile is
    /// a discrete fixed point (V≈0) and the isentropic relaxation
    /// (docs/gentle_evaporation_downflow.md, model_column) holds. Default false
    /// keeps every existing scenario / the 5950-test baseline byte-for-byte; set by
    /// model_column. (Candidate to enable globally — it is the physically
    /// correct reconstruction and removes the spurious chromospheric downflow.)
    bool well_balanced = false;

    /// LEGACY TWO-FLUID ONLY. Log-space MUSCL reconstruction in rhs_explicit_state.
    /// The release reconstruction is always log-space by construction.
    /// In a gravitationally
    /// stratified atmosphere ρ and p are EXPONENTIAL in height, so a linear TVD slope
    /// misrepresents them and minmod clips the steep gradient to first order, injecting
    /// Rusanov numerical diffusion that is worst at the dense lower boundary. When true,
    /// the reconstruction limits log ρ / log p (the strictly-positive, exponentially-
    /// stratified primitive slots {RHO_I, RHO_N, P_I, P_N, P_E}) instead — an isothermal
    /// hydrostatic column is then piecewise-LINEAR, the minmod ratio r→1 passes the full
    /// slope, and the reconstructed face jump (hence the diffusion) collapses. The SIGNED
    /// velocities {V, U} stay linear (log undefined; no large dynamic range). Faces are
    /// exp'd back before prim2cons so the scheme stays conservative; applied AFTER the
    /// well-balanced φ_g pressure correction (which must see linear p). Default false
    /// keeps every existing scenario / the test baseline byte-for-byte; always on in
    /// model_column. Composes with well_balanced.
    bool log_reconstruct = false;

    /// Slope limiter for the explicit MUSCL reconstruction (rhs_explicit_state).
    /// Default false ⇒ the symmetric minmod limiter φ(r)=max(0,min(1,r)) used
    /// everywhere (writeup eq 15; slimL=slimR so one flux_lim(r) feeds both faces
    /// of a cell). True ⇒ the MC3 / Koren limiter (BATSRUS ModFaceValue 'mc3'), the
    /// ASYMMETRIC third-order (κ=1/3) monotonized-central limiter: less diffusive
    /// than minmod, so it clips the under-resolved TR gradient far less (the source
    /// of the O(Δs) base-drainage truncation in the iso runs). Because it is
    /// asymmetric, the '+' and '−' reconstruction terms take DIFFERENT limited
    /// slopes — see flux_lim_mc3_plus / flux_lim_mc3_minus. limiter_beta is the
    /// Koren β (2 ⇒ classic Koren; 1 ⇒ minmod-like). Default off keeps every
    /// existing scenario / the 5958-test baseline byte-for-byte; always on in
    /// model_column (β=2). Composes with log_reconstruct and well_balanced (it only
    /// changes how the same differences are limited).
    bool  mc3_limiter = false;
    float limiter_beta = 2.0f;

    /// REFERENCE numerical flux (comparison only, not the release default): the
    /// corrector dissipation computed from the 3x3 mixture Roe characteristic
    /// decomposition of U = (rho, rho u, E), linearized about the EQUILIBRIUM Saha
    /// system, so its acoustic eigenvalues carry the equilibrium index Gamma1. It is
    /// a general-EOS Roe-type linearization built from the tabulated EOS
    /// derivatives, not a Roe average constructed to satisfy the exact jump
    /// condition, so the classical Roe properties are not claimed for it. False
    /// leaves the face-local Rusanov/LLF solver, which is also the automatic
    /// per-face fallback whenever a Roe average is not admissible; Rusanov's
    /// acoustic-scale dissipation -1/2 a Delta U is far too large for this
    /// very-low-Mach evaporation problem and leaves a persistent TR velocity ripple.
    /// Selected by ISO_RIEMANN=roe-local; mutually exclusive with
    /// swmf_godunov_flux.
    bool roe_characteristic_flux = false;

    /// RELEASE numerical flux. True (the model_column release setting) is the
    /// SWMF-style Godunov flux built from the exact ideal-gas Riemann solver of
    /// `single_fluid/exact_rs.hpp`, with the non-ideal EOS content carried as the
    /// passive specific energy offset of SWMF `util/CRASH/src/test_godunov.f90`.
    /// The Riemann problem is solved at the FROZEN-COMPOSITION index gamma = 5/3:
    /// the equilibrium Saha internal energy splits as e_int = p/(5/3-1) + e_ion, so
    /// the translational gas is monatomic and the ionization reservoir rides in the
    /// offset, held fixed across the local wave interaction. It selects the flux
    /// INSTEAD of the Roe linearization, so it takes precedence over
    /// roe_characteristic_flux; model_column sets exactly one of the two, and every
    /// other stage of the timestep (reconstruction, Hancock predictor, source,
    /// boundaries, EOS, conduction) is independent of the choice, which makes the
    /// face Riemann solver the only variable it changes. Every face whose exact
    /// solve fails or yields an inadmissible sample falls back to the same
    /// face-local Rusanov flux the Roe path uses, and is counted in godunov_stats.
    /// It also sets the timestep: mixture_timestep sizes the step from the frozen
    /// signal speed sqrt(5/3 p/rho) whenever this flag is true. model_column_ic
    /// sets it on the Gamma/Saha path; the struct default keeps the two-fluid
    /// research scenarios on Rusanov.
    bool swmf_godunov_flux = false;

    /// RELEASE thermal reconstruction variable. True (the model_column release
    /// setting) limits
    ///   (log rho, u, log p)
    /// and recovers the face temperature by inverting the same authoritative
    /// closure (equilibrium_temperature_from_density_pressure). False limits
    ///   (log rho, u, log T)
    /// and obtains the face pressure afterwards from the nonlinear Saha mapping
    /// p(rho,T) — so two independently limited variables set one mechanical
    /// quantity, and a mechanically smooth (constant-p) state is NOT reproduced
    /// across the partial-ionization transition; it is retained as a reference
    /// configuration only. Both keep density and pressure positive by construction
    /// and feed the identical equilibrium face builder, flux and source stages, so
    /// this flag changes ONLY which pair of variables is authoritative at a face.
    bool pressure_reconstruct = false;

    /// Equilibrium-reference ("δ-form") well-balancing. A general NON-isothermal
    /// hydrostatic atmosphere (e.g. Model C7, with a temperature minimum) is only
    /// reproduced to O(Δs) by the reconstruction, so the explicit MUSCL step applies
    /// a small spurious force even at rest ⇒ a residual O(Δs) base drainage
    /// (docs/boundary_conditions_plan.md). This device removes it for ANY
    /// stratification: capture the scheme's explicit-RHS residual at a frozen
    /// reference equilibrium once, R_eq = RHS(eq_state), and subtract it every step.
    /// Then RHS(eq_state) − R_eq ≡ 0, so eq_state is an EXACT discrete steady state
    /// (u = 0 held to round-off, resolution-independent), while real deviations from
    /// equilibrium (conduction-driven evaporation) evolve normally minus that fixed
    /// O(Δs) correction. eq_state is the (u = 0, HSE) IC set by the scenario;
    /// eq_residual is computed lazily on the first step with subtraction disabled via
    /// the empty() guard. Always on in model_column; both solvers support it, each in
    /// its own conserved-state width.
    bool eq_wb = false;
    Vec  eq_state;  ///< frozen reference equilibrium (conserved), set by scenario IC
    Vec  eq_residual;  ///< cached explicit-RHS residual at eq_state (empty until computed)

    /// model_c7 "New explanation" upper BC (docs/gentle_evaporation_downflow.md),
    /// in its HYBRID form: when set, model_c7_update_bc imposes at the TR-base top
    /// face a hydrostatic ghost pressure (centered HSE ⇒ V = 0 is a boundary fixed
    /// point, so no spurious boundary downflow), a ghost temperature extrapolated
    /// CONTINUOUSLY from the live top cell (T_g0 = a·T_top, T_g1 = b·T_g0, with
    /// a = b = 1 by default; C7_TJUMP_A / C7_TJUMP_B are an experimentation lever),
    /// and an EOS ghost density at that (p, T) — replacing the Mach-capped RTV
    /// reservoir's ρ,T Neumann ghosts.
    ///
    /// The downward conductive flux q(T) keeps entering through the imposed Stage-D
    /// NEUMANN flux in both modes: impose_outer_heat_flux is deliberately left ON,
    /// because a Dirichlet coronal ghost temperature over-conducts on this coarse
    /// grid (κ_e ∝ T^{5/2}) and over-pressurizes the top cell. The ghost temperature
    /// here sets only the EOS ghost density and the hydro flux.
    ///
    /// Set true only by the bare `model_c7` scenario, which also turns on
    /// well_balanced; default false keeps model_flare / analytic_canopy /
    /// model_gentle — which share model_c7_ic / model_c7_update_bc — unchanged.
    bool c7_tr_jump_bc = false;

    /// BOTH solvers. Enable field-aligned heat conduction: the implicit physical
    /// conduction stage of the release timestep, and Stage D of the two-fluid
    /// advance_Euler_state. Default ON so the existing scenarios / test suite are
    /// byte-for-byte unchanged. A relaxation phase
    /// (docs/gentle_evaporation_downflow.md, model_column Stage 1) sets this false
    /// to switch OFF all heat flux, so the relaxation holds its profile with
    /// V ≈ 0; Stage 2 flips it back on, with the imposed top temperature supplying
    /// the downward conductive flux.
    bool enable_conduction = true;

    /// TWO-FLUID ONLY. Mesh-independent isotropic numerical-diffusion coefficient
    /// C_num [m/s]. The two-fluid Stage D constructs the actual diffusivity
    /// independently at each face as chi_num,f = C_num*Delta_s_f, using the real
    /// centre-to-centre face distance, then adds K_num,f = chi_num,f*C_V,f to the
    /// physical face conductivity, so it is also disabled when conduction is off.
    /// The release conduction operator is structurally PHYSICAL-only: it has no
    /// such term and never reads this field. Default 0.
    float numerical_diffusivity_per_length = 0.0f;

    /// Upper-boundary coronal conductive heat flux q(T) [W/m²], imposed as a
    /// Neumann (fixed-flux) outer BC on the charged-fluid conduction row (Stage D)
    /// when impose_outer_heat_flux is true. This is the RTV (1978) downward
    /// coronal conductive flux F_c = −κ₀T^{5/2}dT/ds re-radiated at the TR base —
    /// the heating that holds the chromosphere top against Route-B ionization
    /// (χ_H) cooling. A fixed *flux* (not a fixed hot temperature) is used because
    /// on the coarse grid a Dirichlet T-wall would over-conduct via κ_e∝T^{5/2}
    /// and over-pressurize the top cell. Positive = heating into the domain.
    /// Default off so existing tests keep the Dirichlet conduction BC.
    bool  impose_outer_heat_flux = false;
    float outer_heat_flux        = 0.0f;

    /// Conduction-only OUTER temperature override. Normally the outer thermal
    /// datum is whatever temperature the HYDRO outer ghost happens to carry, so one
    /// ghost state serves two roles at once: the hydro Riemann/reconstruction state
    /// AND the thermal wall that drives conduction. When
    /// outer_conduction_temperature_override is true the conduction rows instead
    /// use the explicit outer_conduction_temperature [K], letting a scenario relax
    /// the hydro ghost temperature (e.g. zero-gradient extrapolation from the live
    /// top cell) while the conductive driving stays the unchanged fixed hot wall.
    ///
    /// IMPORTANT — this also relocates the datum. The imposed temperature is a
    /// PHYSICAL-BOUNDARY-FACE value, not a ghost-centre value, so the release
    /// conduction operator uses the half-cell distance d = 1/2 ds_top and evaluates
    /// kappa at the face state reconstructed from the imposed external pressure and
    /// that face temperature. Without the override the datum stays at the ghost
    /// CENTRE, distance ds_iph_i(ns-1), with the ghost's own n_e / n_HI. No effect
    /// when impose_outer_heat_flux is true (that path has no outer Dirichlet datum).
    /// Default false ⇒ byte-identical to the pre-existing behaviour.
    bool  outer_conduction_temperature_override = false;
    float outer_conduction_temperature          = 0.0f;

    /// Inner-boundary conduction BC. Default false ⇒ Dirichlet: the innermost
    /// conduction face couples cell 0 to the inner ghost temperature (the
    /// photospheric reservoir's fixed T — a heat sink at the base). True ⇒ Neumann
    /// (zero conductive flux through the base, dT/ds = 0): the inner face coupling
    /// is dropped from the Stage-D tridiagonal so the base is insulating and the
    /// photosphere heats freely. Set by model_column from ISO_INNER_T_NEUMANN.
    bool  inner_conduction_neumann = false;

    /// Time-ramp on the imposed coronal conductive flux q(T) — the gentle
    /// conduction-driven evaporation driver (docs/gentle_evaporation_plan.md v2).
    /// The scenario stores the steady quiet value in outer_heat_flux_base; the
    /// boundary update sets outer_heat_flux = base · enhance(sim_time), where
    /// enhance(t) rises with a cosine ramp from 1 to outer_heat_flux_enhance over
    /// [t_on, t_on+ramp]. enhance = 1 (default) ⇒ no ramp (steady q, unchanged
    /// behavior for model_c7 / model_flare). Raising it above 1 increases the
    /// downward coronal conductive flux into the TR/upper chromosphere → the
    /// gentle upflow (Antiochos & Sturrock 1978). Set by model_gentle.
    float outer_heat_flux_base    = 0.0f;
    float outer_heat_flux_t_on    = 0.0f;
    float outer_heat_flux_ramp    = 1.0f;
    float outer_heat_flux_enhance = 1.0f;

    /// Corona-as-boundary domain top [km] for model_c7_ic(extended=true). The
    /// extended C7 table (Avrett & Loeser 2008 Table 26) reaches 68 Mm / 1.59 MK;
    /// this lifts the upper boundary into a real resolved coronal VOLUME (so the
    /// evaporated mass has somewhere to accumulate — the coronal EM rise that *is*
    /// gentle evaporation), with q(T) imposed at that high coronal top. 0 (default)
    /// ⇒ use the full table top. Set by model_gentle from GENTLE_TOP_KM.
    float corona_top_km           = 0.0f;

    /// Absolute height [km] of the domain base, used only to label the output
    /// heights chromo_main writes (cumulative ds + base). Defaults to the C7 base
    /// (1003 km) so every existing scenario's output is unchanged; model_column
    /// sets it to its actual base (0 km for the photosphere-anchored C7 IC).
    float out_base_km             = 1003.0f;

    /// TWO-FLUID ONLY. Transition-Region Adaptive Conduction (TRAC; Johnston et al.
    /// 2019, 2020). The single-fluid release conduction operator has no TRAC factor
    /// and ignores these fields entirely.
    /// When enabled, an adaptive cutoff temperature T_c is recomputed each step
    /// (the highest temperature where the TR is under-resolved, L_R/L_T > 1/2,
    /// bounded by [trac_T_chrom, 0.2 T_peak]); below T_c the parallel conductivity
    /// κ_e is broadened by ε=(T_c/T)^{5/2} and the optically-thin radiative loss
    /// (and any volumetric heating) divided by ε. This broadens the unresolved TR
    /// so a coarse grid resolves it, while keeping κΛ and κQ — hence the TR-
    /// integrated radiation/heating and the coronal density response — invariant
    /// (Johnston 2020 §2.2), which is the Eq.-27 jump-condition balance enforced
    /// analytically. Default off so existing tests are a clean baseline.
    bool  enable_trac   = false;
    float trac_T_chrom  = 2.0e4f;  ///< TRAC-region base temperature T_b [K]
    float trac_cutoff_T = 2.0e4f;  ///< adaptive cutoff T_c [K], recomputed each step
    float trac_Tc_max_frac = 0.2f;  ///< upper bound on T_c as a fraction of T_peak
                                    // (Johnston 2020 Eq. 9: 0.2). Raise toward 1 on a
                                    // corona-less domain (e.g. iso_t22k, T_peak≈22 kK)
                                    // so the cap does not pin T_c to the T_b floor.

    /// Effective hydrogen photoionization rate from the ground state [s^-1],
    /// modeling the Lyα-excitation + Balmer-continuum two-step channel
    /// (Carlsson & Stein 2002 establish this as the dominant chromospheric
    /// ionization path; the Lyman continuum itself contributes negligibly).
    /// C&S 2002 report a chromospheric ionization/recombination timescale of
    /// 10^3–10^5 s; we adopt the geometric mean ~10^4 s as a single-rate
    /// closure suitable for the 2-state (ground+continuum) reduction used
    /// here. See Leenaarts et al. 2007 for the "fixed radiative rates"
    /// implementation paradigm in non-equilibrium chromospheric MHD codes.
    ///
    /// The scalar `photoionization_rate` is the uniform fallback used when
    /// the cell-centered `photoionization_rate_i` array is empty (i.e. no
    /// height-dependent calibration has been performed by the scenario IC).
    /// Scenarios that want a height-dependent profile (e.g. model_c7_ic,
    /// which calibrates against C7's tabulated f to make C7 a fixed point
    /// of the Stage E quadratic) populate `photoionization_rate_i` of
    /// length `ns`; the Stage E driver picks that up automatically via
    /// `physics.hpp::photoionization_rate_P(grid)`.
    float photoionization_rate = 1.0e-4f;
    Vec   photoionization_rate_i;

    ///@}

    /** @name Flare beam heating (explosive chromospheric evaporation) */
    ///@{
    /// Transient nonthermal-electron-beam volumetric heating that drives
    /// chromospheric evaporation (Fisher, Canfield & McClymont 1985). The beam,
    /// injected from the corona at the top, deposits an energy flux beam_flux
    /// [W/m²] in the upper chromosphere over a layer beam_deposition_height [km]
    /// below the outer boundary, weighted ∝ local total density (the densest
    /// reachable layer absorbs the most — thick-target stopping; the true
    /// ~10²⁰ cm⁻² stopping column is sub-grid on the coarse C7 mesh). The flux is
    /// gated in time by a flat-top window [beam_t_on, beam_t_on+beam_duration]
    /// with cosine ramps of width beam_ramp. Explosive evaporation needs
    /// beam_flux above the Fisher threshold F_crit ≈ 7×10⁶ W m⁻² AND the
    /// radiative sink (Stage R) active, since the threshold is the competition
    /// between deposited heating and the upper-chromosphere radiative loss.
    /// Heating is added to the electron/ion thermal pool (p_i). TWO-FLUID ONLY:
    /// the release timestep has no beam stage and ignores this flag. Default off so
    /// existing scenarios/tests are unaffected.
    bool  enable_beam_heating     = false;
    float beam_flux               = 0.0f;  ///< F_e [W/m²]
    float beam_t_on               = 0.0f;  ///< beam onset time [s]
    float beam_duration           = 0.0f;  ///< flat-top duration τ [s]
    float beam_ramp               = 1.0f;  ///< cosine ramp half-width [s]
    /// Deposition height window [km] (absolute, C7 base = 1003 km). The beam
    /// stops in the upper chromosphere BELOW the TR base (~2150 km); depositing
    /// there (not in the TR cells at the open top boundary) is both physical
    /// (thick-target stopping in the dense chromosphere) and keeps the explosive
    /// overpressure interior so the evaporated upflow develops through the TR.
    float beam_h_lo_km            = 1600.0f;
    float beam_h_hi_km            = 2000.0f;

    /// Self-consistent thick-target beam (Emslie 1978). When true, the fixed
    /// [h_lo,h_hi] deposition window is ignored: the beam is injected at the loop
    /// APEX (the maximum of φ_g — loop top for closed/full, coronal top for open),
    /// streams down each leg, and loses energy collisionally. The volumetric
    /// heating at column depth N (from the apex) is Q ∝ n_tot·(1+N/N_c)^{-δ/2},
    /// normalized so ∫Q ds = beam_flux (thick target: all energy absorbed). N_c is
    /// the stopping column of cutoff-energy electrons (N_c ∝ E_c²). As the loop
    /// fills with evaporated plasma the column rises, so the stopping depth moves UP
    /// the loop — the density feedback a fixed window cannot capture. Default off so
    /// model_flare / the test suite keep the fixed-window deposition.
    bool  beam_thick_target       = false;
    float beam_E_cut_keV          = 20.0f;  ///< low-energy cutoff E_c [keV]
    float beam_delta              = 5.0f;  ///< injected power-law spectral index δ (>2)

    ///@}

    /** @name Ambient coronal footpoint heating (gentle evaporation) */
    ///@{
    /// Steady volumetric heating H(s) [W/m³] standing in for the (sub-grid) coronal
    /// heating mechanism that sustains a loop against conduction + radiation — the
    /// Rosner–Tucker–Vaiana (1978) static energy balance
    ///   d/ds(κ_e T^{5/2} dT/ds) + H(s) − n²Λ(T) = 0.
    /// make_loop_dat.py supplies the hydrostatic line; this term supplies H(s), the
    /// missing ingredient that makes a resolved chromosphere→TR→corona column a true
    /// steady state rather than a draining, prescribed-T atmosphere. Footpoint-
    /// anchored exponential (Aschwanden & Schrijver 2002 §3):
    ///   H(s) = coronal_heat_E0 · exp(−d(s)/s_H) · enhance(t),
    /// where d(s) is the field-aligned distance to the nearest footpoint
    /// (one-sided d = s − s₀ for open / half-loop; two-sided d = min over both ends
    /// for a full loop, set by coronal_heat_two_sided). The uniform-heating limit
    /// is recovered as s_H → ∞ (Martens 2010: T(s) is only weakly sensitive to the
    /// shape). UNLIKE the transient beam this is ALWAYS on for a steady run; gentle
    /// evaporation (Antiochos & Sturrock 1978) is driven by a slow time RAMP of the
    /// amplitude (coronal_heat_enhance > 1), kept sub-threshold so v ≪ c_s.
    ///
    /// Heating is deposited into the charged thermal pool exactly like the beam
    /// (apply_coronal_heating_stage): shared with neutrals by heat capacity in the
    /// single-T baseline, into the electron pool when enable_Te (electrons conduct).
    /// TWO-FLUID ONLY: the release timestep has no volumetric heating stage and
    /// ignores this flag.
    /// Default off so existing scenarios / the test suite are byte-for-byte unchanged.
    bool  enable_coronal_heating  = false;
    float coronal_heat_E0         = 0.0f;  ///< E_H0 [W/m³], footpoint heating amplitude
    float coronal_heat_sH         = 1.0e7f;  ///< heating scale length s_H [m]
    float coronal_heat_s0         = 0.0f;  ///< footpoint anchor s₀ [m] (arc length, cell-0 face = 0)
    bool  coronal_heat_two_sided  = false;  ///< true: decay from BOTH ends (full loop)
    /// Phase-3 ramp: multiply E_H0 by a smooth flat-top window that rises from 1 to
    /// coronal_heat_enhance over [t_on, t_on+ramp]. enhance = 1 (default) ⇒ no ramp
    /// (steady relaxation); enhance > 1 raises the heating to drive the upflow.
    float coronal_heat_t_on       = 0.0f;  ///< ramp onset time [s]
    float coronal_heat_ramp       = 1.0f;  ///< cosine ramp half-width [s]
    float coronal_heat_enhance    = 1.0f;  ///< peak amplitude multiplier after the ramp

    /// Positivity / vacuum floor (apply_flare_floor_stage + the rhs predictor and
    /// reconstructed-face floors). Clips ρ,p positive, caps |V|, and — where the
    /// ionization fraction f ≥ 0.99 — slaves the trace neutral onto the charged
    /// fluid (the exact single-fluid limit). Flare runs switch this on implicitly
    /// via enable_beam_heating; STEADY scenarios whose hot upper layers fully
    /// ionize (ρ_n→0 → 0/0 velocity decode and a runaway spectral radius → NaN,
    /// e.g. the C7 column once the TR ionizes) opt in via enable_vacuum_floor.
    /// Default false keeps the test suite and untouched scenarios byte-for-byte.
    bool  enable_vacuum_floor     = false;
    /// True when the positivity/vacuum floors should run this step.
    bool  floors_active() const { return enable_beam_heating || enable_vacuum_floor; }

    /// When true, the outer boundary uses a transparent (Neumann) outflow on the
    /// velocity instead of the slow Mach-0.05 cap, so a supersonic evaporation
    /// upflow can leave the domain. Set by flare scenarios. Default false keeps
    /// the C7 quiet-Sun Mach-capped outflow. NOTE: a Neumann outflow is ill-posed
    /// for INflow (a downflow at the top feeds mass in unphysically); use the
    /// Mach cap (below) for runs whose residual is a downflow.
    bool  outer_free_outflow      = false;

    /// Mach-number fraction for the capped outer outflow (model_c7_update_bc):
    /// V_ghost is clamped to ±outer_mach_cap·c_s(T_TR). Default 0.05 (the quiet-Sun
    /// value the C7 tests pin). A gentle-evaporation run raises it (e.g. 0.5) so the
    /// subsonic upflow exits without the ill-posed free-outflow inflow — the steady
    /// relaxation flow (≪ cap) is unaffected.
    float outer_mach_cap          = 0.05f;

    /// When true, the outer face is a reflecting symmetry plane rather than an
    /// outflow: the outer ghosts mirror the interior with V,U → −V,−U, so the net
    /// mass/momentum flux across the face is zero. This is the loop-apex condition
    /// for a CLOSED field line (one symmetric leg, footpoint→apex), where the
    /// evaporated plasma is confined and fills the loop. Open field lines leave
    /// this false and keep the (Neumann) coronal outflow. Set by pfss_ic from the
    /// [META] topology= key. Default false → unchanged outflow for all other
    /// scenarios. apply_open_bcs honors this on the outer block.
    bool  outer_reflecting        = false;

    /// Current simulation time [s], refreshed by the driver (chromo_main) each
    /// step before advance_Euler_state so time-dependent terms (beam window) can
    /// read it. Not used by steady scenarios.
    float sim_time                = 0.0f;

    /// Per-Grid reusable workspaces for the release solver, reused between steps so
    /// the hot loops allocate nothing and no function-static state is needed.
    /// Mutable because the stages that use them take a `const Grid&`.
    ///@{
    mutable MixtureConductionScratch mixture_conduction_scratch;
    mutable MixtureRhsScratch mixture_rhs_scratch;
    ///@}

    /// Diagnostic-only: when true, the NEXT mixture_rhs_explicit evaluation copies
    /// its final reconstructed face states, limiter values and numerical mass
    /// flux into face_flux_capture. Pure output — no solver stage ever reads the
    /// capture back, so the flag cannot change the numerical result. Default false.
    /// The driver (chromo_main, CHROMO_FACE_FLUX_DIAG) sets it per step.
    bool capture_face_flux = false;
    mutable MixtureFaceFluxCapture face_flux_capture;  ///< destination of that capture

    /// Diagnostic-only tally of the release swmf_godunov flux, accumulated over
    /// the run. Untouched, and left at zero, when that flux is not selected.
    mutable GodunovFluxStats godunov_stats;

    /// Diagnostic-only: when true, every mixture_apply_conduction solve copies
    /// its FINAL CONVERGED outer-face conduction quantities into
    /// outer_conduction_capture. Pure output; no stage reads it back. Default false.
    /// The driver (chromo_main, CHROMO_OUTER_COND_DIAG) sets it once per run.
    bool capture_outer_conduction = false;
    mutable OuterConductionCapture outer_conduction_capture;  ///< destination of that capture

    ///@}

    /** @name Cell-centered and face arrays (length ns) */
    ///@{
    Vec ds_i;  ///< cell width along the field line [m]; the canonical mesh array
    Vec B_i;    ///< magnetic field strength at the cell centre [T]
    Vec B_imh;  ///< magnetic field strength at the lower face i-1/2 [T]
    Vec B_iph;  ///< magnetic field strength at the upper face i+1/2 [T]
    /// Arc-length derivative of 1/B at the cell centre [T^-1 m^-1]. With
    /// A(s) B(s) = const this is what drives the flux-tube pressure-area source;
    /// it is identically zero for the uniform-B scenarios.
    Vec dinvB_ds_i;
    Vec phi_g_imh;  ///< gravitational potential at the lower face i-1/2 [J kg^-1]
    Vec phi_g_iph;  ///< gravitational potential at the upper face i+1/2 [J kg^-1]

    ///@}

    /** @name Release ghost buffers (length num_of_mixture_eq) */
    ///@{
    /// The two hydrodynamic ghost layers of the release solver, in the same
    /// (rho, rho u, E) rows as the interior state. The conductive physical-face
    /// temperature boundary is a SEPARATE datum (outer_conduction_temperature);
    /// ghost states and thermal boundary data are never conflated.
    ///@{
    Vec mix_outer_boundary0, mix_outer_boundary1;  ///< outer ghosts, nearest first
    Vec mix_inner_boundary0, mix_inner_boundary1;  ///< inner ghosts, nearest first
    ///@}

    ///@}

    /** @name Legacy two-fluid ghost buffers (length num_of_eq) */
    ///@{
    /// The two ghost layers on each side of the historical two-fluid solver, in the
    /// same seven `cons::` rows as the interior state, nearest layer first.
    ///@{
    Vec outer_boundary0_i, outer_boundary1_i;
    Vec inner_boundary0_i, inner_boundary1_i;
    ///@}

    ///@}

    /** @name Packed broadcasts (length n_state) */
    ///@{
    /// Cached versions of B_i, B_imh, B_iph, dinvB_ds_i, ds_i, and a scratch
    /// dt buffer broadcast across all num_of_eq equations. Recomputed by
    /// broadcast() whenever the cell-centered fields change.
    ///@{
    Vec B_state, B_state_imh, B_state_iph;  ///< B_i, B_imh, B_iph broadcast [T]
    Vec dinvB_ds_state;  ///< dinvB_ds_i broadcast [T^-1 m^-1]
    Vec ds_state;        ///< ds_i broadcast [m]
    Vec dt_state;        ///< scratch timestep buffer written by the integrators [s]
    ///@}

    ///@}

    /** @name Static-mesh metric caches */
    ///@{
    /// Support for a metric-aware, static non-uniform finite-volume mesh (static
    /// local refinement — NOT AMR: no runtime regridding). Every quantity here is
    /// derived from the canonical cell-width array `ds_i`, taking the cell CENTER
    /// at the geometric midpoint of its two faces:
    ///   * `s_face`   — face arc length [m] from the inner boundary (s_face[0]=0),
    ///                  length ns+1.
    ///   * `s_i`      — cell-center arc length [m], length ns.
    ///   * `ds_iph_i` — center-to-center distance to the i+1 neighbour [m]. At the
    ///                  outer boundary the ghost width mirrors the last cell, so
    ///                  ds_iph_i[ns-1] = ds_i[ns-1] (Neumann), matching the legacy
    ///                  0.5*(ds_i + ip1(ds_i,SLICE)) formula used in rhs/conduction.
    ///   * `ds_imh_i` — center-to-center distance to the i-1 neighbour [m], with the
    ///                  inner ghost width mirrored: ds_imh_i[0] = ds_i[0].
    /// On a UNIFORM mesh ds_iph_i = ds_imh_i = ds_i, `uniform_mesh` is true, and
    /// every spacing-sensitive operator keeps its legacy code path (byte-for-byte
    /// identical results). broadcast() sets `uniform_mesh` by scanning ds_i and
    /// rejects (throws) any non-positive / non-finite width.
    ///@{
    Vec  s_i;        ///< cell-centre arc length from the inner boundary [m], length ns
    Vec  s_face;     ///< face arc length, s_face[0] = 0 [m], length ns+1
    Vec  ds_iph_i;   ///< centre-to-centre distance to the i+1 neighbour [m]
    Vec  ds_imh_i;   ///< centre-to-centre distance to the i-1 neighbour [m]
    ///@}
    /// True when every ds_i is equal, in which case each spacing-sensitive operator
    /// takes its uniform-mesh code path and results are byte-for-byte identical to
    /// a build without non-uniform support. Set by broadcast().
    bool uniform_mesh = true;

    /// Packed (n_state) center-to-center distances, used by the non-uniform MUSCL
    /// reconstruction and conduction paths. Populated by broadcast() on every mesh
    /// (they equal ds_state on a uniform mesh).
    ///@{
    Vec ds_iph_state, ds_imh_state;
    ///@}
    ///@}

    /// Allocate arrays for ns cells, fill physical constants, zero the fields.
    void init(arma::uword ns_in, float CFL_in);

    /// Re-allocate the ns-sized arrays for a new cell count WITHOUT touching the
    /// physical constants, γ, CFL, or any runtime toggle. Used by the static-mesh
    /// scenarios: peek_ns() returns the coarse-equivalent count, then the IC builds
    /// the (possibly refined) face grid and resizes so the allocated count matches
    /// the generated mesh. No-op when ns_new == ns.
    void resize(arma::uword ns_new);

    /// Recompute the packed `_state` broadcasts and the static-mesh metric caches
    /// from ds_i, B_i, B_imh, B_iph, dinvB_ds_i. Call after editing any of those
    /// fields. Throws std::runtime_error on a non-positive / non-finite ds_i.
    ///
    /// Everything broadcast() derives depends ONLY on the static mesh (ds_i) and
    /// the magnetic geometry (B_i, B_imh, B_iph, dinvB_ds_i) — never on the ghost
    /// buffers or the conserved state. Several scenarios nevertheless call it from
    /// their per-step boundary refresh, where it re-derived ~50 n_state-sized
    /// temporaries per timestep for an unchanged mesh. broadcast() therefore
    /// fingerprints those five source arrays and rebuilds only when one of them
    /// actually changed (grid initialization, resize, mesh construction, magnetic
    /// geometry edits); an ordinary boundary update becomes an O(ns) compare.
    /// force_rebuild_metrics() bypasses the fingerprint.
    void broadcast();

    /// Unconditional rebuild of the static-mesh / packed-geometry caches.
    void force_rebuild_metrics();

    /// Bumped by every force_rebuild_metrics(). Downstream caches derived from
    /// the static mesh (e.g. MixtureRhsScratch's reconstruction weights) compare
    /// against it to know when they must be rebuilt.
    std::uint64_t metrics_generation() const { return metrics_generation_; }

private:
    /// Fingerprint copies of the five geometry arrays the caches were last derived
    /// from. broadcast() compares against these and rebuilds only on a real change,
    /// so a per-step boundary refresh costs an O(ns) comparison.
    ///@{
    Vec metrics_ds_i, metrics_B_i, metrics_B_imh, metrics_B_iph,
        metrics_dinvB_ds_i;
    ///@}
    bool metrics_valid = false;  ///< false until the caches have been built once
    std::uint64_t metrics_generation_ = 0;  ///< bumped by every force_rebuild_metrics()

    /// True when the fingerprint above still matches the live geometry arrays.
    bool static_metrics_current() const;
};

/** @} */

} // namespace chromosphere
