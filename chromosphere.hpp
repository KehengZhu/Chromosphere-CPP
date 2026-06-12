/*!
 * Chromosphere model — public API.
 *
 * All solver state lives in a Grid struct that the caller owns; functions
 * take `const Grid&` (or `Grid&` for those that mutate scratch buffers).
 * No global mutable state — multiple Grids can coexist.
 *
 * This header contains code derived from CoMFi (MIT License).
 */

#pragma once

#include <armadillo>

namespace chromosphere {

/// Armadillo float column vector — the workhorse type for cell-wise data.
typedef arma::Col<float> Vec;

// ============================================================================
// State indices
// ============================================================================

const arma::uword num_of_eq = 7;

// Conserved-variable indices (writeup eq 61).
//
// Three-temperature extension (docs/electron_temperature_plan.md): a 7th
// variable E_E carries the electron internal energy. To keep the change a
// clean, energy-conserving generalization, E_I retains its original meaning —
// the TOTAL charged-fluid energy (protons + electrons + bulk KE + gravity) —
// so the conservative MUSCL/Rusanov flux, gravity/area source, and sound speed
// are byte-for-byte unchanged. The proton temperature is then a *derived*
// quantity, T_i = (p_total − p_e)/(n_i k_B), and the electron temperature is
// T_e = p_e/(n_i k_B) with p_e = ⅔ E_E. This "total-energy + electron-energy"
// split is the standard two-temperature MHD formulation (AWSoM/BATSRUS;
// Sokolov 2021) and guarantees ENABLE_TE=0 reproduces the single-temperature
// baseline exactly (E_E never feeds back into E_I / momentum / ρ).
namespace cons {
    const arma::uword RHO_I = 0; // ion mass density       ρ_i
    const arma::uword RHO_N = 1; // neutral mass density   ρ_n
    const arma::uword MOM_I = 2; // ion momentum density   ρ_i V
    const arma::uword MOM_N = 3; // neutral momentum dens. ρ_n U
    const arma::uword E_I   = 4; // TOTAL charged energy   e_i (protons+electrons+KE+φ)
    const arma::uword E_N   = 5; // neutral total energy   e_n
    const arma::uword E_E   = 6; // electron internal energy ε_e = (3/2) p_e
}

// Primitive-variable indices. P_I is the TOTAL charged pressure (p_proton+p_e);
// P_E is the electron partial pressure. Proton pressure = P_I − P_E (derived).
namespace prim {
    const arma::uword RHO_I = 0;
    const arma::uword RHO_N = 1;
    const arma::uword V     = 2;
    const arma::uword U     = 3;
    const arma::uword P_I   = 4;  // total charged pressure p_proton + p_e
    const arma::uword P_N   = 5;
    const arma::uword P_E   = 6;  // electron partial pressure p_e
}

// Used by ip1/im1/...: SLICE means "treat input as one scalar field of length ns",
// CUBE means "treat input as the packed state of length ns*num_of_eq".
const arma::uword SLICE = 1;
const arma::uword CUBE  = num_of_eq;

// Outer/inner ghost-cell handling for ip2/im2: copy interior (Neumann) when true,
// otherwise pull from the ghost buffers in Grid.
const bool USE_NEUMANN_BC = false;

// ============================================================================
// Grid: owns all solver state.
// ============================================================================

struct Grid {
    // --- sizes -------------------------------------------------------------
    arma::uword ns      = 0;
    arma::uword n_state = 0;            // ns * num_of_eq
    float       CFL     = 0.25f;

    // --- physical constants (SI) ------------------------------------------
    float gamma_mono = 5.0f / 3.0f;
    float m_i        = 1.6726219e-27f;
    float m_n        = 1.6726219e-27f;
    float m_e        = 9.10938356e-31f;
    float g          = 0.27395e3f;
    float mu_0       = 4.0f * static_cast<float>(arma::datum::pi) * 1.0e-7f;
    float k_b        = 1.380649e-23f;
    float q_e        = 1.602176634e-19f;
    float chi_H_J    = 2.179872361e-18f;     // hydrogen ionization potential, 13.6 eV in J

    // --- runtime toggles ---------------------------------------------------
    // Enable Stage E (hydrogen ionization / recombination) in advance_Euler_state.
    // Default off so the existing test suite remains a clean regression baseline;
    // flip to true to activate the writeup §5.3 ionization stage.
    bool enable_ionization = false;

    // Sub-options for the Stage E ionization network (only active when
    // enable_ionization = true). Both default ON: the full network is
    //   ionization:    P_phot (photoionization) + S_CR (multilevel collisional)
    //                  + S_i (Voronov ground-state electron-impact)
    //   recombination: α_r (case-B radiative) + α_c = κ_c n_e (three-body)
    // No single channel dominates across the parameter range this code spans
    // (cool weakly-ionized chromosphere ~6 kK through the dense flare
    // condensation to the ~10⁷ K, fully-ionized evaporated plasma), so the
    // complete network is the correct default; each channel is the leading
    // term in some regime (see util/visualization/c7_ioniz_recomb_rates.png):
    //
    // enable_direct_collisional_ionization: S_i (Voronov 1997 ground-state
    //   electron-impact rate). Subdominant to the multilevel S_CR by 2–3 orders
    //   in the cool chromosphere, but it is the *physically robust* ionization
    //   channel above ~10⁵ K (transition region / flare), where the cool-plasma
    //   S_CR fit is an extrapolation (S_i overtakes S_CR near 1.1×10⁵ K).
    //
    // enable_threebody_recombination: α_c = κ_c n_e (Hinnov & Hirschberg 1962 /
    //   Stevefelt 1975). Its volumetric rate ∝ n_e³ makes it the dominant
    //   recombination sink wherever the gas is compressed — the dense base and,
    //   critically, the shock-compressed flare condensation, where dropping it
    //   destabilizes the run. Enabling it makes the Stage E solve cubic rather
    //   than quadratic. Negligible only in the tenuous upper chromosphere.
    bool enable_direct_collisional_ionization = true;
    bool enable_threebody_recombination       = true;

    // Enable optically-thick chromospheric radiative cooling (CL2012 recipe,
    // writeup §2.4) as Stage R in the operator-split integrator. Default off
    // so the existing test suite remains a clean regression baseline; flip
    // to true to activate.
    bool enable_radiative_cooling = false;

    // Single-fluid limit ("neutrals off"). When true, the ion-neutral drag
    // (Stage B) and T_i↔T_n equilibration (Stage C) rigidly slave the neutral
    // fluid to the ions each step: V→U→V_cm (mass-weighted bulk velocity) and
    // T_i→T_n→T_cm (heat-capacity-weighted bulk temperature). This is the exact
    // infinite-coupling limit of the two-fluid equations — neutrals still carry
    // mass and follow the ionization network, but have no independent dynamics,
    // i.e. a single-fluid partially-ionized treatment. Default false = full
    // two-fluid ("neutrals on"). Toggled at runtime by SINGLE_FLUID=1.
    bool single_fluid = false;

    // Separate electron temperature T_e ≠ T_i (docs/electron_temperature_plan.md).
    // The 7th conserved variable E_E (electron internal energy) is ALWAYS carried;
    // this flag only decides whether the electrons evolve independently:
    //   * enable_Te = false (default): the single-temperature baseline. Every
    //     stage runs the original code path on the combined charged pool
    //     (P_I = total charged pressure, "T_i" = T_charged = p_total/(2 n_i k_B)),
    //     and at the end of each step E_E is slaved to ½ the charged thermal
    //     energy so that T_e ≡ T_i is reported. E_E never feeds back, so this
    //     reproduces the pre-T_e physics to round-off (cf. SINGLE_FLUID).
    //   * enable_Te = true: the three-temperature model (T_e, T_i, T_n). The
    //     beam deposits into electrons, Spitzer κ_e conducts on T_e, the
    //     optically-thin radiation and the χ_H ionization cost are electron-pool
    //     terms, and a 3-way point-implicit ν_ei/ν_en/ν_in relaxation (Stage C)
    //     equilibrates the three temperatures. Decoupling (T_e ≠ T_i) appears at
    //     the beam-driven onset and the tenuous loop-top where ν_ei is slow
    //     (Bradshaw 2006; Manchester 2012). Toggled at runtime by ENABLE_TE=1.
    bool enable_Te = false;

    // Isotropic numerical thermal diffusivity χ_num [m²/s] added to the Stage D
    // conduction operator (Pandey et al. 2024 §3.3): an explicit ∂T/∂t = χ ∂²T/∂s²
    // term, realized as an additive face conductivity K_num = χ·C (C = heat
    // capacity per volume), that damps under-resolved TR / base gradients and
    // reflected boundary waves which Spitzer conduction alone leaves unstable.
    // Pandey's benchmark is χ ≳ 30×10⁸ m² s⁻¹ at Δz = 40 km; finer grids need
    // larger χ. Default 0 so the existing test suite stays a clean baseline;
    // scenarios scale it to Δs (model_c7 sets it in model_c7_ic).
    float numerical_diffusivity = 0.0f;

    // Upper-boundary coronal conductive heat flux q(T) [W/m²], imposed as a
    // Neumann (fixed-flux) outer BC on the charged-fluid conduction row (Stage D)
    // when impose_outer_heat_flux is true. This is the RTV (1978) downward
    // coronal conductive flux F_c = −κ₀T^{5/2}dT/ds re-radiated at the TR base —
    // the heating that holds the chromosphere top against Route-B ionization
    // (χ_H) cooling. A fixed *flux* (not a fixed hot temperature) is used because
    // on the coarse grid a Dirichlet T-wall would over-conduct via κ_e∝T^{5/2}
    // and over-pressurize the top cell. Positive = heating into the domain.
    // Default off so existing tests keep the Dirichlet conduction BC.
    bool  impose_outer_heat_flux = false;
    float outer_heat_flux        = 0.0f;

    // Transition-Region Adaptive Conduction (TRAC; Johnston et al. 2019, 2020).
    // When enabled, an adaptive cutoff temperature T_c is recomputed each step
    // (the highest temperature where the TR is under-resolved, L_R/L_T > 1/2,
    // bounded by [trac_T_chrom, 0.2 T_peak]); below T_c the parallel conductivity
    // κ_e is broadened by ε=(T_c/T)^{5/2} and the optically-thin radiative loss
    // (and any volumetric heating) divided by ε. This broadens the unresolved TR
    // so a coarse grid resolves it, while keeping κΛ and κQ — hence the TR-
    // integrated radiation/heating and the coronal density response — invariant
    // (Johnston 2020 §2.2), which is the Eq.-27 jump-condition balance enforced
    // analytically. Default off so existing tests are a clean baseline.
    bool  enable_trac   = false;
    float trac_T_chrom  = 2.0e4f;   // TRAC-region base temperature T_b [K]
    float trac_cutoff_T = 2.0e4f;   // adaptive cutoff T_c [K], recomputed each step

    // Effective hydrogen photoionization rate from the ground state [s^-1],
    // modeling the Lyα-excitation + Balmer-continuum two-step channel
    // (Carlsson & Stein 2002 establish this as the dominant chromospheric
    // ionization path; the Lyman continuum itself contributes negligibly).
    // C&S 2002 report a chromospheric ionization/recombination timescale of
    // 10^3–10^5 s; we adopt the geometric mean ~10^4 s as a single-rate
    // closure suitable for the 2-state (ground+continuum) reduction used
    // here. See Leenaarts et al. 2007 for the "fixed radiative rates"
    // implementation paradigm in non-equilibrium chromospheric MHD codes.
    //
    // The scalar `photoionization_rate` is the uniform fallback used when
    // the cell-centered `photoionization_rate_i` array is empty (i.e. no
    // height-dependent calibration has been performed by the scenario IC).
    // Scenarios that want a height-dependent profile (e.g. model_c7_ic,
    // which calibrates against C7's tabulated f to make C7 a fixed point
    // of the Stage E quadratic) populate `photoionization_rate_i` of
    // length `ns`; the Stage E driver picks that up automatically via
    // `physics.hpp::photoionization_rate_P(grid)`.
    float photoionization_rate = 1.0e-4f;
    Vec   photoionization_rate_i;

    // --- flare beam heating (explosive chromospheric evaporation) ----------
    // Transient nonthermal-electron-beam volumetric heating that drives
    // chromospheric evaporation (Fisher, Canfield & McClymont 1985). The beam,
    // injected from the corona at the top, deposits an energy flux beam_flux
    // [W/m²] in the upper chromosphere over a layer beam_deposition_height [km]
    // below the outer boundary, weighted ∝ local total density (the densest
    // reachable layer absorbs the most — thick-target stopping; the true
    // ~10²⁰ cm⁻² stopping column is sub-grid on the coarse C7 mesh). The flux is
    // gated in time by a flat-top window [beam_t_on, beam_t_on+beam_duration]
    // with cosine ramps of width beam_ramp. Explosive evaporation needs
    // beam_flux above the Fisher threshold F_crit ≈ 7×10⁶ W m⁻² AND the
    // radiative sink (Stage R) active, since the threshold is the competition
    // between deposited heating and the upper-chromosphere radiative loss.
    // Heating is added to the electron/ion thermal pool (p_i). Default off so
    // existing scenarios/tests are unaffected.
    bool  enable_beam_heating     = false;
    float beam_flux               = 0.0f;    // F_e [W/m²]
    float beam_t_on               = 0.0f;    // beam onset time [s]
    float beam_duration           = 0.0f;    // flat-top duration τ [s]
    float beam_ramp               = 1.0f;    // cosine ramp half-width [s]
    // Deposition height window [km] (absolute, C7 base = 1003 km). The beam
    // stops in the upper chromosphere BELOW the TR base (~2150 km); depositing
    // there (not in the TR cells at the open top boundary) is both physical
    // (thick-target stopping in the dense chromosphere) and keeps the explosive
    // overpressure interior so the evaporated upflow develops through the TR.
    float beam_h_lo_km            = 1600.0f;
    float beam_h_hi_km            = 2000.0f;

    // Self-consistent thick-target beam (Emslie 1978). When true, the fixed
    // [h_lo,h_hi] deposition window is ignored: the beam is injected at the loop
    // APEX (the maximum of φ_g — loop top for closed/full, coronal top for open),
    // streams down each leg, and loses energy collisionally. The volumetric
    // heating at column depth N (from the apex) is Q ∝ n_tot·(1+N/N_c)^{-δ/2},
    // normalized so ∫Q ds = beam_flux (thick target: all energy absorbed). N_c is
    // the stopping column of cutoff-energy electrons (N_c ∝ E_c²). As the loop
    // fills with evaporated plasma the column rises, so the stopping depth moves UP
    // the loop — the density feedback a fixed window cannot capture. Default off so
    // model_flare / the test suite keep the fixed-window deposition.
    bool  beam_thick_target       = false;
    float beam_E_cut_keV          = 20.0f;   // low-energy cutoff E_c [keV]
    float beam_delta              = 5.0f;    // injected power-law spectral index δ (>2)

    // When true, the outer boundary uses a transparent (Neumann) outflow on the
    // velocity instead of the slow Mach-0.05 cap, so a supersonic evaporation
    // upflow can leave the domain. Set by flare scenarios. Default false keeps
    // the C7 quiet-Sun Mach-capped outflow.
    bool  outer_free_outflow      = false;

    // When true, the outer face is a reflecting symmetry plane rather than an
    // outflow: the outer ghosts mirror the interior with V,U → −V,−U, so the net
    // mass/momentum flux across the face is zero. This is the loop-apex condition
    // for a CLOSED field line (one symmetric leg, footpoint→apex), where the
    // evaporated plasma is confined and fills the loop. Open field lines leave
    // this false and keep the (Neumann) coronal outflow. Set by pfss_ic from the
    // [META] topology= key. Default false → unchanged outflow for all other
    // scenarios. apply_open_bcs honors this on the outer block.
    bool  outer_reflecting        = false;

    // Current simulation time [s], refreshed by the driver (chromo_main) each
    // step before advance_Euler_state so time-dependent terms (beam window) can
    // read it. Not used by steady scenarios.
    float sim_time                = 0.0f;

    // --- cell-centered & face arrays (length ns) --------------------------
    Vec ds_i;
    Vec B_i, B_imh, B_iph;
    Vec dinvB_ds_i;
    Vec phi_g_imh, phi_g_iph;

    // --- ghost buffers (length num_of_eq) ---------------------------------
    Vec outer_boundary0_i, outer_boundary1_i;
    Vec inner_boundary0_i, inner_boundary1_i;

    // --- packed broadcasts of the above (length n_state) ------------------
    // Cached versions of B_i, B_imh, B_iph, dinvB_ds_i, ds_i, and a scratch
    // dt buffer broadcast across all num_of_eq equations. Recomputed by
    // broadcast() whenever the cell-centered fields change.
    Vec B_state, B_state_imh, B_state_iph;
    Vec dinvB_ds_state, ds_state, dt_state;

    /// Allocate arrays for ns cells, fill physical constants, zero the fields.
    void init(arma::uword ns_in, float CFL_in);

    /// Recompute the packed `_state` broadcasts from ds_i, B_i, B_imh, B_iph,
    /// dinvB_ds_i. Call after editing any of those fields.
    void broadcast();
};

// ============================================================================
// Packed-state helpers
// ============================================================================

/// Extract one scalar field (length ns) from a packed state.
Vec get_scalar(const Grid& grid, const Vec& xn_state, arma::uword index);

/// Insert a scalar field (length ns) into the given slot of a packed state
/// (length n_state, other slots zero).
Vec scalar_to(const Grid& grid, const Vec& xn_i, arma::uword index);

/// Index shifts. `nk = CUBE` treats the input as packed state; `nk = SLICE`
/// treats it as one scalar field of length ns.
Vec ip1(const Grid& grid, const Vec& xn, arma::uword nk = CUBE);
Vec im1(const Grid& grid, const Vec& xn, arma::uword nk = CUBE);
Vec ip2(const Grid& grid, const Vec& xn, arma::uword nk = CUBE);
Vec im2(const Grid& grid, const Vec& xn, arma::uword nk = CUBE);

// ============================================================================
// Equation of state
// ============================================================================

Vec cons2prim(const Grid& grid, const Vec& cons_state);
Vec prim2cons(const Grid& grid, const Vec& prim_state);

// ============================================================================
// Numerics: flux, source, spectral radius, MUSCL limiter
// ============================================================================

/// Minmod limiter (writeup eq 15): φ(r) = max(0, min(1, r)).
Vec flux_lim(const Vec& r);

/// Cell-centered flux F(U).
Vec cal_flux_state(const Grid& grid, const Vec& xn_state);

/// Per-cell spectral radius ρ(∂F/∂U) broadcast across all num_of_eq equations.
Vec cal_spectral_radius_state(const Grid& grid, const Vec& xn_state);

/// Explicit source term (pressure-area + gravity).
Vec cal_source_state(const Grid& grid, const Vec& xn_state);

// ============================================================================
// RHS of the conservation laws (writeup §4)
// ============================================================================

/// Explicit RHS: TVD-MUSCL + Rusanov flux differencing + cal_source_state.
Vec rhs_explicit_state(const Grid& grid, const Vec& xn_state);

/// Implicit RHS: ion-neutral drag, collisional + frictional heating, and
/// field-aligned heat conduction (writeup §4.3–4.4).
Vec rhs_implicit_state(const Grid& grid, const Vec& xn_state);

// ============================================================================
// Time stepping
// ============================================================================

/// Per-cell maximum signal speed (extracted from the spectral radius broadcast).
Vec cal_max_v_i(const Grid& grid, const Vec& xn_state);

/// Uniform CFL-limited timestep.
Vec cal_dt_i(const Grid& grid, const Vec& xn_state);

/// Semi-implicit (backward-Euler) integrator. Explicit MUSCL+Rusanov step
/// for R_E, then point-implicit drag + frictional heating, point-implicit
/// ion–neutral temperature equilibration, a tridiagonal heat-conduction
/// solve per species, and (if grid.enable_ionization) the point-implicit
/// ionization Stage E (writeup §3.7, §5.3). Mutates grid.dt_state.
Vec advance_Euler_state(Grid& grid, const Vec& xn_state, const Vec& dt_i);

/// Stage E: backward-Euler ionization / recombination on a single cell, scalar
/// quadratic solve in ionization fraction f = ρ_i / (ρ_i + ρ_n). Operates on
/// primitive state in place; intended for the operator-split integrator but
/// exposed for direct testing. (writeup §5.3.)
void apply_ionization_stage(const Grid& grid, Vec& prim_state, float dt);

/// Stage R (writeup §2.4): backward-Euler optically-thick chromospheric
/// radiative cooling on the folded ion+electron thermal pressure. Carlsson
/// & Leenaarts 2012 recipe summed over H I + Ca II + Mg II; tables in
/// physics.hpp::cl2012. Activated by Grid::enable_radiative_cooling.
void apply_radiative_cooling_stage(const Grid& grid, Vec& prim_state, float dt);

/// Flare beam-heating stage: deposits grid.beam_flux into the upper-chromosphere
/// electron/ion thermal pool (p_i) over a finite layer, gated by the temporal
/// window around grid.sim_time, to drive chromospheric evaporation (Fisher et
/// al. 1985). No-op unless grid.enable_beam_heating. Mutates p_i in place.
void apply_beam_heating_stage(const Grid& grid, Vec& prim_state, float dt);

/// TRAC adaptive cutoff temperature T_c (Johnston et al. 2020 Eq. 8): the
/// maximum temperature among grid cells where the charged-fluid TR is
/// under-resolved (L_R/L_T > 1/2, with L_T = T/|dT/ds|, L_R = Δs), clamped to
/// [grid.trac_T_chrom, 0.2 T_peak]. Used by the conduction and cooling stages to
/// broaden the unresolved TR. Exposed for testing.
float compute_trac_cutoff_T(const Grid& grid, const Vec& prim_state);

/// Pure-explicit forward-Euler step — only the MUSCL+Rusanov R_E predictor of
/// advance_Euler_state, with the implicit drag / temperature / conduction
/// stages skipped (equivalent to R_I ≡ 0). Mutates grid.dt_state.
Vec advance_Euler_explicit_state(Grid& grid, const Vec& xn_state, const Vec& dt_i);

/// Explicit RK4. Mutates grid.dt_state.
Vec advance_RK4(Grid& grid, const Vec& xn_state, const Vec& dt_i);

// ============================================================================
// Debug
// ============================================================================

void print_xn(const Grid& grid, const Vec& xn);

} // namespace chromosphere
