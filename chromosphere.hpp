/*!
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

/// Armadillo float column vector — the workhorse type for cell-wise data.
typedef arma::Col<float> Vec;

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
    const arma::uword RHO    = 0;
    const arma::uword MOM    = 1;
    const arma::uword ENERGY = 2;
}
const arma::uword num_of_mixture_eq = 3;

// ============================================================================
// Historical two-fluid state indices — NON-RELEASE research solver
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

/// Per-Grid reusable storage for the release nonlinear conduction solve. This is
/// scratch only, never a thermodynamic cache; keeping it on the owning Grid
/// preserves reentrancy and avoids function-static state.
struct MixtureConductionScratch {
    std::vector<double> rho, e_old, temperature, target;
    // `e_at_T` is e_int(rho, temperature) of the CURRENT Newton iterate. The
    // ionization fraction, the heat capacity, the conductivity inputs and this
    // residual energy all come from one fused Saha evaluation per cell per pass
    // (they used to be three independent Saha solves).
    std::vector<double> conductivity, capacity, n_e, n_hi, e_at_T;
    std::vector<double> g_left, g_right;
    std::vector<double> a, b, c, rhs, delta;

    void resize(std::size_t n);
};

struct Grid;

/// Equilibrium thermodynamics decoded from one immutable packed release state
/// U = (rho, rho u, E). Physical cells and two explicitly decoded ghost layers
/// share the same state-carried gravitational-potential convention.
struct MixtureField {
    const Grid* source_grid = nullptr;
    const Vec* source_state = nullptr;
    const float* source_memory = nullptr;
    arma::uword source_elements = 0;
    std::uint64_t state_generation = 0;
    std::array<float, 4*num_of_mixture_eq> boundary_signature{};
    std::vector<MixtureThermo> cells;
    /// (ns+4) x 3 reconstruction variables (log rho, u, log p | log T).
    arma::Col<double> extended_primitive;
    /// Decoded temperature of every extended cell (ns+4, ghosts included), in the
    /// same indexing as extended_primitive. It is NOT a reconstruction variable:
    /// it exists only so the pressure-based face builder can seed its (rho,p)->T
    /// Newton iteration with the temperature of the cell each one-sided face
    /// state was extrapolated from.
    arma::Col<double> extended_temperature;

    void require_matches(const Grid& grid, const Vec& state) const;
};

/// One side of every interface, as the release flux needs it. All three-row
/// arrays use the (ns, num_of_mixture_eq) packing of the conserved state.
struct MixtureFaceArrays {
    arma::Col<double> conserved;    // ns*3, (rho, rho u, E)
    arma::Col<double> flux;         // ns*3, (rho u, rho u^2 + p, (E+p) u)
    arma::Col<double> velocity;     // ns
    arma::Col<double> pressure;     // ns
    arma::Col<double> sound_speed;  // ns
    arma::Col<double> dp_deint;     // ns, (dp/de_int)_rho

    void resize(arma::uword ns);
};

/// Reusable storage for the release MUSCL predictor/reconstruction. The numeric
/// arrays are deliberately untyped scratch; thermodynamic validity remains in
/// the state-bound MixtureField objects.
struct MixtureRhsScratch {
    MixtureField predicted;
    std::array<arma::Col<double>, 24> work;
    std::array<MixtureFaceArrays, 4> faces;
    arma::Col<double> flux_iph, flux_imh, source;
    Vec predicted_state, rhs;

    // Non-uniform MUSCL reconstruction weights and slope-ratio metrics. They are
    // functions of the STATIC mesh alone (ds_i via ds_iph_i/ds_imh_i), so they are
    // built once per mesh generation instead of once per timestep.
    //   W1/W3/W4       half-width weights of the owning cell
    //   metric_r/_ip1/_im1  center-to-center gradient-ratio corrections
    arma::Col<double> W1, W3, W4, metric_r, metric_r_ip1, metric_r_im1;
    std::uint64_t weights_generation = 0;
    bool weights_valid = false;
};

/// READ-ONLY capture of the release reconstruction and numerical mass flux for
/// ONE explicit RHS evaluation. Filled by mixture_rhs_explicit only when
/// Grid::capture_face_flux is true, and never read back by any solver stage — so
/// enabling it cannot change a numerical result. Index i refers to the UPPER face
/// i+1/2 of cell i, which is exactly the face the continuity row differences
/// (flux_iph); the cell-centred slots are the reconstruction INPUTS w of the same
/// evaluation, so a capture record is self-contained.
struct MixtureFaceFluxCapture {
    // cell-centred (length ns): the decoded state the reconstruction saw
    std::vector<double> rho_cell, v_cell, T_cell;
    // one-sided reconstructed face states at i+1/2 (L = from cell i, R = cell i+1)
    std::vector<double> rho_L, rho_R, v_L, v_R, T_L, T_R, cs_L, cs_R;
    // Authoritative total pressure of the same corrector face states, rebuilt
    // through the production face builder on capture steps only. p_L vs p_R is
    // the reconstruction-induced pressure mismatch the (log rho, V, log p)
    // reconstruction targets; p_cell is the decoded cell-centre value.
    std::vector<double> p_cell, p_L, p_R;
    // Acoustic spectral radius a = max(|V_L|+c_L, |V_R|+c_R), retained as a
    // diagnostic/fallback value in Roe-local mode, plus the central/dissipative split.
    std::vector<double> a_face, f_central, f_diff, f_total;
    // limiter inputs/outputs actually used for this face, per reconstruction slot:
    // r/phi_plus build the L state from cell i, r_ip1/phi_minus the R state from i+1
    std::vector<double> r_rho, phi_plus_rho, r_ip1_rho, phi_minus_rho;
    std::vector<double> r_v, phi_plus_v, r_T, phi_plus_T;
    bool valid = false;

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
    double T_top = 0.0, T_wall = 0.0;   // [K] top-cell centre / imposed thermal datum
    double kappa_face = 0.0;            // [W m^-1 K^-1] physical conductivity
    double ds_face = 0.0;               // [m] top-cell centre to the thermal datum
    double area_ratio = 0.0;            // B_i(ns-1)/B_iph(ns-1)
    double q_face = 0.0;                // [W m^-2], into the top cell
    bool   imposed_neumann = false;     // impose_outer_heat_flux replaced the above
    bool   valid = false;
};

// ============================================================================
// Grid: owns all solver state.
// ============================================================================

struct Grid {
    // --- sizes -------------------------------------------------------------
    arma::uword ns              = 0;
    arma::uword n_mixture_state = 0;    // ns * num_of_mixture_eq  (release)
    arma::uword n_state         = 0;    // ns * num_of_eq          (legacy two-fluid)
    float       CFL             = 0.25f;

    // --- physical constants (SI) ------------------------------------------
    // Fixed adiabatic index of the LEGACY two-fluid solver only. The release
    // solver derives every thermodynamic index from the Saha/Gamma1 closure and
    // never reads gamma_mono.
    float gamma_mono = 5.0f / 3.0f;
    // CRASH Gamma1 table. A non-empty table selects the RELEASE single-fluid
    // equilibrium-mixture solver (mixture.hpp); an empty one selects the legacy
    // fixed-gamma two-fluid solver.
    EosGammaTable eos_gamma_table;
    bool eos_gamma_debug_clamp = false;
    // Per-cell temperature hint for the EOS inversion — the last T decoded in this
    // cell by ANY stage. Purely an accelerator: it seeds Newton, never changes the
    // converged root or the residual tolerance, and a stale or absent hint only
    // costs the ordinary safeguarded solve. Measured on the h1600 ns=1000 column:
    // a seeded call converges in 1.13 iterations with ZERO bisections, an unseeded
    // one takes 8.68 iterations and 4.05 bisections. Mutable so the const-Grid
    // decode paths (flux, conduction) can refresh it.
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
    // Adiabatic-index factors derived from gamma_mono so the equation of state
    // and every energy↔pressure / heat-capacity conversion tracks a single γ.
    // The thermal energy density of an ideal gas is ε = p/(γ−1), so the pressure
    // is p = (γ−1)[E − ½ρv² − ρφ] and the heat capacity per volume is C = nk/(γ−1).
    // At the γ=5/3 default these reduce to the historical literals: gm1()=2/3,
    // inv_gm1()=3/2, half_gm1()=1/3. NOTE: these are the *thermodynamic* index
    // only — Spitzer κ∝T^{5/2} and the ion–neutral collisional energy-exchange
    // coefficient (3 k_B α/(m_i+m_n)) are kinetic-theory factors, intentionally
    // NOT tied to gamma_mono.
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
    float chi_H_J    = static_cast<float>(eos_constants::chi_h); // 13.6 eV

    // --- runtime toggles ---------------------------------------------------
    // LEGACY TWO-FLUID ONLY. Enable Stage E (hydrogen ionization / recombination)
    // in advance_Euler_state. The release solver has no finite-rate ionization: its
    // ionization state is the Saha equilibrium value of the EOS closure.
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
    // term in some regime (see visualization/c7_ioniz_recomb_rates.png):
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

    // TWO-FLUID ONLY. Enable optically-thick chromospheric radiative cooling
    // (CL2012 recipe, writeup §2.4) as Stage R in the operator-split two-fluid
    // integrator. The single-fluid release timestep has no radiative stage and
    // ignores this flag entirely. Default off so the existing test suite remains
    // a clean regression baseline; flip to true to activate.
    bool enable_radiative_cooling = false;

    // LEGACY TWO-FLUID ONLY. Single-fluid limit ("neutrals off"). When true, the ion-neutral drag
    // (Stage B) and T_i↔T_n equilibration (Stage C) rigidly slave the neutral
    // fluid to the ions each step: V→U→V_cm (mass-weighted bulk velocity) and
    // T_i→T_n→T_cm (heat-capacity-weighted bulk temperature). This is the exact
    // infinite-coupling limit of the two-fluid equations — neutrals still carry
    // mass and follow the ionization network, but have no independent dynamics,
    // i.e. a single-fluid partially-ionized treatment. Default false = full
    // two-fluid ("neutrals on"). Toggled at runtime by SINGLE_FLUID=1.
    bool single_fluid = false;

    // LEGACY TWO-FLUID ONLY. Separate electron temperature T_e ≠ T_i
    // (docs/electron_temperature_plan.md). The release solver is a common-temperature
    // equilibrium mixture and never evolves a separate electron energy.
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

    // LEGACY TWO-FLUID ONLY. Well-balanced gravity in the explicit MUSCL
    // reconstruction. The release reconstruction decodes each cell with its own
    // phi_g, so it has no shifted-potential bias to correct.
    // The TVD-MUSCL
    // step limits PRIMITIVE variables, obtaining the neighbour pressures by
    // cons2prim of the spatially-shifted conserved states (ip1/im1/ip2/im2 in
    // rhs_explicit_state). But cons2prim subtracts the LOCAL-index φ_g, while the
    // shifted cell's total energy carries that cell's ρφ_g — so each neighbour
    // pressure is biased by ⅔ρ·Δφ_g = ⅔ρg·Δs, which turns the hydrostatic slope
    // dp/ds = −ρg into −⅓ρg (the reconstructed pressure gradient is only ⅓ of the
    // true one). A hydrostatic atmosphere then feels a spurious, resolution-
    // independent (γ−1)g ≈ ⅔g downward force and cannot stay static. (Stage-D
    // conduction already shifts φ_g correctly — rhs_implicit_state — so only the
    // explicit reconstruction is affected.) When true, rhs_explicit_state corrects
    // the shifted-state pressures with the shifted φ_g, so a hydrostatic profile is
    // a discrete fixed point (V≈0) and the isentropic relaxation
    // (docs/gentle_evaporation_downflow.md, model_column) holds. Default false
    // keeps every existing scenario / the 5950-test baseline byte-for-byte; set by
    // model_column. (Candidate to enable globally — it is the physically
    // correct reconstruction and removes the spurious chromospheric downflow.)
    bool well_balanced = false;

    // LEGACY TWO-FLUID ONLY. Log-space MUSCL reconstruction in rhs_explicit_state.
    // The release reconstruction is always log-space by construction.
    // In a gravitationally
    // stratified atmosphere ρ and p are EXPONENTIAL in height, so a linear TVD slope
    // misrepresents them and minmod clips the steep gradient to first order, injecting
    // Rusanov numerical diffusion that is worst at the dense lower boundary. When true,
    // the reconstruction limits log ρ / log p (the strictly-positive, exponentially-
    // stratified primitive slots {RHO_I, RHO_N, P_I, P_N, P_E}) instead — an isothermal
    // hydrostatic column is then piecewise-LINEAR, the minmod ratio r→1 passes the full
    // slope, and the reconstructed face jump (hence the diffusion) collapses. The SIGNED
    // velocities {V, U} stay linear (log undefined; no large dynamic range). Faces are
    // exp'd back before prim2cons so the scheme stays conservative; applied AFTER the
    // well-balanced φ_g pressure correction (which must see linear p). Default false
    // keeps every existing scenario / the test baseline byte-for-byte; always on in
    // model_column. Composes with well_balanced.
    bool log_reconstruct = false;

    // Slope limiter for the explicit MUSCL reconstruction (rhs_explicit_state).
    // Default false ⇒ the symmetric minmod limiter φ(r)=max(0,min(1,r)) used
    // everywhere (writeup eq 15; slimL=slimR so one flux_lim(r) feeds both faces
    // of a cell). True ⇒ the MC3 / Koren limiter (BATSRUS ModFaceValue 'mc3'), the
    // ASYMMETRIC third-order (κ=1/3) monotonized-central limiter: less diffusive
    // than minmod, so it clips the under-resolved TR gradient far less (the source
    // of the O(Δs) base-drainage truncation in the iso runs). Because it is
    // asymmetric, the '+' and '−' reconstruction terms take DIFFERENT limited
    // slopes — see flux_lim_mc3_plus / flux_lim_mc3_minus. limiter_beta is the
    // Koren β (2 ⇒ classic Koren; 1 ⇒ minmod-like). Default off keeps every
    // existing scenario / the 5958-test baseline byte-for-byte; always on in
    // model_column (β=2). Composes with log_reconstruct and well_balanced (it only
    // changes how the same differences are limited).
    bool  mc3_limiter = false;
    float limiter_beta = 2.0f;

    // RELEASE numerical flux. True (the model_column release setting) computes the
    // corrector dissipation from the 3x3 mixture Roe characteristic decomposition
    // of U = (rho, rho u, E); false is the face-local Rusanov/LLF reference solver,
    // which is also the automatic per-face fallback whenever a Roe average is not
    // admissible. Rusanov is kept for regression and controlled comparison only —
    // its acoustic-scale dissipation -1/2 a Delta U is far too large for this
    // very-low-Mach evaporation problem and leaves a persistent TR velocity ripple.
    // model_column_ic sets it; the struct default keeps the legacy solver untouched.
    bool roe_characteristic_flux = false;

    // RELEASE thermal reconstruction variable. True (the model_column release
    // setting) limits
    //   (log rho, u, log p)
    // and recovers the face temperature by inverting the same authoritative
    // closure (equilibrium_temperature_from_density_pressure). False limits
    //   (log rho, u, log T)
    // and obtains the face pressure afterwards from the nonlinear Saha mapping
    // p(rho,T) — so two independently limited variables set one mechanical
    // quantity, and a mechanically smooth (constant-p) state is NOT reproduced
    // across the partial-ionization transition; it is retained as a reference
    // configuration only. Both keep density and pressure positive by construction
    // and feed the identical equilibrium face builder, flux and source stages, so
    // this flag changes ONLY which pair of variables is authoritative at a face.
    bool pressure_reconstruct = false;

    // Equilibrium-reference ("δ-form") well-balancing. A general NON-isothermal
    // hydrostatic atmosphere (e.g. Model C7, with a temperature minimum) is only
    // reproduced to O(Δs) by the reconstruction, so the explicit MUSCL step applies
    // a small spurious force even at rest ⇒ a residual O(Δs) base drainage
    // (docs/boundary_conditions_plan.md). This device removes it for ANY
    // stratification: capture the scheme's explicit-RHS residual at a frozen
    // reference equilibrium once, R_eq = RHS(eq_state), and subtract it every step.
    // Then RHS(eq_state) − R_eq ≡ 0, so eq_state is an EXACT discrete steady state
    // (u = 0 held to round-off, resolution-independent), while real deviations from
    // equilibrium (conduction-driven evaporation) evolve normally minus that fixed
    // O(Δs) correction. eq_state is the (u = 0, HSE) IC set by the scenario;
    // eq_residual is computed lazily on the first step with subtraction disabled via
    // the empty() guard. Always on in model_column; both solvers support it, each in
    // its own conserved-state width.
    bool eq_wb = false;
    Vec  eq_state;      // frozen reference equilibrium (conserved), set by scenario IC
    Vec  eq_residual;   // cached explicit-RHS residual at eq_state (empty until computed)

    // model_c7 "New explanation" upper BC (docs/gentle_evaporation_downflow.md):
    // when set, model_c7_update_bc imposes at the TR-base top face a hydrostatic
    // ghost pressure (centered HSE ⇒ no spurious boundary downflow), a FIXED
    // temperature jump T_ghost = a·T_ref, b·a·T_ref (the downward conductive flux
    // q(T) then enters via the Stage-D Dirichlet ghost-T, not the imposed Neumann
    // flux), and an EOS ghost density at that (p, T) — replacing the Mach-capped
    // RTV reservoir. Set true only by the bare `model_c7` scenario (which also
    // turns on well_balanced and turns off impose_outer_heat_flux); default false
    // keeps model_flare / analytic_canopy / model_gentle — which share
    // model_c7_ic / model_c7_update_bc — and the 5950-test baseline unchanged.
    bool c7_tr_jump_bc = false;

    // BOTH solvers. Enable field-aligned heat conduction: the implicit physical
    // conduction stage of the release timestep, and Stage D of the two-fluid
    // advance_Euler_state. Default ON so the existing scenarios / test suite are
    // byte-for-byte unchanged. A relaxation phase
    // (docs/gentle_evaporation_downflow.md, model_column Stage 1) sets this false
    // to switch OFF all heat flux, so the relaxation holds its profile with
    // V ≈ 0; Stage 2 flips it back on, with the imposed top temperature supplying
    // the downward conductive flux.
    bool enable_conduction = true;

    // TWO-FLUID ONLY. Mesh-independent isotropic numerical-diffusion coefficient
    // C_num [m/s]. The two-fluid Stage D constructs the actual diffusivity
    // independently at each face as chi_num,f = C_num*Delta_s_f, using the real
    // centre-to-centre face distance, then adds K_num,f = chi_num,f*C_V,f to the
    // physical face conductivity, so it is also disabled when conduction is off.
    // The release conduction operator is structurally PHYSICAL-only: it has no
    // such term and never reads this field. Default 0.
    float numerical_diffusivity_per_length = 0.0f;

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

    // Conduction-only OUTER temperature override. Normally the outer thermal
    // datum is whatever temperature the HYDRO outer ghost happens to carry, so one
    // ghost state serves two roles at once: the hydro Riemann/reconstruction state
    // AND the thermal wall that drives conduction. When
    // outer_conduction_temperature_override is true the conduction rows instead
    // use the explicit outer_conduction_temperature [K], letting a scenario relax
    // the hydro ghost temperature (e.g. zero-gradient extrapolation from the live
    // top cell) while the conductive driving stays the unchanged fixed hot wall.
    //
    // IMPORTANT — this also relocates the datum. The imposed temperature is a
    // PHYSICAL-BOUNDARY-FACE value, not a ghost-centre value, so the release
    // conduction operator uses the half-cell distance d = 1/2 ds_top and evaluates
    // kappa at the face state reconstructed from the imposed external pressure and
    // that face temperature. Without the override the datum stays at the ghost
    // CENTRE, distance ds_iph_i(ns-1), with the ghost's own n_e / n_HI. No effect
    // when impose_outer_heat_flux is true (that path has no outer Dirichlet datum).
    // Default false ⇒ byte-identical to the pre-existing behaviour.
    bool  outer_conduction_temperature_override = false;
    float outer_conduction_temperature          = 0.0f;

    // Inner-boundary conduction BC. Default false ⇒ Dirichlet: the innermost
    // conduction face couples cell 0 to the inner ghost temperature (the
    // photospheric reservoir's fixed T — a heat sink at the base). True ⇒ Neumann
    // (zero conductive flux through the base, dT/ds = 0): the inner face coupling
    // is dropped from the Stage-D tridiagonal so the base is insulating and the
    // photosphere heats freely. Set by model_column from ISO_INNER_T_NEUMANN.
    bool  inner_conduction_neumann = false;

    // Time-ramp on the imposed coronal conductive flux q(T) — the gentle
    // conduction-driven evaporation driver (docs/gentle_evaporation_plan.md v2).
    // The scenario stores the steady quiet value in outer_heat_flux_base; the
    // boundary update sets outer_heat_flux = base · enhance(sim_time), where
    // enhance(t) rises with a cosine ramp from 1 to outer_heat_flux_enhance over
    // [t_on, t_on+ramp]. enhance = 1 (default) ⇒ no ramp (steady q, unchanged
    // behavior for model_c7 / model_flare). Raising it above 1 increases the
    // downward coronal conductive flux into the TR/upper chromosphere → the
    // gentle upflow (Antiochos & Sturrock 1978). Set by model_gentle.
    float outer_heat_flux_base    = 0.0f;
    float outer_heat_flux_t_on    = 0.0f;
    float outer_heat_flux_ramp    = 1.0f;
    float outer_heat_flux_enhance = 1.0f;

    // Corona-as-boundary domain top [km] for model_c7_ic(extended=true). The
    // extended C7 table (Avrett & Loeser 2008 Table 26) reaches 68 Mm / 1.59 MK;
    // this lifts the upper boundary into a real resolved coronal VOLUME (so the
    // evaporated mass has somewhere to accumulate — the coronal EM rise that *is*
    // gentle evaporation), with q(T) imposed at that high coronal top. 0 (default)
    // ⇒ use the full table top. Set by model_gentle from GENTLE_TOP_KM.
    float corona_top_km           = 0.0f;

    // Absolute height [km] of the domain base, used only to label the output
    // heights chromo_main writes (cumulative ds + base). Defaults to the C7 base
    // (1003 km) so every existing scenario's output is unchanged; model_column
    // sets it to its actual base (0 km for the photosphere-anchored C7 IC).
    float out_base_km             = 1003.0f;

    // TWO-FLUID ONLY. Transition-Region Adaptive Conduction (TRAC; Johnston et al.
    // 2019, 2020). The single-fluid release conduction operator has no TRAC factor
    // and ignores these fields entirely.
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
    float trac_Tc_max_frac = 0.2f;  // upper bound on T_c as a fraction of T_peak
                                    // (Johnston 2020 Eq. 9: 0.2). Raise toward 1 on a
                                    // corona-less domain (e.g. iso_t22k, T_peak≈22 kK)
                                    // so the cap does not pin T_c to the T_b floor.

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
    // Heating is added to the electron/ion thermal pool (p_i). TWO-FLUID ONLY:
    // the release timestep has no beam stage and ignores this flag. Default off so
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

    // --- ambient coronal (footpoint) heating (gentle conduction-driven evaporation) ---
    // Steady volumetric heating H(s) [W/m³] standing in for the (sub-grid) coronal
    // heating mechanism that sustains a loop against conduction + radiation — the
    // Rosner–Tucker–Vaiana (1978) static energy balance
    //   d/ds(κ_e T^{5/2} dT/ds) + H(s) − n²Λ(T) = 0.
    // make_loop_dat.py supplies the hydrostatic line; this term supplies H(s), the
    // missing ingredient that makes a resolved chromosphere→TR→corona column a true
    // steady state rather than a draining, prescribed-T atmosphere. Footpoint-
    // anchored exponential (Aschwanden & Schrijver 2002 §3):
    //   H(s) = coronal_heat_E0 · exp(−d(s)/s_H) · enhance(t),
    // where d(s) is the field-aligned distance to the nearest footpoint
    // (one-sided d = s − s₀ for open / half-loop; two-sided d = min over both ends
    // for a full loop, set by coronal_heat_two_sided). The uniform-heating limit
    // is recovered as s_H → ∞ (Martens 2010: T(s) is only weakly sensitive to the
    // shape). UNLIKE the transient beam this is ALWAYS on for a steady run; gentle
    // evaporation (Antiochos & Sturrock 1978) is driven by a slow time RAMP of the
    // amplitude (coronal_heat_enhance > 1), kept sub-threshold so v ≪ c_s.
    //
    // Heating is deposited into the charged thermal pool exactly like the beam
    // (apply_coronal_heating_stage): shared with neutrals by heat capacity in the
    // single-T baseline, into the electron pool when enable_Te (electrons conduct).
    // TWO-FLUID ONLY: the release timestep has no volumetric heating stage and
    // ignores this flag.
    // Default off so existing scenarios / the test suite are byte-for-byte unchanged.
    bool  enable_coronal_heating  = false;
    float coronal_heat_E0         = 0.0f;     // E_H0 [W/m³], footpoint heating amplitude
    float coronal_heat_sH         = 1.0e7f;   // heating scale length s_H [m]
    float coronal_heat_s0         = 0.0f;     // footpoint anchor s₀ [m] (arc length, cell-0 face = 0)
    bool  coronal_heat_two_sided  = false;    // true: decay from BOTH ends (full loop)
    // Phase-3 ramp: multiply E_H0 by a smooth flat-top window that rises from 1 to
    // coronal_heat_enhance over [t_on, t_on+ramp]. enhance = 1 (default) ⇒ no ramp
    // (steady relaxation); enhance > 1 raises the heating to drive the upflow.
    float coronal_heat_t_on       = 0.0f;     // ramp onset time [s]
    float coronal_heat_ramp       = 1.0f;     // cosine ramp half-width [s]
    float coronal_heat_enhance    = 1.0f;     // peak amplitude multiplier after the ramp

    // Positivity / vacuum floor (apply_flare_floor_stage + the rhs predictor and
    // reconstructed-face floors). Clips ρ,p positive, caps |V|, and — where the
    // ionization fraction f ≥ 0.99 — slaves the trace neutral onto the charged
    // fluid (the exact single-fluid limit). Flare runs switch this on implicitly
    // via enable_beam_heating; STEADY scenarios whose hot upper layers fully
    // ionize (ρ_n→0 → 0/0 velocity decode and a runaway spectral radius → NaN,
    // e.g. the C7 column once the TR ionizes) opt in via enable_vacuum_floor.
    // Default false keeps the test suite and untouched scenarios byte-for-byte.
    bool  enable_vacuum_floor     = false;
    // True when the positivity/vacuum floors should run this step.
    bool  floors_active() const { return enable_beam_heating || enable_vacuum_floor; }

    // When true, the outer boundary uses a transparent (Neumann) outflow on the
    // velocity instead of the slow Mach-0.05 cap, so a supersonic evaporation
    // upflow can leave the domain. Set by flare scenarios. Default false keeps
    // the C7 quiet-Sun Mach-capped outflow. NOTE: a Neumann outflow is ill-posed
    // for INflow (a downflow at the top feeds mass in unphysically); use the
    // Mach cap (below) for runs whose residual is a downflow.
    bool  outer_free_outflow      = false;

    // Mach-number fraction for the capped outer outflow (model_c7_update_bc):
    // V_ghost is clamped to ±outer_mach_cap·c_s(T_TR). Default 0.05 (the quiet-Sun
    // value the C7 tests pin). A gentle-evaporation run raises it (e.g. 0.5) so the
    // subsonic upflow exits without the ill-posed free-outflow inflow — the steady
    // relaxation flow (≪ cap) is unaffected.
    float outer_mach_cap          = 0.05f;

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

    mutable MixtureConductionScratch mixture_conduction_scratch;
    mutable MixtureRhsScratch mixture_rhs_scratch;

    // Diagnostic-only: when true, the NEXT mixture_rhs_explicit evaluation copies
    // its final reconstructed face states, limiter values and numerical mass
    // flux into face_flux_capture. Pure output — no solver stage ever reads the
    // capture back, so the flag cannot change the numerical result. Default false.
    // The driver (chromo_main, CHROMO_FACE_FLUX_DIAG) sets it per step.
    bool capture_face_flux = false;
    mutable MixtureFaceFluxCapture face_flux_capture;

    // Diagnostic-only: when true, every mixture_apply_conduction solve copies
    // its FINAL CONVERGED outer-face conduction quantities into
    // outer_conduction_capture. Pure output; no stage reads it back. Default false.
    // The driver (chromo_main, CHROMO_OUTER_COND_DIAG) sets it once per run.
    bool capture_outer_conduction = false;
    mutable OuterConductionCapture outer_conduction_capture;

    // --- cell-centered & face arrays (length ns) --------------------------
    Vec ds_i;
    Vec B_i, B_imh, B_iph;
    Vec dinvB_ds_i;
    Vec phi_g_imh, phi_g_iph;

    // --- release ghost buffers (length num_of_mixture_eq) -----------------
    // The two hydrodynamic ghost layers of the release solver, in the same
    // (rho, rho u, E) rows as the interior state. The conductive physical-face
    // temperature boundary is a SEPARATE datum (outer_conduction_temperature);
    // ghost states and thermal boundary data are never conflated.
    Vec mix_outer_boundary0, mix_outer_boundary1;
    Vec mix_inner_boundary0, mix_inner_boundary1;

    // --- legacy two-fluid ghost buffers (length num_of_eq) ----------------
    Vec outer_boundary0_i, outer_boundary1_i;
    Vec inner_boundary0_i, inner_boundary1_i;

    // --- packed broadcasts of the above (length n_state) ------------------
    // Cached versions of B_i, B_imh, B_iph, dinvB_ds_i, ds_i, and a scratch
    // dt buffer broadcast across all num_of_eq equations. Recomputed by
    // broadcast() whenever the cell-centered fields change.
    Vec B_state, B_state_imh, B_state_iph;
    Vec dinvB_ds_state, ds_state, dt_state;

    // --- static-mesh metric caches (rebuilt by broadcast() from ds_i) ------
    // Support for a metric-aware, static non-uniform finite-volume mesh (static
    // local refinement — NOT AMR: no runtime regridding). Every quantity here is
    // derived from the canonical cell-width array `ds_i`, taking the cell CENTER
    // at the geometric midpoint of its two faces:
    //   * `s_face`   — face arc length [m] from the inner boundary (s_face[0]=0),
    //                  length ns+1.
    //   * `s_i`      — cell-center arc length [m], length ns.
    //   * `ds_iph_i` — center-to-center distance to the i+1 neighbour [m]. At the
    //                  outer boundary the ghost width mirrors the last cell, so
    //                  ds_iph_i[ns-1] = ds_i[ns-1] (Neumann), matching the legacy
    //                  0.5*(ds_i + ip1(ds_i,SLICE)) formula used in rhs/conduction.
    //   * `ds_imh_i` — center-to-center distance to the i-1 neighbour [m], with the
    //                  inner ghost width mirrored: ds_imh_i[0] = ds_i[0].
    // On a UNIFORM mesh ds_iph_i = ds_imh_i = ds_i, `uniform_mesh` is true, and
    // every spacing-sensitive operator keeps its legacy code path (byte-for-byte
    // identical results). broadcast() sets `uniform_mesh` by scanning ds_i and
    // rejects (throws) any non-positive / non-finite width.
    Vec  s_i, s_face;
    Vec  ds_iph_i, ds_imh_i;
    bool uniform_mesh = true;

    // Packed (n_state) center-to-center distances, used by the non-uniform MUSCL
    // reconstruction and conduction paths. Populated by broadcast() on every mesh
    // (they equal ds_state on a uniform mesh).
    Vec ds_iph_state, ds_imh_state;

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
    /// Snapshot of the geometry the caches above were derived from.
    Vec metrics_ds_i, metrics_B_i, metrics_B_imh, metrics_B_iph,
        metrics_dinvB_ds_i;
    bool metrics_valid = false;
    std::uint64_t metrics_generation_ = 0;

    bool static_metrics_current() const;
};

} // namespace chromosphere
