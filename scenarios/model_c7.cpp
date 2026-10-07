/*!
 * @file scenarios/model_c7.cpp
 * @brief The quiet-Sun Model C7 scenario IC/BC and the Route-B photoionization
 *        closure.
 * @ingroup scenarios
 *
 * The tabulated C7 atmosphere itself (MODEL_C7, c7_full_profile(),
 * c7_full_temperature_pchip()) lives in scenarios/c7_profile.cpp. model_flare
 * and analytic_canopy reuse model_c7_ic(); the two-fluid column reuses
 * c7_route_b_photoionization().
 */
#include "model_c7.hpp"
#include "core/physics.hpp"
#include "scenario.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

namespace chromosphere {

// Photospheric inner-ghost reservoir values for model_c7. kInnerRhoNPinned
// and kInnerTnPinned are captured from the C7 base in model_c7_ic.
//   * ρ_n, T_n       → Dirichlet (set from IC's C7 base)
//   * n_i, T_i, V, U → Neumann (zero gradient from cell 0)
namespace {
// --- inner (photospheric) reservoir, captured from the C7 base in model_c7_ic
Real kInnerRhoIPinned = 0.0f;
Real kInnerRhoNPinned = 0.0f;
Real kInnerTnPinned   = 0.0f;
// Discrete-HSE inner-ghost pressures: built so that V=U=0 is an exact fixed
// point of the discretized momentum balance (Pandey et al. 2024 §3.1). Captured
// in model_c7_ic from cell 0's pressure plus the hydrostatic increment ρg·Δs
// across the inner half-face (the ghost sits below cell 0, so p_ghost > p_0).
Real kInnerPiGhost = 0.0f;
Real kInnerPnGhost = 0.0f;

// --- outer (lower-TR) reservoir, captured from the top of the extended C7
// table (Avrett & Loeser 2008 Table 26 at h ≈ 2153 km). This is the fixed
// TR-base state the upper boundary imposes: RTV (1978) constant-pressure
// T_0 ≈ 2×10⁴ K with ρ(T) from P_TR, acting as the coronal heat reservoir that
// conductively opposes the chromosphere top cooling (the q(T) closure). It is
// FIXED (not slaved to the cooling interior) — that is the essential change
// from the old (T, 2T) cascade, which tracked the collapsing top and so
// supplied no net heat.
Real kOuterTtr   = 0.0f;   // TR-base temperature  [K]
Real kOuterNeTr  = 0.0f;   // TR-base electron / proton density [m^-3]
Real kOuterNnTr  = 0.0f;   // TR-base neutral density (slaved, tiny) [m^-3]

// "New explanation" upper BC (docs/studies/evaporation/gentle_evaporation_downflow.md), HYBRID form:
// the hydrostatic-pressure + EOS-density ghost is applied, but the coronal heat
// keeps entering as model_c7's imposed Neumann flux q(T) — NOT a Dirichlet
// temperature jump. On the 100-cell C7 grid the top cell sits at chromospheric
// ~6.6 kK (the thin TR rise to ~23 kK is unresolved), so any Dirichlet jump hot
// enough to deliver a physical ~30 W/m² over-conducts (κ_e∝T^{5/2}) and throws
// off spurious ±km/s flows — exactly the failure the imposed-flux BC was built to
// avoid. The ghost temperature is therefore CONTINUOUS from the top cell,
// T_ghost0 = a·T_top, T_ghost1 = b·T_ghost0 with a = b = 1 by default (a pure
// hydrostatic/EOS extrapolation); a, b are kept env-tunable (C7_TJUMP_A,
// C7_TJUMP_B) only as an experimentation lever and scale the ghost density via
// the EOS — they do NOT drive conduction. Only used when grid.c7_tr_jump_bc is set.
Real kC7Ajump = 1.0f;
Real kC7Bjump = 1.0f;
}

// Natural cubic spline interpolation, port of interp1 in ModChromosphereTest.f90
static Vec interp1_spline(const Vec& xData, const Vec& yData, const Vec& xVal) {
    const arma::uword n = xData.n_elem;
    Vec h(n - 1);
    Vec alpha(n, arma::fill::zeros);
    Vec a(n), b(n - 1), c(n, arma::fill::zeros), d(n - 1);
    Vec l(n, arma::fill::zeros), mu(n, arma::fill::zeros), z(n, arma::fill::zeros);

    for (arma::uword i = 0; i < n - 1; ++i) h(i) = xData(i + 1) - xData(i);
    for (arma::uword i = 1; i < n - 1; ++i)
        alpha(i) = 3.0f / h(i) * (yData(i + 1) - yData(i))
                 - 3.0f / h(i - 1) * (yData(i) - yData(i - 1));

    l(0) = 1.0f;
    for (arma::uword i = 1; i < n - 1; ++i) {
        l(i)  = 2.0f * (xData(i + 1) - xData(i - 1)) - h(i - 1) * mu(i - 1);
        mu(i) = h(i) / l(i);
        z(i)  = (alpha(i) - h(i - 1) * z(i - 1)) / l(i);
    }
    l(n - 1) = 1.0f;

    for (long j = (long)n - 2; j >= 0; --j) {
        c(j) = z(j) - mu(j) * c(j + 1);
        b(j) = (yData(j + 1) - yData(j)) / h(j) - h(j) * (c(j + 1) + 2.0f * c(j)) / 3.0f;
        d(j) = (c(j + 1) - c(j)) / (3.0f * h(j));
        a(j) = yData(j);
    }

    Vec yVal(xVal.n_elem);
    for (arma::uword i = 0; i < xVal.n_elem; ++i) {
        arma::uword j = 0;
        while (j < n - 1 && xVal(i) > xData(j + 1)) ++j;
        if (j >= n - 1) j = n - 2;
        const Real dx = xVal(i) - xData(j);
        yVal(i) = a(j) + b(j) * dx + c(j) * dx * dx + d(j) * dx * dx * dx;
    }
    return yVal;
}

// Route-B photoionization closure (docs/studies/eos-ionization/photoionization_c7_inversion_plan.md;
// writeup §3.1), shared by model_c7 and model_column. Given cell-centered
// (T, n_i, n_n) and heights, returns the per-cell photoionization rate P_phot:
// invert the LOCAL ionization-equilibrium balance of EXACTLY the active Stage-E
// network so the C7 profile is a fixed point,
//   P_phot = (n_i/n_n)(α_r + α_c) n_i − (S_CR + S_i) n_i ,   clipped ≥ 0,
// with the α_c / S_i terms present only when their flags are set. Above the
// equilibrium↔NEQ crossover (h ≈ 1500 km) blend to the frozen FAL-C rate (Chae
// 2021) — forcing equilibrium there is what blew up the earlier calibration.
// n_n is floored at 1 m^-3 so the (n_i/n_n) ratio stays finite where the corona
// ionizes n_n → 0 (the w_eq → 0 weight then discards the irrelevant value); on
// the C7 chromosphere n_n ≥ 1e8, so the floor is a no-op there.
Vec c7_route_b_photoionization(const Grid& grid, const Vec& T_c,
                               const Vec& ni_c, const Vec& nn_c,
                               const Vec& h_cell) {
    const Vec nn_safe = arma::clamp(nn_c, 1.0f, arma::datum::inf);
    Vec a_tot = recombination_rate_alpha(grid, T_c);              // α_r, always on
    if (grid.enable_threebody_recombination)
        a_tot += recombination_rate_alpha_c(grid, T_c, ni_c);     // α_c = κ_c n_e
    Vec S_tot = ionization_rate_S_CR(grid, T_c);                  // S_CR, always on
    if (grid.enable_direct_collisional_ionization)
        S_tot += ionization_rate_S(grid, T_c);                    // S_i, Voronov
    Vec P_inv = (ni_c / nn_safe) % (a_tot % ni_c) - S_tot % ni_c;
    P_inv = arma::clamp(P_inv, 0.0f, arma::datum::inf);
    const Vec  P_falc    = photoionization_rate_chae(h_cell);
    const Real H_EQ_NEQ = 1500.0f;   // km, equilibrium ↔ NEQ crossover (Chae marker)
    const Real BLEND_W  = 120.0f;    // km, blend width
    const Vec w_eq = 0.5f * (1.0f - arma::tanh((h_cell - H_EQ_NEQ) / BLEND_W));
    return w_eq % P_inv + (1.0f - w_eq) % P_falc;
}

Vec model_c7_ic(Grid& grid, bool extended, bool tr_jump_bc) {
    // Domain top: standard column (rows 0..N_C7_BASE-1, h ≈ 2.153 Mm) or the real
    // Table-26 corona extension (rows 0..N_C7-1, up to h ≈ 68 Mm) when the
    // corona-as-boundary gentle-evaporation config asks for it. When extended and
    // grid.corona_top_km > 0, truncate the table at the first row reaching that
    // height — so the upper boundary sits high in a resolved coronal volume at a
    // chosen altitude (default set by model_gentle ≈ 10 Mm) rather than at the very
    // top of the table.
    int n_rows = extended ? N_C7 : N_C7_BASE;
    if (extended && grid.corona_top_km > 0.0f) {
        int nr = N_C7_BASE;   // never go below the standard column
        while (nr < N_C7 && MODEL_C7[nr - 1][0] < grid.corona_top_km) ++nr;
        n_rows = nr;
    }
    // Phase-1 geometry choices (quiet-sun relaxation):
    //   * Uniform B = 1, ∂(1/B)/∂s = 0 — single straight field line; the
    //     flux-tube expansion term in the source/flux drops out trivially.
    //     A non-uniform B(s) (e.g. Spruit 1981 thin-flux-tube) is the
    //     follow-on; the integrator already supports it.
    //   * Gravity ON: φ_g(s) = g · (h(s) - h_base) populated below.
    //
    // HSE strategy (Option A from the Phase-2 plan): we leave C7's tabulated
    // (n, T) as the IC even though C7 was built with NLTE radiative balance,
    // not ideal-gas HSE. Under the code's EOS the C7 profile is off HSE by a
    // nearly uniform 20–26%: |dp/dh + ρg| / |ρg| ≈ 0.25 across the chromosphere
    // (dp/dh is more negative than -ρg, i.e. C7 is "lighter at the top" than
    // a pure-HSE atmosphere with the same T(h) would be). The run absorbs
    // this imbalance during the first few acoustic crossings; the residual
    // drift away from C7 is dominated by missing radiative cooling, not by
    // the IC choice.
    grid.B_imh.ones();
    grid.B_iph.ones();
    grid.B_i.ones();
    grid.dinvB_ds_i.zeros();

    // Tabulated profile (n_rows: chromosphere + TR, Table 26).
    Vec h_tab(n_rows), ne_tab(n_rows), nn_tab(n_rows), T_tab(n_rows);
    for (arma::uword i = 0; i < (arma::uword)n_rows; ++i) {
        h_tab(i)  = MODEL_C7[i][0];
        ne_tab(i) = MODEL_C7[i][1];
        nn_tab(i) = MODEL_C7[i][2];
        T_tab(i)  = MODEL_C7[i][3];
    }

    // Faces: indices -1..ns+3 (size ns+5), uniformly over [h_tab(0), h_tab(top)].
    const arma::uword nF = grid.ns + 5;
    Vec h_F(nF);
    const Real h0 = h_tab(0), h1 = h_tab(n_rows - 1);
    for (arma::uword k = 0; k < nF; ++k) {
        const Real frac = (Real)k / (Real)(nF - 1);
        h_F(k) = h0 + (h1 - h0) * frac;
    }

    // φ_g at each face, measured relative to the C7 base height (h0 = 1003 km)
    // so the absolute potential stays comparable to thermal energy in float32.
    // Only ∂φ_g/∂s matters physically; the constant offset is gauge.
    Vec phi_g_F(nF);
    for (arma::uword k = 0; k < nF; ++k) {
        phi_g_F(k) = grid.g * (h_F(k) - h0) * 1000.0f; // km -> m
    }
    for (arma::uword i = 0; i < grid.ns; ++i) {
        grid.phi_g_imh(i) = phi_g_F(i);
        grid.phi_g_iph(i) = phi_g_F(i + 1);
    }

    Vec ne_F = interp1_spline(h_tab, ne_tab, h_F);
    Vec nn_F = interp1_spline(h_tab, nn_tab, h_F);
    Vec T_F  = interp1_spline(h_tab, T_tab,  h_F);
    // The natural cubic spline can overshoot to negative values across the steep
    // TR where n_HI drops ~2 orders of magnitude over a few km; clamp to small
    // positive floors so densities/temperatures stay physical on every face.
    ne_F = arma::clamp(ne_F, 1.0e10f, arma::datum::inf);
    nn_F = arma::clamp(nn_F, 1.0e8f,  arma::datum::inf);
    T_F  = arma::clamp(T_F,  3.0e3f,  arma::datum::inf);

    // Extended-cell quantities: averages of adjacent faces, indices 0..ns+2.
    const arma::uword nE = grid.ns + 3;
    Vec len_E(nE), ne_E(nE), nn_E(nE), T_E(nE);
    for (arma::uword k = 0; k < nE; ++k) {
        len_E(k) = (h_F(k + 1) - h_F(k)) * 1000.0f; // km -> m
        ne_E(k)  = 0.5f * (ne_F(k + 1) + ne_F(k));
        nn_E(k)  = 0.5f * (nn_F(k + 1) + nn_F(k));
        T_E(k)   = 0.5f * (T_F(k + 1) + T_F(k));
    }

    // Cell-centered interior: first ns of the extended array.
    for (arma::uword i = 0; i < grid.ns; ++i) grid.ds_i(i) = len_E(i);

    // Build the conserved-variable state vector. Total energy includes the
    // ρ·φ_g term so that cons2prim (state.cpp:108) recovers the right pressure.
    Vec xn = arma::zeros<Vec>(grid.n_state);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const Real n_e   = ne_E(i);
        const Real n_n   = nn_E(i);
        const Real T     = T_E(i);
        const Real phi_g = 0.5f * (grid.phi_g_imh(i) + grid.phi_g_iph(i));
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::RHO_I)) = n_e * grid.m_i;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::RHO_N)) = n_n * grid.m_n;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::MOM_I)) = 0.0f;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::MOM_N)) = 0.0f;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::E_I))   = grid.inv_gm1() * grid.k_b * n_e * 2.0f * T + n_e * grid.m_i * phi_g;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::E_N))   = grid.inv_gm1() * grid.k_b * n_n * T       + n_n * grid.m_n * phi_g;
        // Electron internal energy ε_e = p_e/(γ−1) with T_e = T at IC (no KE/φ).
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::E_E))   = grid.inv_gm1() * grid.k_b * n_e * T;
    }

    // Outer ghost cells: the fixed lower-TR reservoir. Both ghosts hold the
    // same TR-base state taken from the top of the extended C7 table (h ≈ 2153
    // km, T ≈ 2.3×10⁴ K) — RTV (1978) constant-pressure TR base. This is the
    // coronal heat reservoir; model_c7_update_bc keeps it FIXED each step (it
    // does NOT track the cooling interior), which is what supplies the net
    // downward conductive flux q(T) that holds the chromosphere top against
    // ionization cooling. φ_g is the upper face of the last interior cell.
    {
        const Real phi_g_outer = grid.phi_g_iph(grid.ns - 1);
        kOuterTtr  = T_tab(n_rows - 1);  // standard 2.310×10⁴ K; extended 2.925×10⁵ K
        kOuterNeTr = ne_tab(n_rows - 1); // electron/proton density at the top
        kOuterNnTr = nn_tab(n_rows - 1); // neutral density at the top (slaved, ~fully ionized)
        for (arma::uword k = 0; k < 2; ++k) {
            Vec& ob = (k == 0 ? grid.outer_boundary0_i : grid.outer_boundary1_i);
            ob(cons::RHO_I) = kOuterNeTr * grid.m_i;
            ob(cons::RHO_N) = kOuterNnTr * grid.m_n;
            ob(cons::MOM_I) = 0.0f;
            ob(cons::MOM_N) = 0.0f;
            ob(cons::E_I)   = grid.inv_gm1() * grid.k_b * kOuterNeTr * 2.0f * kOuterTtr + kOuterNeTr * grid.m_i * phi_g_outer;
            ob(cons::E_N)   = grid.inv_gm1() * grid.k_b * kOuterNnTr *        kOuterTtr + kOuterNnTr * grid.m_n * phi_g_outer;
            ob(cons::E_E)   = grid.inv_gm1() * grid.k_b * kOuterNeTr *        kOuterTtr;
        }
    }

    // Inner ghost cells replicate the innermost extended cell. φ_g taken as
    // the lower face of cell 0 (which is 0 by gauge choice), matching the
    // im1 path in rhs.cpp:178 and integrators.cpp:189.
    //
    // The (n_n, T) snapshot here is what model_c7_update_bc later uses as
    // the Dirichlet-pinned photospheric reservoir state for the inner ghost.
    {
        const Real n_e = ne_E(0);
        const Real n_n = nn_E(0);
        const Real T   = T_E(0);
        const Real ds0 = len_E(0);
        const Real phi_g_inner = grid.phi_g_imh(0);
        grid.inner_boundary0_i(cons::RHO_I) = n_e * grid.m_i;
        grid.inner_boundary0_i(cons::RHO_N) = n_n * grid.m_n;
        grid.inner_boundary0_i(cons::MOM_I) = 0.0f;
        grid.inner_boundary0_i(cons::MOM_N) = 0.0f;
        grid.inner_boundary0_i(cons::E_I)   = grid.inv_gm1() * grid.k_b * n_e * 2.0f * T + n_e * grid.m_i * phi_g_inner;
        grid.inner_boundary0_i(cons::E_N)   = grid.inv_gm1() * grid.k_b * n_n * T       + n_n * grid.m_n * phi_g_inner;
        grid.inner_boundary0_i(cons::E_E)   = grid.inv_gm1() * grid.k_b * n_e * T;
        grid.inner_boundary1_i = grid.inner_boundary0_i;

        kInnerRhoIPinned = n_e * grid.m_i;
        kInnerRhoNPinned = n_n * grid.m_n;
        kInnerTnPinned   = T;

        // Discrete-HSE ghost pressures (Pandey 2024 §3.1): the inner ghost sits
        // a cell-spacing ds0 below cell 0, so a hydrostatic reservoir carries the
        // extra weight ρg·ds0 → p_ghost = p_0 + ρ_0 g ds0. Imposing these (rather
        // than a zero-gradient copy of p_0) makes the pressure gradient across
        // the inner face balance gravity discretely, so V=U=0 is an exact fixed
        // point of the momentum update and the base does not launch a wave.
        const Real p_i_0 = 2.0f * n_e * grid.k_b * T;
        const Real p_n_0 =        n_n * grid.k_b * T;
        kInnerPiGhost = p_i_0 + (n_e * grid.m_i) * grid.g * ds0;
        kInnerPnGhost = p_n_0 + (n_n * grid.m_n) * grid.g * ds0;
    }

    // Height-dependent photoionization rate P_phot(h): Route B closure
    // (docs/studies/eos-ionization/photoionization_c7_inversion_plan.md; writeup §3.1).
    //
    // Equilibrium region (dense lower/mid chromosphere): set P_phot by
    // inverting the *local ionization-equilibrium balance of exactly the Stage E
    // network* on C7's tabulated (n_e, n_n, T), so C7 is a self-consistent
    // fixed point of the rate network Stage E integrates under whatever flags
    // are active:
    //   P_phot(z) = (n_e/n_n)(α_r + α_c) n_e − (S_i + S_CR) n_e ,   clipped ≥ 0,
    // where the α_c (three-body) and S_i (Voronov) terms are present only when
    // their flags are set. The default network is P_phot + S_CR (ionization) and
    // α_r (recombination); S_i is 2–3 orders below S_CR everywhere and α_c is
    // within 1–2 orders of α_r only below ~700 km (outside C7's 1003 km floor),
    // so dropping them barely perturbs the inverted P_phot — but the inversion
    // must still match the active network exactly, or C7 drifts off its own IC.
    // This closes the over-ionization gap where it physically should (the
    // ionization timescale is short and equilibrium holds). The inverted rate
    // runs a few × above the ground-state FAL-C rate of Chae (2021) — the
    // expected case-B-vs-ground-state recombination difference, not a defect.
    //
    // NEQ region (rarefied upper chromosphere, z ≳ 1500 km): the recombination
    // timescale exceeds the dynamical time, so the gas is genuinely out of
    // ionization equilibrium (real over-ionization). There we keep the frozen
    // FAL-C radiation-field rate (photoionization_rate_chae) and let the gas
    // depart from C7 — forcing C7 as a fixed point here is what blew up the
    // old calibration. A tanh window in height blends the two.
    {
        Vec h_cell(grid.ns), Tc(grid.ns), nec(grid.ns), nnc(grid.ns);
        for (arma::uword i = 0; i < grid.ns; ++i) {
            h_cell(i) = 0.5f * (h_F(i) + h_F(i + 1));   // km, tau_500=1 scale
            Tc(i)  = T_E(i);
            nec(i) = ne_E(i);
            nnc(i) = nn_E(i);
        }
        grid.photoionization_rate_i = c7_route_b_photoionization(grid, Tc, nec, nnc, h_cell);
    }

    // Upper-boundary coronal conductive flux q(T) (RTV 1978). The static
    // conduction-dominated corona delivers F_c = (2/7) κ₀ T_cor^{7/2} / L to the
    // TR base, where it is re-radiated; we impose that flux as the Stage-D
    // Neumann outer BC. κ₀ = 10⁻¹¹ W m⁻¹ K⁻⁷ᐟ² (Johnston 2020 SI), and a quiet-Sun
    // corona T_cor ≈ 8×10⁵ K over a half-length L ≈ 15 Mm gives F_c ≈ 90 W m⁻².
    //
    // q(T) heating and radiative cooling are a physical pair: the coronal flux is
    // re-radiated at/below the TR base. Impose q(T) only when Stage R cooling is
    // active — otherwise the injected heat has no sink and the top runs away. So
    // the complete physical model_c7 is run WITH cooling; the no-cooling mode is
    // a diagnostic baseline (no coronal heating imposed).
    {
        const Real kappa0 = 1.0e-11f;   // W m^-1 K^-7/2
        if (extended) {
            // Corona-as-boundary (top at h ≈ 2.628 Mm, T ≈ 0.29 MK): impose the
            // conductive flux C7 ITSELF carries across the new top, computed from
            // its own gradient q = κ₀ T^{5/2} dT/ds using the top two Table-26 rows
            // (≈130 W m⁻² downward) — self-consistent with the model, no assumed
            // T_cor / L. This is the steady quiet value; the gentle-evaporation
            // driver (model_gentle) ramps it up via outer_heat_flux_enhance.
            const Real T_top = T_tab(n_rows - 1);
            const Real dTds  = (T_tab(n_rows - 1) - T_tab(n_rows - 2))
                              / ((h_tab(n_rows - 1) - h_tab(n_rows - 2)) * 1.0e3f);  // km→m
            grid.outer_heat_flux = kappa0 * std::pow(T_top, 2.5f) * dTds;
        } else {
            const Real T_cor  = 8.0e5f;     // K, quiet-Sun coronal apex
            const Real L_cor  = 4.5e7f;     // m, coronal loop half-length (quiet-Sun)
            // ≈ 29 W m⁻². With the optically-thin TR radiative sink present, this
            // sits in the stable window q ∈ [~30, 90] W m⁻²: the TR/corona radiates
            // the imposed flux and the chromosphere proper stays cool. Below ~10 W m⁻²
            // the weak heating drops the upper chromosphere into the radiative-loss-peak
            // thermal-instability (condensation) regime and the run fails.
            grid.outer_heat_flux = (2.0f / 7.0f) * kappa0
                                 * std::pow(T_cor, 3.5f) / L_cor;
        }
        // Steady base for the q(T) ramp (model_gentle); enhance = 1 ⇒ q ≡ base.
        grid.outer_heat_flux_base   = grid.outer_heat_flux;
        grid.impose_outer_heat_flux = grid.enable_radiative_cooling;
    }

    // Isotropic numerical diffusivity (Pandey 2024 §3.3), grid-scaled. Pandey's
    // benchmark is χ ≳ 30×10⁸ m² s⁻¹ at Δz = 40 km; here it only needs to damp the
    // base wave and any TR-gradient ringing (the coronal heating is carried by
    // q(T), not by χ), so a modest value suffices. Scale ∝ Δs so coarser grids
    // diffuse over a comparable number of cells.
    {
        grid.numerical_diffusivity_per_length = 2.0e3f; // C_num [m/s]
    }

    // TRAC (Johnston et al. 2020): broaden the under-resolved TR so the coarse
    // 100-cell grid recovers the grid-converged coronal density response while
    // conserving the TR-integrated radiation/heating (the Eq.-27 jump balance).
    // Paired with the radiative sink (q(T)/Λ(T) need a resolved TR to balance).
    grid.enable_trac   = grid.enable_radiative_cooling;
    grid.trac_T_chrom  = 2.0e4f;

    // Vacuum / positivity floor (the same machinery the flare scenario uses).
    // As the C7 column relaxes, the upper chromosphere/TR fully ionizes and the
    // neutral density n_n collapses toward zero there; the 0/0 neutral-velocity
    // decode then makes the spectral radius blow up → NaN (the run otherwise dies
    // around step ~2000). The floor clips ρ,p positive and, where f = n_i/(n_i+n_n)
    // ≥ 0.99, slaves the trace neutral onto the charged fluid (U=V, T_n=T_e) — the
    // single-fluid limit the outer BC already imposes. It is a no-op below the TR
    // where the gas is partially neutral and well-conditioned.
    grid.enable_vacuum_floor = true;

    // "New explanation" upper BC (docs/studies/evaporation/gentle_evaporation_downflow.md), HYBRID
    // form, gated to the bare model_c7 scenario (model_flare / analytic_canopy /
    // model_gentle keep the default RTV reservoir):
    //   * well-balanced explicit reconstruction (removes the spurious (γ−1)g
    //     downforce that the old run damped into the persistent ~1 km/s TR
    //     downflow — see the doc's "load-bearing discovery"),
    //   * the top face imposes a HYDROSTATIC ghost pressure + EOS ghost density
    //     (model_c7_update_bc, c7_tr_jump_bc branch) so V=0 is a boundary fixed
    //     point and there is no boundary downflow.
    // The coronal heat keeps entering via the imposed Neumann flux q(T)
    // (impose_outer_heat_flux left ON with cooling) — NOT a Dirichlet temperature
    // jump, which over-conducts on the coarse C7 grid. So we do NOT touch
    // impose_outer_heat_flux here.
    if (tr_jump_bc) {
        auto env_f = [](const char* key, Real fallback) -> Real {
            if (const char* e = std::getenv(key)) {
                try { return std::stof(std::string(e)); } catch (...) {}
            }
            return fallback;
        };
        kC7Ajump = env_f("C7_TJUMP_A", 1.0f);   // ghost-T multiplier on the live top T
        kC7Bjump = env_f("C7_TJUMP_B", 1.0f);   // (default 1 ⇒ continuous extrapolation)
        grid.c7_tr_jump_bc = true;
        grid.well_balanced = true;
    }

    grid.broadcast();
    return xn;
}

void model_c7_update_bc(Grid& grid, const Vec& xn) {
    // Gentle-evaporation driver (docs/studies/evaporation/gentle_evaporation_plan.md v2): ramp the
    // imposed coronal conductive flux q(T) = q_base · enhance(t). enhance(t) holds
    // at 1 until t_on, then cosine-rises to outer_heat_flux_enhance over
    // [t_on, t_on+ramp]. enhance = 1 (default) ⇒ q stays at the steady base, so
    // model_c7 / model_flare are byte-for-byte unchanged. The stronger downward
    // flux heats the upper TR → conduction-driven upflow (Antiochos & Sturrock 1978).
    if (grid.outer_heat_flux_base > 0.0f) {
        Real amp = 1.0f;
        if (grid.outer_heat_flux_enhance != 1.0f) {
            const Real t   = grid.sim_time;
            const Real ton = grid.outer_heat_flux_t_on;
            const Real r   = (grid.outer_heat_flux_ramp > 1.0e-6f) ? grid.outer_heat_flux_ramp : 1.0e-6f;
            Real w;
            if (t <= ton)          w = 0.0f;
            else if (t >= ton + r) w = 1.0f;
            else w = 0.5f * (1.0f - std::cos(static_cast<Real>(arma::datum::pi) * (t - ton) / r));
            amp = 1.0f + (grid.outer_heat_flux_enhance - 1.0f) * w;
        }
        grid.outer_heat_flux = grid.outer_heat_flux_base * amp;
    }

    // Physically-grounded boundary closures (docs/studies/boundaries/boundary_conditions_plan.md):
    //   * Outer face — fixed lower-TR reservoir (RTV 1978 + Bradshaw & Emslie
    //     2020 + the single-fluid slaving of Gómez Míguez et al. 2024).
    //   * Inner face — discrete-HSE photospheric reservoir, V=U=0 (Pandey 2024).
    // apply_open_bcs lays down a default ghost layer; we overwrite both faces.
    apply_open_bcs(grid, xn);

    const auto  sz        = arma::size(grid.ns, num_of_eq);
    const Real m_i       = grid.m_i;
    const Real m_n       = grid.m_n;
    const Real k_b       = grid.k_b;

    // ====================================================================
    // Outer face — lower transition region (h ≈ 2153 km, T ≈ 2.3×10⁴ K).
    // ====================================================================
    if (grid.c7_tr_jump_bc) {
        // "New explanation" upper BC (docs/studies/evaporation/gentle_evaporation_downflow.md), HYBRID
        // form. The boundary downflow is removed by the well-balanced
        // reconstruction + a hydrostatic/EOS ghost; the coronal HEAT keeps coming
        // from the imposed Neumann flux q(T) (not a Dirichlet T-jump, which over-
        // conducts on the coarse C7 grid). Three steps:
        //   (1) PRESSURE — hydrostatic ghost so the centered gradient at the
        //       boundary cell balances gravity, (p_{nl-1} − p_ghost)/(2Δs) = ρ_nl g
        //       ⇒ V=0 is a fixed point (paired with well_balanced, no boundary
        //       downflow). Applied per fluid (ions carry the electron pressure).
        //   (2) TEMPERATURE — CONTINUOUS from the top cell, T_ghost0 = a·T_top,
        //       T_ghost1 = b·T_ghost0 with a=b=1 by default (no jump). The ghost-T
        //       sets only the EOS density and the hydro flux; it does NOT drive
        //       conduction (the imposed flux does), so there is no over-conduction.
        //   (3) DENSITY — from the EOS at (p_ghost, T_ghost): p_i = 2 n_i k T
        //       (protons+electrons), p_n = n_n k T.
        // Velocity: Mach-capped OUTFLOW (lets a conduction-driven upflow leave; the
        // cap |V| ≤ outer_mach_cap·c_s prevents the ill-posed-inflow runaway).
        const Real phi_g_out = grid.phi_g_iph(grid.ns - 1);
        const Real g         = grid.g;

        // Decode an interior cell's per-fluid pressures / densities (cons2prim).
        auto decode = [&](arma::uword i, Real& p_i, Real& p_n,
                          Real& rho_i, Real& rho_n) {
            rho_i = xn(arma::sub2ind(sz, i, cons::RHO_I));
            rho_n = xn(arma::sub2ind(sz, i, cons::RHO_N));
            const Real Vi    = xn(arma::sub2ind(sz, i, cons::MOM_I)) / rho_i;
            const Real Un    = xn(arma::sub2ind(sz, i, cons::MOM_N)) / rho_n;
            const Real E_i   = xn(arma::sub2ind(sz, i, cons::E_I));
            const Real E_n   = xn(arma::sub2ind(sz, i, cons::E_N));
            const Real phi_g = 0.5f * (grid.phi_g_imh(i) + grid.phi_g_iph(i));
            p_i = grid.gm1() * E_i - grid.half_gm1() * rho_i * Vi * Vi - grid.gm1() * rho_i * phi_g;
            p_n = grid.gm1() * E_n - grid.half_gm1() * rho_n * Un * Un - grid.gm1() * rho_n * phi_g;
        };

        const arma::uword nl = grid.ns - 1;
        const Real ds = grid.ds_i(nl);
        Real p_i_top, p_n_top, rho_i_top, rho_n_top;
        decode(nl, p_i_top, p_n_top, rho_i_top, rho_n_top);
        Real p_i_2, p_n_2, rho_i_2, rho_n_2;
        decode(grid.ns >= 2 ? nl - 1 : nl, p_i_2, p_n_2, rho_i_2, rho_n_2);

        // (2) ghost temperature, continuous from the top cell (a=b=1 default).
        const Real T_top = p_i_top / (2.0f * (rho_i_top / m_i) * k_b);  // charged-fluid T
        const Real T_g0  = kC7Ajump * T_top;
        const Real T_g1  = kC7Bjump * T_g0;

        // (1)+(3) per-fluid hydrostatic ghost pressure → EOS density at (p, T).
        auto build = [&](Real p_i_lo, Real p_n_lo, Real rho_i_w, Real rho_n_w,
                         Real T_g, Real& rho_i_g, Real& rho_n_g,
                         Real& p_i_g, Real& p_n_g) {
            p_i_g = p_i_lo - 2.0f * ds * rho_i_w * g;
            p_n_g = p_n_lo - 2.0f * ds * rho_n_w * g;
            if (p_i_g < 1.0e-6f)  p_i_g = 1.0e-6f;     // keep positive at the tenuous top
            if (p_n_g < 1.0e-12f) p_n_g = 1.0e-12f;
            rho_i_g = (p_i_g / (2.0f * k_b * T_g)) * m_i;   // p_i = 2 n_i k T
            rho_n_g = (p_n_g / (k_b * T_g)) * m_n;          // p_n = n_n k T
        };
        Real rho_i_g0, rho_n_g0, p_i_g0, p_n_g0;
        build(p_i_2,  p_n_2,  rho_i_top, rho_n_top, T_g0, rho_i_g0, rho_n_g0, p_i_g0, p_n_g0);
        Real rho_i_g1, rho_n_g1, p_i_g1, p_n_g1;
        build(p_i_top, p_n_top, rho_i_g0, rho_n_g0, T_g1, rho_i_g1, rho_n_g1, p_i_g1, p_n_g1);

        // velocity: Mach-capped outflow (lifted for flare's free outflow).
        const Real c_s   = std::sqrt(2.0f * grid.gamma_mono * k_b * kOuterTtr / m_i);
        const Real V_cap = grid.outer_mach_cap * c_s;
        Real V_g = xn(arma::sub2ind(sz, nl, cons::MOM_I)) / rho_i_top;
        if (!grid.outer_free_outflow) {
            if (V_g >  V_cap) V_g =  V_cap;
            if (V_g < -V_cap) V_g = -V_cap;
        }

        auto pack = [&](Vec& ob, Real rho_i, Real rho_n, Real p_i, Real p_n) {
            ob(cons::RHO_I) = rho_i;
            ob(cons::RHO_N) = rho_n;
            ob(cons::MOM_I) = rho_i * V_g;
            ob(cons::MOM_N) = rho_n * V_g;       // neutral locked to charge fluid
            ob(cons::E_I)   = grid.inv_gm1() * p_i + 0.5f * rho_i * V_g * V_g + rho_i * phi_g_out;
            ob(cons::E_N)   = grid.inv_gm1() * p_n + 0.5f * rho_n * V_g * V_g + rho_n * phi_g_out;
            ob(cons::E_E)   = 0.5f * grid.inv_gm1() * p_i;   // ε_e = p_e/(γ−1), p_e = ½ p_i ⇒ T_e = T_i
        };
        pack(grid.outer_boundary0_i, rho_i_g0, rho_n_g0, p_i_g0, p_n_g0);
        pack(grid.outer_boundary1_i, rho_i_g1, rho_n_g1, p_i_g1, p_n_g1);
    } else {
        // Default: fixed lower-TR reservoir (RTV 1978). The TR base is a constant-
        // pressure layer: pressure is continuous across it while temperature jumps
        // to the coronal value over the (here unresolved) TR. We do NOT impose a
        // hot dense ghost wall — on the coarse grid a fixed coronal T would over-
        // conduct (κ_e∝T^{5/2}) and over-pressurize the top cell, driving a
        // spurious supersonic downflow. Instead:
        //   * ρ, T  — Neumann (continuous pressure across the constant-pressure TR;
        //             apply_open_bcs already laid these down).
        //   * V, U  — mass-flux outflow capped at a slow Mach number M ≤ 0.05 of the
        //             TR sound speed c_s = √(2γ k_B T/m_p) (Bradshaw & Emslie 2020;
        //             |V| ≲ 1 km/s, sign/magnitude as in Hansteen 1993). U = V locks
        //             the trace neutral to the charged fluid (Gómez Míguez 2024,
        //             single-fluid TR) — no ill-conditioned U_n = ρU_n/ρ_n.
        // The TR thermal structure / coronal support enters entirely through the
        // imposed downward conductive flux q(T) (grid.outer_heat_flux, set in
        // model_c7_ic and applied as the Stage-D Neumann BC) — the load-bearing
        // closure that balances Route-B ionization cooling.
        const Real phi_g_out = grid.phi_g_iph(grid.ns - 1);
        const Real c_s   = std::sqrt(2.0f * grid.gamma_mono * k_b * kOuterTtr / m_i);
        const Real V_cap = grid.outer_mach_cap * c_s;

        auto cap_ghost = [&](Vec& ob) {
            const Real rho_i = ob(cons::RHO_I);
            const Real rho_n = ob(cons::RHO_N);
            const Real n_i   = rho_i / m_i;
            const Real n_n   = rho_n / m_n;
            const Real V_old = ob(cons::MOM_I) / rho_i;
            const Real U_old = ob(cons::MOM_N) / rho_n;
            // Decode Neumann T_i, T_n from the energy apply_open_bcs packed.
            const Real T_i_in = (ob(cons::E_I) - 0.5f * rho_i * V_old * V_old
                                  - rho_i * phi_g_out) / (2.0f * grid.inv_gm1() * k_b * n_i);
            const Real T_n_in = (ob(cons::E_N) - 0.5f * rho_n * U_old * U_old
                                  - rho_n * phi_g_out) / (grid.inv_gm1() * k_b * n_n);
            // Mach-0.05 cap for the quiet-Sun TR (Bradshaw & Emslie 2020). For
            // flare scenarios (grid.outer_free_outflow) the cap is lifted so a
            // supersonic evaporation upflow can leave the domain through a
            // transparent (Neumann) outer face — clamping it would choke the
            // ~100s km/s evaporated outflow.
            Real V_g = V_old;
            if (!grid.outer_free_outflow) {
                if (V_g >  V_cap) V_g =  V_cap;
                if (V_g < -V_cap) V_g = -V_cap;
            }
            const Real U_g = V_g;                // neutral locked to charge fluid
            ob(cons::MOM_I) = rho_i * V_g;
            ob(cons::MOM_N) = rho_n * U_g;
            ob(cons::E_I)   = grid.inv_gm1() * k_b * n_i * 2.0f * T_i_in
                              + 0.5f * rho_i * V_g * V_g + rho_i * phi_g_out;
            ob(cons::E_N)   = grid.inv_gm1() * k_b * n_n *        T_n_in
                              + 0.5f * rho_n * U_g * U_g + rho_n * phi_g_out;
            ob(cons::E_E)   = grid.inv_gm1() * k_b * n_i *        T_i_in;   // T_e = T_charged ghost
        };
        cap_ghost(grid.outer_boundary0_i);
        cap_ghost(grid.outer_boundary1_i);
    }   // end outer-face mode branch

    // ====================================================================
    // Inner face — photosphere, discrete-HSE V=0 reservoir (Pandey 2024).
    // ====================================================================
    // Dirichlet ρ_i, ρ_n at the C7 base; V=U=0; and hydrostatic ghost pressures
    // p_ghost = p_0 + ρ_0 g Δs (captured in model_c7_ic) so the pressure
    // gradient across the inner face balances gravity discretely and V=U=0 is an
    // exact fixed point — removing the spurious base wave the old "pin-to-IC"
    // ghost launched. Numerical diffusion (Stage D isotropic floor) damps any
    // residual ringing rather than reflecting it.
    {
        const Real phi_g_inner = grid.phi_g_imh(0);
        const Real rho_i_g = kInnerRhoIPinned;
        const Real rho_n_g = kInnerRhoNPinned;

        Vec& ob0 = grid.inner_boundary0_i;
        ob0(cons::RHO_I) = rho_i_g;
        ob0(cons::RHO_N) = rho_n_g;
        ob0(cons::MOM_I) = 0.0f;                  // V = 0 reservoir
        ob0(cons::MOM_N) = 0.0f;                  // U = 0 reservoir
        ob0(cons::E_I)   = grid.inv_gm1() * kInnerPiGhost + rho_i_g * phi_g_inner;
        ob0(cons::E_N)   = grid.inv_gm1() * kInnerPnGhost + rho_n_g * phi_g_inner;
        ob0(cons::E_E)   = 0.5f * grid.inv_gm1() * kInnerPiGhost;   // p_e = ½ p_total ⇒ T_e = T_i
        grid.inner_boundary1_i = ob0;
    }

    grid.broadcast();
}

} // namespace chromosphere
