#include "model_c7.hpp"
#include "physics.hpp"
#include "scenario.hpp"

#include <cmath>

namespace chromosphere {

// Photospheric inner-ghost reservoir values for model_c7. kInnerRhoNPinned
// and kInnerTnPinned are captured from the C7 base in model_c7_ic.
//   * ρ_n, T_n       → Dirichlet (set from IC's C7 base)
//   * n_i, T_i, V, U → Neumann (zero gradient from cell 0)
namespace {
// --- inner (photospheric) reservoir, captured from the C7 base in model_c7_ic
float kInnerRhoIPinned = 0.0f;
float kInnerRhoNPinned = 0.0f;
float kInnerTnPinned   = 0.0f;
// Discrete-HSE inner-ghost pressures: built so that V=U=0 is an exact fixed
// point of the discretized momentum balance (Pandey et al. 2024 §3.1). Captured
// in model_c7_ic from cell 0's pressure plus the hydrostatic increment ρg·Δs
// across the inner half-face (the ghost sits below cell 0, so p_ghost > p_0).
float kInnerPiGhost = 0.0f;
float kInnerPnGhost = 0.0f;

// --- outer (lower-TR) reservoir, captured from the top of the extended C7
// table (Avrett & Loeser 2008 Table 26 at h ≈ 2153 km). This is the fixed
// TR-base state the upper boundary imposes: RTV (1978) constant-pressure
// T_0 ≈ 2×10⁴ K with ρ(T) from P_TR, acting as the coronal heat reservoir that
// conductively opposes the chromosphere top cooling (the q(T) closure). It is
// FIXED (not slaved to the cooling interior) — that is the essential change
// from the old (T, 2T) cascade, which tracked the collapsing top and so
// supplied no net heat.
float kOuterTtr   = 0.0f;   // TR-base temperature  [K]
float kOuterNeTr  = 0.0f;   // TR-base electron / proton density [m^-3]
float kOuterNnTr  = 0.0f;   // TR-base neutral density (slaved, tiny) [m^-3]
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
        const float dx = xVal(i) - xData(j);
        yVal(i) = a(j) + b(j) * dx + c(j) * dx * dx + d(j) * dx * dx * dx;
    }
    return yVal;
}

// Model C7 chromosphere + lower transition region profile.
// Columns: height [km], n_e [m^-3], n_n (= n_HI) [m^-3], T [K].
//
// Rows 0–14 are the original chromospheric table (writeup Table 1); rows 15–31
// extend it up through the lower TR to h ≈ 2153 km, T ≈ 23000 K — just above the
// κ_e/κ_n conductivity crossover (2149 km) — so the upper boundary can sit at a
// physically-grounded TR base (RTV T_0 ≈ 2×10⁴ K) rather than mid-chromosphere.
// All rows are Avrett & Loeser (2008) Table 26 (Model C7, ApJS 175, 229, p. 248):
// n_e ← their n_e column, n_n ← their n_HI column, T ← their T column, with the
// cm⁻³→m⁻³ (×10⁶) conversion applied. n_HI drops ~2 orders of magnitude across
// the thin TR (2138→2153 km) as hydrogen ionizes.
static const int N_C7 = 32;
static const float MODEL_C7[N_C7][4] = {
    {1.003e+03f, 1.903e+17f, 2.693e+19f, 6.225e+03f},
    {1.032e+03f, 2.021e+17f, 2.179e+19f, 6.315e+03f},
    {1.065e+03f, 2.091e+17f, 1.722e+19f, 6.400e+03f},
    {1.101e+03f, 2.104e+17f, 1.340e+19f, 6.474e+03f},
    {1.143e+03f, 1.995e+17f, 1.009e+19f, 6.531e+03f},
    {1.214e+03f, 1.760e+17f, 6.322e+18f, 6.576e+03f},
    {1.299e+03f, 1.489e+17f, 3.641e+18f, 6.598e+03f},
    {1.398e+03f, 1.438e+17f, 1.880e+18f, 6.610e+03f},
    {1.520e+03f, 1.423e+17f, 7.956e+17f, 6.623e+03f},
    {1.617e+03f, 1.267e+17f, 4.012e+17f, 6.633e+03f},
    {1.722e+03f, 1.027e+17f, 1.916e+17f, 6.643e+03f},
    {1.820e+03f, 7.969e+16f, 9.843e+16f, 6.652e+03f},
    {1.894e+03f, 6.450e+16f, 6.005e+16f, 6.660e+03f},
    {1.946e+03f, 5.526e+16f, 4.371e+16f, 6.667e+03f},
    {1.989e+03f, 4.826e+16f, 3.447e+16f, 6.674e+03f},
    // --- lower transition region (Avrett & Loeser 2008 Table 26, rows 87→68) ---
    {2.024e+03f, 4.307e+16f, 2.916e+16f, 6.680e+03f},
    {2.055e+03f, 3.883e+16f, 2.583e+16f, 6.686e+03f},
    {2.083e+03f, 3.498e+16f, 2.351e+16f, 6.694e+03f},
    {2.098e+03f, 3.293e+16f, 2.258e+16f, 6.700e+03f},
    {2.110e+03f, 3.128e+16f, 2.204e+16f, 6.706e+03f},
    {2.120e+03f, 2.985e+16f, 2.176e+16f, 6.718e+03f},
    {2.130e+03f, 2.826e+16f, 2.157e+16f, 6.768e+03f},
    {2.134e+03f, 2.741e+16f, 2.118e+16f, 6.870e+03f},
    {2.138e+03f, 2.608e+16f, 1.954e+16f, 7.248e+03f},
    {2.142e+03f, 2.432e+16f, 1.693e+16f, 7.950e+03f},
    {2.145e+03f, 2.325e+16f, 1.165e+16f, 9.115e+03f},
    {2.147e+03f, 2.274e+16f, 4.413e+15f, 1.098e+04f},
    {2.149e+03f, 2.046e+16f, 1.749e+15f, 1.320e+04f},
    {2.150e+03f, 1.761e+16f, 8.787e+14f, 1.576e+04f},
    {2.151e+03f, 1.550e+16f, 5.412e+14f, 1.814e+04f},
    {2.152e+03f, 1.383e+16f, 3.809e+14f, 2.051e+04f},
    {2.153e+03f, 1.237e+16f, 2.769e+14f, 2.310e+04f}
};

Vec model_c7_ic(Grid& grid) {
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

    // Tabulated profile (N_C7 rows: chromosphere + lower TR, Table 26).
    Vec h_tab(N_C7), ne_tab(N_C7), nn_tab(N_C7), T_tab(N_C7);
    for (arma::uword i = 0; i < (arma::uword)N_C7; ++i) {
        h_tab(i)  = MODEL_C7[i][0];
        ne_tab(i) = MODEL_C7[i][1];
        nn_tab(i) = MODEL_C7[i][2];
        T_tab(i)  = MODEL_C7[i][3];
    }

    // Faces: indices -1..ns+3 (size ns+5), uniformly over [h_tab(0), h_tab(top)].
    const arma::uword nF = grid.ns + 5;
    Vec h_F(nF);
    const float h0 = h_tab(0), h1 = h_tab(N_C7 - 1);
    for (arma::uword k = 0; k < nF; ++k) {
        const float frac = (float)k / (float)(nF - 1);
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
        const float n_e   = ne_E(i);
        const float n_n   = nn_E(i);
        const float T     = T_E(i);
        const float phi_g = 0.5f * (grid.phi_g_imh(i) + grid.phi_g_iph(i));
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::RHO_I)) = n_e * grid.m_i;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::RHO_N)) = n_n * grid.m_n;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::MOM_I)) = 0.0f;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::MOM_N)) = 0.0f;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::E_I))   = 1.5f * grid.k_b * n_e * 2.0f * T + n_e * grid.m_i * phi_g;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::E_N))   = 1.5f * grid.k_b * n_n * T       + n_n * grid.m_n * phi_g;
    }

    // Outer ghost cells: the fixed lower-TR reservoir. Both ghosts hold the
    // same TR-base state taken from the top of the extended C7 table (h ≈ 2153
    // km, T ≈ 2.3×10⁴ K) — RTV (1978) constant-pressure TR base. This is the
    // coronal heat reservoir; model_c7_update_bc keeps it FIXED each step (it
    // does NOT track the cooling interior), which is what supplies the net
    // downward conductive flux q(T) that holds the chromosphere top against
    // ionization cooling. φ_g is the upper face of the last interior cell.
    {
        const float phi_g_outer = grid.phi_g_iph(grid.ns - 1);
        kOuterTtr  = T_tab(N_C7 - 1);    // 2.310×10⁴ K
        kOuterNeTr = ne_tab(N_C7 - 1);   // 1.237×10¹⁶ m⁻³
        kOuterNnTr = nn_tab(N_C7 - 1);   // 2.769×10¹⁴ m⁻³ (slaved, ~fully ionized)
        for (arma::uword k = 0; k < 2; ++k) {
            Vec& ob = (k == 0 ? grid.outer_boundary0_i : grid.outer_boundary1_i);
            ob(cons::RHO_I) = kOuterNeTr * grid.m_i;
            ob(cons::RHO_N) = kOuterNnTr * grid.m_n;
            ob(cons::MOM_I) = 0.0f;
            ob(cons::MOM_N) = 0.0f;
            ob(cons::E_I)   = 1.5f * grid.k_b * kOuterNeTr * 2.0f * kOuterTtr + kOuterNeTr * grid.m_i * phi_g_outer;
            ob(cons::E_N)   = 1.5f * grid.k_b * kOuterNnTr *        kOuterTtr + kOuterNnTr * grid.m_n * phi_g_outer;
        }
    }

    // Inner ghost cells replicate the innermost extended cell. φ_g taken as
    // the lower face of cell 0 (which is 0 by gauge choice), matching the
    // im1 path in rhs.cpp:178 and integrators.cpp:189.
    //
    // The (n_n, T) snapshot here is what model_c7_update_bc later uses as
    // the Dirichlet-pinned photospheric reservoir state for the inner ghost.
    {
        const float n_e = ne_E(0);
        const float n_n = nn_E(0);
        const float T   = T_E(0);
        const float ds0 = len_E(0);
        const float phi_g_inner = grid.phi_g_imh(0);
        grid.inner_boundary0_i(cons::RHO_I) = n_e * grid.m_i;
        grid.inner_boundary0_i(cons::RHO_N) = n_n * grid.m_n;
        grid.inner_boundary0_i(cons::MOM_I) = 0.0f;
        grid.inner_boundary0_i(cons::MOM_N) = 0.0f;
        grid.inner_boundary0_i(cons::E_I)   = 1.5f * grid.k_b * n_e * 2.0f * T + n_e * grid.m_i * phi_g_inner;
        grid.inner_boundary0_i(cons::E_N)   = 1.5f * grid.k_b * n_n * T       + n_n * grid.m_n * phi_g_inner;
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
        const float p_i_0 = 2.0f * n_e * grid.k_b * T;
        const float p_n_0 =        n_n * grid.k_b * T;
        kInnerPiGhost = p_i_0 + (n_e * grid.m_i) * grid.g * ds0;
        kInnerPnGhost = p_n_0 + (n_n * grid.m_n) * grid.g * ds0;
    }

    // Height-dependent photoionization rate P_phot(h): Route B closure
    // (docs/photoionization_c7_inversion_plan.md; writeup §3.1).
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

        // Equilibrium inversion over EXACTLY the network Stage E will integrate,
        // so C7 is its fixed point under whatever rate flags are active:
        //   P_phot = (n_e/n_n)(α_r + α_c) n_e − (S_i + S_CR) n_e ,   clipped ≥ 0,
        // with the α_c (three-body) and S_i (Voronov) terms included only when
        // their Grid flags are set (default network drops both — see Stage E in
        // src/integrators.cpp and the flag docs in chromosphere.hpp).
        Vec a_tot  = recombination_rate_alpha(grid, Tc);          // α_r, always on
        if (grid.enable_threebody_recombination)
            a_tot += recombination_rate_alpha_c(grid, Tc, nec);   // α_c = κ_c n_e
        Vec S_tot  = ionization_rate_S_CR(grid, Tc);              // S_CR, always on
        if (grid.enable_direct_collisional_ionization)
            S_tot += ionization_rate_S(grid, Tc);                 // S_i, Voronov
        Vec P_inv = (nec / nnc) % (a_tot % nec) - S_tot % nec;
        P_inv = arma::clamp(P_inv, 0.0f, arma::datum::inf);

        // Frozen FAL-C rate for the NEQ region, and the tanh blend.
        const Vec P_falc = photoionization_rate_chae(h_cell);
        const float H_EQ_NEQ = 1500.0f;   // km, equilibrium ↔ NEQ crossover (Chae marker)
        const float BLEND_W  = 120.0f;    // km, blend width
        const Vec w_eq = 0.5f * (1.0f - arma::tanh((h_cell - H_EQ_NEQ) / BLEND_W));
        grid.photoionization_rate_i = w_eq % P_inv + (1.0f - w_eq) % P_falc;
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
        const float kappa0 = 1.0e-11f;   // W m^-1 K^-7/2
        const float T_cor  = 8.0e5f;     // K, quiet-Sun coronal apex
        const float L_cor  = 4.5e7f;     // m, coronal loop half-length (quiet-Sun)
        // ≈ 29 W m⁻². With the optically-thin TR radiative sink now present, this
        // sits in the stable window q ∈ [~30, 90] W m⁻²: the TR/corona radiates
        // the imposed flux and the chromosphere proper stays cool. Below ~10 W m⁻²
        // the weak heating drops the upper chromosphere into the radiative-loss-peak
        // thermal-instability (condensation) regime and the run fails; this value
        // keeps the TR as high as the stable window allows (base near 1.77×10³ km).
        grid.outer_heat_flux        = (2.0f / 7.0f) * kappa0
                                      * std::pow(T_cor, 3.5f) / L_cor;
        grid.impose_outer_heat_flux = grid.enable_radiative_cooling;
    }

    // Isotropic numerical diffusivity (Pandey 2024 §3.3), grid-scaled. Pandey's
    // benchmark is χ ≳ 30×10⁸ m² s⁻¹ at Δz = 40 km; here it only needs to damp the
    // base wave and any TR-gradient ringing (the coronal heating is carried by
    // q(T), not by χ), so a modest value suffices. Scale ∝ Δs so coarser grids
    // diffuse over a comparable number of cells.
    {
        const float ds_mean = arma::mean(grid.ds_i);   // ≈ 1.1×10⁴ m
        grid.numerical_diffusivity = 2.0e3f * ds_mean;  // ≈ 2.2×10⁷ m² s⁻¹
    }

    // TRAC (Johnston et al. 2020): broaden the under-resolved TR so the coarse
    // 100-cell grid recovers the grid-converged coronal density response while
    // conserving the TR-integrated radiation/heating (the Eq.-27 jump balance).
    // Paired with the radiative sink (q(T)/Λ(T) need a resolved TR to balance).
    grid.enable_trac   = grid.enable_radiative_cooling;
    grid.trac_T_chrom  = 2.0e4f;

    grid.broadcast();
    return xn;
}

void model_c7_update_bc(Grid& grid, const Vec& xn) {
    // Physically-grounded boundary closures (docs/boundary_conditions_plan.md):
    //   * Outer face — fixed lower-TR reservoir (RTV 1978 + Bradshaw & Emslie
    //     2020 + the single-fluid slaving of Gómez Míguez et al. 2024).
    //   * Inner face — discrete-HSE photospheric reservoir, V=U=0 (Pandey 2024).
    // apply_open_bcs lays down a default ghost layer; we overwrite both faces.
    apply_open_bcs(grid, xn);

    const auto  sz        = arma::size(grid.ns, num_of_eq);
    const float m_i       = grid.m_i;
    const float m_n       = grid.m_n;
    const float k_b       = grid.k_b;

    // ====================================================================
    // Outer face — lower transition region (h ≈ 2153 km, T ≈ 2.3×10⁴ K).
    // ====================================================================
    // The TR base is a constant-pressure layer (RTV 1978): pressure is continuous
    // across it while temperature jumps to the coronal value over the (here
    // unresolved) TR. We therefore do NOT impose a hot dense ghost wall — on the
    // coarse grid a fixed coronal T would over-conduct (κ_e∝T^{5/2}) and over-
    // pressurize the top cell, driving a spurious supersonic downflow. Instead:
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
    {
        const float phi_g_out = grid.phi_g_iph(grid.ns - 1);
        const float c_s   = std::sqrt(2.0f * grid.gamma_mono * k_b * kOuterTtr / m_i);
        const float V_cap = 0.05f * c_s;

        auto cap_ghost = [&](Vec& ob) {
            const float rho_i = ob(cons::RHO_I);
            const float rho_n = ob(cons::RHO_N);
            const float n_i   = rho_i / m_i;
            const float n_n   = rho_n / m_n;
            const float V_old = ob(cons::MOM_I) / rho_i;
            const float U_old = ob(cons::MOM_N) / rho_n;
            // Decode Neumann T_i, T_n from the energy apply_open_bcs packed.
            const float T_i_in = (ob(cons::E_I) - 0.5f * rho_i * V_old * V_old
                                  - rho_i * phi_g_out) / (3.0f * k_b * n_i);
            const float T_n_in = (ob(cons::E_N) - 0.5f * rho_n * U_old * U_old
                                  - rho_n * phi_g_out) / (1.5f * k_b * n_n);
            // Mach-0.05 cap for the quiet-Sun TR (Bradshaw & Emslie 2020). For
            // flare scenarios (grid.outer_free_outflow) the cap is lifted so a
            // supersonic evaporation upflow can leave the domain through a
            // transparent (Neumann) outer face — clamping it would choke the
            // ~100s km/s evaporated outflow.
            float V_g = V_old;
            if (!grid.outer_free_outflow) {
                if (V_g >  V_cap) V_g =  V_cap;
                if (V_g < -V_cap) V_g = -V_cap;
            }
            const float U_g = V_g;                // neutral locked to charge fluid
            ob(cons::MOM_I) = rho_i * V_g;
            ob(cons::MOM_N) = rho_n * U_g;
            ob(cons::E_I)   = 1.5f * k_b * n_i * 2.0f * T_i_in
                              + 0.5f * rho_i * V_g * V_g + rho_i * phi_g_out;
            ob(cons::E_N)   = 1.5f * k_b * n_n *        T_n_in
                              + 0.5f * rho_n * U_g * U_g + rho_n * phi_g_out;
        };
        cap_ghost(grid.outer_boundary0_i);
        cap_ghost(grid.outer_boundary1_i);
    }

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
        const float phi_g_inner = grid.phi_g_imh(0);
        const float rho_i_g = kInnerRhoIPinned;
        const float rho_n_g = kInnerRhoNPinned;

        Vec& ob0 = grid.inner_boundary0_i;
        ob0(cons::RHO_I) = rho_i_g;
        ob0(cons::RHO_N) = rho_n_g;
        ob0(cons::MOM_I) = 0.0f;                  // V = 0 reservoir
        ob0(cons::MOM_N) = 0.0f;                  // U = 0 reservoir
        ob0(cons::E_I)   = 1.5f * kInnerPiGhost + rho_i_g * phi_g_inner;
        ob0(cons::E_N)   = 1.5f * kInnerPnGhost + rho_n_g * phi_g_inner;
        grid.inner_boundary1_i = ob0;
    }

    grid.broadcast();
}

} // namespace chromosphere
