#include "chromosphere.hpp"

#include <stdexcept>

namespace chromosphere {

void GammaConductionScratch::resize(std::size_t n) {
    rho.resize(n); e_old.resize(n); temperature.resize(n); target.resize(n);
    conductivity.resize(n); capacity.resize(n); n_e.resize(n); n_hi.resize(n);
    g_left.resize(n); g_right.resize(n);
    a.resize(n); b.resize(n); c.resize(n); rhs.resize(n); delta.resize(n);
}

void Grid::init(arma::uword ns_in, float CFL_in) {
    ns      = ns_in;
    n_state = ns * num_of_eq;
    CFL     = CFL_in;

    // Physical constants are already set by the default member initializers
    // in the header. (Re-stated here for clarity.)
    gamma_mono = 5.0f / 3.0f;
    m_i        = static_cast<float>(eos_constants::m_h);
    m_n        = m_i;
    m_e        = static_cast<float>(eos_constants::m_e);
    g          = 0.27395e3f;
    mu_0       = 4.0f * static_cast<float>(arma::datum::pi) * 1.0e-7f;
    k_b        = static_cast<float>(eos_constants::k_b);
    q_e        = 1.602176634e-19f;
    chi_H_J    = static_cast<float>(eos_constants::chi_h);

    ds_i.zeros(ns);
    B_imh.zeros(ns); B_iph.zeros(ns); B_i.zeros(ns);
    dinvB_ds_i.zeros(ns);
    phi_g_imh.zeros(ns);
    phi_g_iph.zeros(ns);

    s_i.zeros(ns);
    s_face.zeros(ns + 1);
    ds_iph_i.zeros(ns);
    ds_imh_i.zeros(ns);
    uniform_mesh = true;

    outer_boundary0_i.zeros(num_of_eq);
    outer_boundary1_i.zeros(num_of_eq);
    inner_boundary0_i.zeros(num_of_eq);
    inner_boundary1_i.zeros(num_of_eq);

    B_state_imh.zeros(n_state);
    B_state_iph.zeros(n_state);
    B_state.zeros(n_state);
    ds_state.zeros(n_state);
    dinvB_ds_state.zeros(n_state);
    dt_state.zeros(n_state);
    ds_iph_state.zeros(n_state);
    ds_imh_state.zeros(n_state);
}

void Grid::resize(arma::uword ns_new) {
    if (ns_new == ns) return;
    ns      = ns_new;
    n_state = ns * num_of_eq;

    ds_i.zeros(ns);
    B_imh.zeros(ns); B_iph.zeros(ns); B_i.zeros(ns);
    dinvB_ds_i.zeros(ns);
    phi_g_imh.zeros(ns);
    phi_g_iph.zeros(ns);

    s_i.zeros(ns);
    s_face.zeros(ns + 1);
    ds_iph_i.zeros(ns);
    ds_imh_i.zeros(ns);

    // Ghost buffers stay length num_of_eq — untouched. Physical constants, γ,
    // CFL and every runtime toggle are deliberately preserved (this is a pure
    // reallocation of the ns-sized fields), so a scenario can size the coarse
    // grid, set its flags, then resize to the refined count.
    B_state_imh.zeros(n_state);
    B_state_iph.zeros(n_state);
    B_state.zeros(n_state);
    ds_state.zeros(n_state);
    dinvB_ds_state.zeros(n_state);
    dt_state.zeros(n_state);
    ds_iph_state.zeros(n_state);
    ds_imh_state.zeros(n_state);
}

void Grid::broadcast() {
    // --- static-mesh metric caches --------------------------------------------
    // Rebuild the center-to-center distances and canonical coordinates from ds_i,
    // taking cell centers at face midpoints. Mirror the boundary cell width into
    // the ghost (Neumann) at both ends so ds_iph_i[ns-1]/ds_imh_i[0] reproduce the
    // legacy 0.5*(ds_i + ip1/im1(ds_i,SLICE)) expressions used by rhs/conduction.
    if (ds_i.n_elem == ns && ns > 0) {
        const float dmin = ds_i.min();
        const float dmax = ds_i.max();
        // Reject a genuinely broken mesh. The all-zero state right after init()
        // (dmax == 0) is treated as "not yet populated" and skipped, not rejected.
        if (dmax > 0.0f) {
            if (!ds_i.is_finite() || dmin <= 0.0f)
                throw std::runtime_error(
                    "Grid::broadcast: ds_i must be strictly positive and finite");
            // Relative spread tolerance: 1 part in 10^6 counts as uniform.
            uniform_mesh = (dmax - dmin) <= 1.0e-6f * dmax;

            s_face(0) = 0.0f;
            for (arma::uword i = 0; i < ns; ++i)
                s_face(i + 1) = s_face(i) + ds_i(i);
            for (arma::uword i = 0; i < ns; ++i) {
                s_i(i) = 0.5f * (s_face(i) + s_face(i + 1));
                const float ds_ip1 = (i + 1 < ns) ? ds_i(i + 1) : ds_i(i);   // mirror
                const float ds_im1 = (i > 0)      ? ds_i(i - 1) : ds_i(i);   // mirror
                ds_iph_i(i) = 0.5f * (ds_i(i) + ds_ip1);
                ds_imh_i(i) = 0.5f * (ds_im1 + ds_i(i));
            }
        }
    }

    B_state_imh.zeros(n_state);
    B_state_iph.zeros(n_state);
    B_state.zeros(n_state);
    ds_state.zeros(n_state);
    dinvB_ds_state.zeros(n_state);
    ds_iph_state.zeros(n_state);
    ds_imh_state.zeros(n_state);
    for (arma::uword k = 0; k < num_of_eq; ++k) {
        B_state_imh    += scalar_to(*this, B_imh,      k);
        B_state_iph    += scalar_to(*this, B_iph,      k);
        B_state        += scalar_to(*this, B_i,        k);
        ds_state       += scalar_to(*this, ds_i,       k);
        dinvB_ds_state += scalar_to(*this, dinvB_ds_i, k);
        ds_iph_state   += scalar_to(*this, ds_iph_i,   k);
        ds_imh_state   += scalar_to(*this, ds_imh_i,   k);
    }
}

} // namespace chromosphere
