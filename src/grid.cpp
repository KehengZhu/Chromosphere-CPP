#include "chromosphere.hpp"

namespace chromosphere {

void Grid::init(arma::uword ns_in, float CFL_in) {
    ns      = ns_in;
    n_state = ns * num_of_eq;
    CFL     = CFL_in;

    // Physical constants are already set by the default member initializers
    // in the header. (Re-stated here for clarity.)
    gamma_mono = 5.0f / 3.0f;
    m_i        = 1.6726219e-27f;
    m_n        = m_i;
    m_e        = 9.10938356e-31f;
    g          = 0.27395e3f;
    mu_0       = 4.0f * static_cast<float>(arma::datum::pi) * 1.0e-7f;
    k_b        = 1.380649e-23f;
    q_e        = 1.602176634e-19f;

    ds_i.zeros(ns);
    B_imh.zeros(ns); B_iph.zeros(ns); B_i.zeros(ns);
    dinvB_ds_i.zeros(ns);
    phi_g_imh.zeros(ns);
    phi_g_iph.zeros(ns);

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
}

void Grid::broadcast() {
    B_state_imh.zeros(n_state);
    B_state_iph.zeros(n_state);
    B_state.zeros(n_state);
    ds_state.zeros(n_state);
    dinvB_ds_state.zeros(n_state);
    for (arma::uword k = 0; k < num_of_eq; ++k) {
        B_state_imh    += scalar_to(*this, B_imh,      k);
        B_state_iph    += scalar_to(*this, B_iph,      k);
        B_state        += scalar_to(*this, B_i,        k);
        ds_state       += scalar_to(*this, ds_i,       k);
        dinvB_ds_state += scalar_to(*this, dinvB_ds_i, k);
    }
}

} // namespace chromosphere
