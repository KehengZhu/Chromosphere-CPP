/*!
 * @file grid.cpp
 * @brief Grid allocation, resize, and the packed/static-mesh cache rebuild.
 * @ingroup grid
 *
 * SHARED between both solvers. Grid::broadcast() fingerprints the five geometry
 * arrays it derives from (`ds_i`, `B_i`, `B_imh`, `B_iph`, `dinvB_ds_i`) and
 * rebuilds the caches only when one of them actually changed, so a per-step
 * boundary refresh that calls it costs an O(ns) comparison instead of ~50
 * `n_state`-sized temporaries.
 */
#include "chromosphere.hpp"

#include <stdexcept>

namespace chromosphere {

void Grid::init(arma::uword ns_in, float CFL_in) {
    ns              = ns_in;
    n_mixture_state = ns * num_of_mixture_eq;
    n_state         = ns * num_of_eq;
    CFL             = CFL_in;

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
    metrics_valid = false;

    mix_outer_boundary0.zeros(num_of_mixture_eq);
    mix_outer_boundary1.zeros(num_of_mixture_eq);
    mix_inner_boundary0.zeros(num_of_mixture_eq);
    mix_inner_boundary1.zeros(num_of_mixture_eq);

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
    eos_T_hint.assign(ns, std::numeric_limits<double>::quiet_NaN());
}

void Grid::resize(arma::uword ns_new) {
    if (ns_new == ns) return;
    ns              = ns_new;
    n_mixture_state = ns * num_of_mixture_eq;
    n_state         = ns * num_of_eq;

    ds_i.zeros(ns);
    B_imh.zeros(ns); B_iph.zeros(ns); B_i.zeros(ns);
    dinvB_ds_i.zeros(ns);
    phi_g_imh.zeros(ns);
    phi_g_iph.zeros(ns);

    s_i.zeros(ns);
    s_face.zeros(ns + 1);
    ds_iph_i.zeros(ns);
    ds_imh_i.zeros(ns);
    metrics_valid = false;

    // Ghost buffers keep their own solver's width — untouched. Physical constants, γ,
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
    eos_T_hint.assign(ns, std::numeric_limits<double>::quiet_NaN());
}

namespace {
// Exact (bit-for-bit) equality of a cached geometry snapshot against the live
// array. A mismatch in length counts as changed.
bool same_vector(const Vec& cached, const Vec& live) {
    if (cached.n_elem != live.n_elem) return false;
    for (arma::uword i = 0; i < live.n_elem; ++i)
        if (!(cached[i] == live[i])) return false;
    return true;
}
} // namespace

bool Grid::static_metrics_current() const {
    return metrics_valid
        && B_state.n_elem == n_state
        && s_face.n_elem == ns + 1
        && same_vector(metrics_ds_i, ds_i)
        && same_vector(metrics_B_i, B_i)
        && same_vector(metrics_B_imh, B_imh)
        && same_vector(metrics_B_iph, B_iph)
        && same_vector(metrics_dinvB_ds_i, dinvB_ds_i);
}

void Grid::broadcast() {
    if (static_metrics_current()) return;
    force_rebuild_metrics();
}

void Grid::force_rebuild_metrics() {
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

    // Packed (n_state) broadcasts. Previously built as
    //   dst.zeros(n_state); for k: dst += scalar_to(*this, src, k);
    // which allocated num_of_eq n_state-sized temporaries per field. The direct
    // fill below writes the same values; the explicit `0.0f +` reproduces the
    // accumulation's one signed-zero difference (0.0f + -0.0f == +0.0f) so the
    // result is bit-for-bit what the old loop produced.
    auto pack = [&](Vec& dst, const Vec& src) {
        dst.set_size(n_state);
        for (arma::uword k = 0; k < num_of_eq; ++k) {
            Vec::elem_type* out = dst.memptr() + k * ns;
            for (arma::uword i = 0; i < ns; ++i) out[i] = src[i];
        }
    };
    pack(B_state_imh,    B_imh);
    pack(B_state_iph,    B_iph);
    pack(B_state,        B_i);
    pack(ds_state,       ds_i);
    pack(dinvB_ds_state, dinvB_ds_i);
    pack(ds_iph_state,   ds_iph_i);
    pack(ds_imh_state,   ds_imh_i);

    metrics_ds_i = ds_i;
    metrics_B_i = B_i;
    metrics_B_imh = B_imh;
    metrics_B_iph = B_iph;
    metrics_dinvB_ds_i = dinvB_ds_i;
    metrics_valid = true;
    ++metrics_generation_;
}

} // namespace chromosphere
