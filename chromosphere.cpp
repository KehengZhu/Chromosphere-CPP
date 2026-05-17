#include "chromosphere.hpp"
#include "armadillo"

// By Keheng: Data structure here is different from the one in the dissertation.
// Here xn is a pseudo 3D array.
// The order is: ns, num_of_eq

namespace chromosphere {

uword ns;       // Number of grid points in the x-direction
uword num_of_elem; // Total number of elements (ns * num_of_eq)
// constants
float gammamono; // Ratio of specific heats (adiabatic index for a monatomic ideal gas)
float alpha_p;   // Divergence error propagation parameter (alpha_p)
float pi;        // The value of π (pi)

// Physical constants
float m_i;       // Mass of proton in kilograms (also m0, normalization constant for mass)
float m_n;       // Normalized mass of neutrals
float m_e;       // Normalized mass of electron
float g;         // Gravitational acceleration in normalized units, default is Sun value of 0.27395 km/s^2
float mu_0;      // Permeability of free space (magnetic constant)
float k_b;       // Boltzmann constant in SI units
float e_;        // Elementary charge in Coulombs

Vec ds_i;          // Grid spacing in the x-direction in normalized domain units
float CFL;       // Courant-Friedrichs-Lewy number

// Stream-aligned quantities
Vec cellX_i, cellY_i, cellZ_i;
Vec dinvB_ds_i, B_imh, B_iph, B_i;
Vec gPotential_imh, gPotential_iph;
Vec outer_boundary0_i, outer_boundary1_i;
Vec inner_boundary0_i, inner_boundary1_i;

// Broadcasted variables
Vec B_iimh, B_iiph, B_ii, dinvB_ds_ii, dt_ii, ds_ii;

//--------------------------------------------------------------------
// Physical functions

// Flux limiter: see (21)
Vec flux_lim(const Vec& r) {
    // minmod
    Vec one(arma::size(r), fill::ones);
    Vec zero(arma::size(r), fill::zeros);
    Vec res = arma::max(zero, arma::min(one, r));
    // if there is nan, the flux limiter value does not matter.
    res(arma::find_nan(res)).zeros();
    return res;
}

// Get the flux from conserved variables
Vec cal_F_ii(const Vec& xn_ii) {
    Vec F_ii(arma::size(xn_ii), fill::zeros);
    const Vec rhoi = get_scalar(xn_ii, CNI);
    const Vec rhon = get_scalar(xn_ii, CNN);
    const Vec rhov = get_scalar(xn_ii, CNV);
    const Vec rhou = get_scalar(xn_ii, CNU);
    const Vec ei = get_scalar(xn_ii, CEI);
    const Vec en = get_scalar(xn_ii, CEN);

    const Vec ni = rhoi/m_i;
    const Vec nn = rhon/m_n;

    const Vec vv = rhov / rhoi;
    const Vec uu = rhou / rhon;
    const Vec phig = 0.5*(gPotential_imh+gPotential_iph);
    const Vec pi = 2.0/3.0*ei - 1.0/3.0*rhoi%vv%vv - 2.0/3.0*rhoi%phig; 
    const Vec pn = 2.0/3.0*en - 1.0/3.0*rhon%uu%uu - 2.0/3.0*rhon%phig;

    F_ii += scalar_to(rhov, CNI);
    F_ii += scalar_to(rhou, CNN);
    F_ii += scalar_to(rhoi%vv%vv + pi, CNV);
    F_ii += scalar_to(rhon%uu%uu + pn, CNU);
    // cout << "Fu = " << rhon%uu%uu + pn << endl;
    F_ii += scalar_to((ei + pi)%vv, CEI);
    F_ii += scalar_to((en + pn)%uu, CEN);
    return F_ii;
}

// Find the max eigenvalue
Vec find_spectral_radius_ii(const Vec& xn_ii) {
    Vec res_ii(size(xn_ii), fill::zeros);
    const Vec rhoi = get_scalar(xn_ii, CNI);
    const Vec rhon = get_scalar(xn_ii, CNN);
    const Vec rhov = get_scalar(xn_ii, CNV);
    const Vec rhou = get_scalar(xn_ii, CNU);
    const Vec ei = get_scalar(xn_ii, CEI);
    const Vec en = get_scalar(xn_ii, CEN);

    const Vec ni = rhoi/m_i;
    const Vec nn = rhon/m_n;

    const Vec vv = rhov / rhoi;
    const Vec uu = rhou / rhon;

    const Vec phig = 0.5*(gPotential_imh+gPotential_iph);
    const Vec pi = 2.0/3.0*ei - 1.0/3.0*rhoi%vv%vv - 2.0/3.0*rhoi%phig; 
    const Vec pn = 2.0/3.0*en - 1.0/3.0*rhon%uu%uu - 2.0/3.0*rhon%phig;
    const Vec Ti = pi/(2.0*ni*k_b);
    const Vec Tn = pn/(nn*k_b);

    const Vec csi = arma::sqrt(gammamono*pi/rhoi);
    const Vec csn = arma::sqrt(gammamono*pn/rhon);

    Vec spc_r_i = arma::max(arma::abs(vv)+arma::abs(csi), 
                          arma::abs(uu)+arma::abs(csn));
    
    // Make spc_r a matrix
    for(int i=0; i<num_of_eq; i++) {
        res_ii += scalar_to(spc_r_i, i);
    }
    
    return res_ii;
}

// Explicit source: pressure-area (geometric) and gravity, per writeup §4.3.
// S_E = [0, 0, p_i B D - rho_i G, p_n B D - rho_n G, 0, 0]^T
// with D = d(1/B)/ds, G = d phi_g/ds.
Vec cal_S_ii(const Vec& xn_ii) {
    Vec S_ii(arma::size(xn_ii), arma::fill::zeros);
    const Vec rhoi = get_scalar(xn_ii, CNI);
    const Vec rhon = get_scalar(xn_ii, CNN);
    const Vec rhov = get_scalar(xn_ii, CNV);
    const Vec rhou = get_scalar(xn_ii, CNU);
    const Vec ei = get_scalar(xn_ii, CEI);
    const Vec en = get_scalar(xn_ii, CEN);

    const Vec vv = rhov / rhoi;
    const Vec uu = rhou / rhon;
    const Vec phig = 0.5*(gPotential_imh+gPotential_iph);
    const Vec pi = 2.0/3.0*ei - 1.0/3.0*rhoi%vv%vv - 2.0/3.0*rhoi%phig;
    const Vec pn = 2.0/3.0*en - 1.0/3.0*rhon%uu%uu - 2.0/3.0*rhon%phig;

    // Gravity: g = -d phi_g/ds evaluated cell-centered from the face potentials.
    const Vec gb = -(gPotential_iph-gPotential_imh)/ds_i;

    S_ii += scalar_to(pi%B_i%dinvB_ds_i + rhoi%gb, CNV);
    S_ii += scalar_to(pn%B_i%dinvB_ds_i + rhon%gb, CNU);
    return S_ii;
}

// By Keheng: this function is the same as Re_MUSCL in the CoMFi code.
// It is only 1.5D though, so z-dir is not considered.
Vec rhs_explicit_ii(const Vec& xn_ii) {
    Vec xn_rhs_ii(size(xn_ii));

// BEGIN GET UL UR ---------------------------------
    // innerbound and outerbound are given as conserved variables.
    Vec cons_xn_iip1 = ip1(xn_ii);
    Vec cons_xn_iim1 = im1(xn_ii);
    Vec cons_xn_iip2 = ip2(xn_ii);
    Vec cons_xn_iim2 = im2(xn_ii);

// Should use primitive variables to limit the slope
    Vec prim_xn_ii = cons2prim(xn_ii);
    Vec prim_xn_iip1 = cons2prim(cons_xn_iip1);
    Vec prim_xn_iim1 = cons2prim(cons_xn_iim1);

    Vec prim_xn_iip2 = cons2prim(cons_xn_iip2);
    Vec prim_xn_iim2 = cons2prim(cons_xn_iim2);

    const Vec dxn_iiph = prim_xn_iip1 - prim_xn_ii;
    const Vec dxn_iimh = prim_xn_ii - prim_xn_iim1;

    // may need a CON_stop function later
    if ((dxn_iiph.has_nan() + dxn_iimh.has_nan()) > 0) {
        cout << "XN" << endl;
        print_xn(xn_ii);
        cout << "XNPRIM" << endl;
        print_xn(prim_xn_ii);
        cout << "dxn_iph has nan: " << dxn_iiph.has_nan() << endl
             << "dxn_imh has nan: " << dxn_iimh.has_nan() << endl;
        exit(1);
    }

    // limiter inputs
    Vec r_ii = dxn_iimh / dxn_iiph;
    Vec r_iip1 = dxn_iiph / (prim_xn_iip2 - prim_xn_iip1);
    Vec r_iim1 = (prim_xn_iim1 - prim_xn_iim2) / dxn_iimh;

    // extrapolated cell edge variables
    // rewrite according to MUSCL paper eq(4.5)
    Vec Rxn_iiph = prim_xn_iip1 - 0.5*flux_lim(r_iip1)%(prim_xn_iip2 - prim_xn_iip1);
    Vec Lxn_iiph = prim_xn_ii     + 0.5*flux_lim(r_ii)%(prim_xn_iip1 - prim_xn_ii);
    Vec Rxn_iimh = prim_xn_ii     - 0.5*flux_lim(r_ii)%(prim_xn_iip1 - prim_xn_ii);
    Vec Lxn_iimh = prim_xn_iim1 + 0.5*flux_lim(r_iim1)%(prim_xn_ii - prim_xn_iim1);

    if ((Lxn_iiph.has_nan() + Rxn_iiph.has_nan()) + (Lxn_iimh.has_nan() + Rxn_iimh.has_nan()) > 0) {
        cout << "Lxn_iph has nan: " << Lxn_iiph.has_nan() << endl
             << "Rxn_iph has nan: " << Rxn_iiph.has_nan() << endl
             << "Lxn_imh has nan: " << Lxn_iimh.has_nan() << endl
             << "Rxn_imh has nan: " << Rxn_iimh.has_nan() << endl
             << "Lxn_iph(nn) zeros: " << endl;
        exit(1);
    }

    Rxn_iiph = prim2cons(Rxn_iiph);
    Lxn_iiph = prim2cons(Lxn_iiph);
    Rxn_iimh = prim2cons(Rxn_iimh);
    Lxn_iimh = prim2cons(Lxn_iimh);

    // Prediction for next step.
    Vec prim_xt_ii = cons2prim(xn_ii - dt_ii/ds_ii%(cal_F_ii(Lxn_iiph) - cal_F_ii(Rxn_iimh)));
    Vec prim_xt_iip1 = cons2prim(ip1(prim2cons(prim_xt_ii)));
    Vec prim_xt_iim1 = cons2prim(im1(prim2cons(prim_xt_ii)));

    // cout << "dF/ds(uu) = " << get_scalar(dt_ii/ds_ii%(cal_F_ii(Lxn_iiph) - cal_F_ii(Rxn_iimh)), CNU) << endl;

    // 2nd time accuracy reconstruction
    Rxn_iiph = 0.5*(prim_xn_iip1+prim_xt_iip1) - 0.5*flux_lim(r_iip1)%(prim_xn_iip2 - prim_xn_iip1);
    Lxn_iiph = 0.5*(prim_xn_ii+prim_xt_ii) + 0.5*flux_lim(r_ii)%(prim_xn_iip1 - prim_xn_ii);
    Rxn_iimh = 0.5*(prim_xn_ii+prim_xt_ii) - 0.5*flux_lim(r_ii)%(prim_xn_iip1 - prim_xn_ii);
    Lxn_iimh = 0.5*(prim_xn_iim1+prim_xt_iim1) + 0.5*flux_lim(r_iim1)%(prim_xn_ii - prim_xn_iim1);

    Rxn_iiph = prim2cons(Rxn_iiph);
    Lxn_iiph = prim2cons(Lxn_iiph);
    Rxn_iimh = prim2cons(Rxn_iimh);
    Lxn_iimh = prim2cons(Lxn_iimh);

// END GET UL UR -----------------------------

    const Vec a_iimh =
        arma::max(find_spectral_radius_ii(Lxn_iimh), find_spectral_radius_ii(Rxn_iimh));
    const Vec a_iiph =
        arma::max(find_spectral_radius_ii(Lxn_iiph), find_spectral_radius_ii(Rxn_iiph));

    const Vec F_iimh =
        0.5 * (cal_F_ii(Lxn_iimh) + cal_F_ii(Rxn_iimh) - a_iimh % (Rxn_iimh - Lxn_iimh));
    const Vec F_iiph =
        0.5 * (cal_F_ii(Lxn_iiph) + cal_F_ii(Rxn_iiph) - a_iiph % (Rxn_iiph - Lxn_iiph));

    if ((a_iimh.has_nan() + a_iiph.has_nan() + F_iimh.has_nan() + F_iiph.has_nan()) > 0) {
        cout << "a_imh has nan: " << a_iimh.has_nan() << endl
             << "a_iph has nan: " << a_iiph.has_nan() << endl
             << "Fximh has nan: " << F_iimh.has_nan() << endl
             << "Fxiph has nan: " << F_iiph.has_nan() << endl;
        exit(1);
    }

    xn_rhs_ii = -1.0*B_ii%(F_iiph/B_iiph - F_iimh/B_iimh)/ds_ii + cal_S_ii(xn_ii);

    return xn_rhs_ii;
}

// Implicit terms (writeup §4.3, R_I): ion-neutral drag, collisional heating,
// frictional heating, and heat conduction. Pressure-area and gravity sit in
// the explicit source cal_S_ii.
Vec rhs_implicit_ii(const Vec& xn_ii) {
    Vec RI_ii(arma::size(xn_ii), fill::zeros);

    // Cell-centered primitives.
    const Vec rhoi = get_scalar(xn_ii, CNI);
    const Vec rhon = get_scalar(xn_ii, CNN);
    const Vec rhov = get_scalar(xn_ii, CNV);
    const Vec rhou = get_scalar(xn_ii, CNU);
    const Vec ei   = get_scalar(xn_ii, CEI);
    const Vec en   = get_scalar(xn_ii, CEN);

    const Vec ni = rhoi/m_i;
    const Vec nn = rhon/m_n;
    const Vec vv = rhov / rhoi;
    const Vec uu = rhou / rhon;
    const Vec phig = 0.5*(gPotential_imh+gPotential_iph);
    const Vec pi = 2.0/3.0*ei - 1.0/3.0*rhoi%vv%vv - 2.0/3.0*rhoi%phig;
    const Vec pn = 2.0/3.0*en - 1.0/3.0*rhon%uu%uu - 2.0/3.0*rhon%phig;
    const Vec Ti = pi/(2.0*ni*k_b);
    const Vec Tn = pn/(nn*k_b);

    const Vec Ke = kappa_e(ni, nn, Ti);
    const Vec Ki(size(Ti), arma::fill::zeros);
    const Vec Kn = kappa_n(ni, nn, Ti, Tn);

    // i+1 neighbor (uses outer ghost at i = ns-1).
    const Vec xn_iip1 = ip1(xn_ii);
    const Vec rhoi_ip1 = get_scalar(xn_iip1, CNI);
    const Vec rhon_ip1 = get_scalar(xn_iip1, CNN);
    const Vec rhov_ip1 = get_scalar(xn_iip1, CNV);
    const Vec rhou_ip1 = get_scalar(xn_iip1, CNU);
    const Vec ei_ip1   = get_scalar(xn_iip1, CEI);
    const Vec en_ip1   = get_scalar(xn_iip1, CEN);

    const Vec ni_ip1 = rhoi_ip1/m_i;
    const Vec nn_ip1 = rhon_ip1/m_n;
    const Vec vv_ip1 = rhov_ip1 / rhoi_ip1;
    const Vec uu_ip1 = rhou_ip1 / rhon_ip1;
    Vec phig_ip1 = ip1(0.5*(gPotential_iph+gPotential_imh), SLICE);
    phig_ip1[ns-1] = gPotential_iph[ns-1];
    const Vec pi_ip1 = 2.0/3.0*ei_ip1 - 1.0/3.0*rhoi_ip1%vv_ip1%vv_ip1 - 2.0/3.0*rhoi_ip1%phig_ip1;
    const Vec pn_ip1 = 2.0/3.0*en_ip1 - 1.0/3.0*rhon_ip1%uu_ip1%uu_ip1 - 2.0/3.0*rhon_ip1%phig_ip1;
    const Vec Ti_ip1 = pi_ip1/(2.0*ni_ip1*k_b);
    const Vec Tn_ip1 = pn_ip1/(nn_ip1*k_b);

    const Vec Ke_ip1 = kappa_e(ni_ip1, nn_ip1, Ti_ip1);
    const Vec Ki_ip1(size(Ti), arma::fill::zeros);
    const Vec Kn_ip1 = kappa_n(ni_ip1, nn_ip1, Ti_ip1, Tn_ip1);

    // i-1 neighbor (uses inner ghost at i = 0).
    const Vec xn_iim1 = im1(xn_ii);
    const Vec rhoi_im1 = get_scalar(xn_iim1, CNI);
    const Vec rhon_im1 = get_scalar(xn_iim1, CNN);
    const Vec rhov_im1 = get_scalar(xn_iim1, CNV);
    const Vec rhou_im1 = get_scalar(xn_iim1, CNU);
    const Vec ei_im1   = get_scalar(xn_iim1, CEI);
    const Vec en_im1   = get_scalar(xn_iim1, CEN);

    const Vec ni_im1 = rhoi_im1/m_i;
    const Vec nn_im1 = rhon_im1/m_n;
    const Vec vv_im1 = rhov_im1 / rhoi_im1;
    const Vec uu_im1 = rhou_im1 / rhon_im1;
    Vec phig_im1 = im1(0.5*(gPotential_iph+gPotential_imh), SLICE);
    phig_im1[0] = gPotential_imh[0];
    const Vec pi_im1 = 2.0/3.0*ei_im1 - 1.0/3.0*rhoi_im1%vv_im1%vv_im1 - 2.0/3.0*rhoi_im1%phig_im1;
    const Vec pn_im1 = 2.0/3.0*en_im1 - 1.0/3.0*rhon_im1%uu_im1%uu_im1 - 2.0/3.0*rhon_im1%phig_im1;
    const Vec Ti_im1 = pi_im1/(2.0*ni_im1*k_b);
    const Vec Tn_im1 = pn_im1/(nn_im1*k_b);

    const Vec Ke_im1 = kappa_e(ni_im1, nn_im1, Ti_im1);
    const Vec Ki_im1(size(Ti), arma::fill::zeros);
    const Vec Kn_im1 = kappa_n(ni_im1, nn_im1, Ti_im1, Tn_im1);

    // --- ion-neutral drag (alpha = rho_i * nu_in) ----------------------
    const Vec nuin = nu_in(nn, Ti, Tn);
    const Vec alpha = rhoi % nuin;
    const Vec w = vv - uu;
    RI_ii += scalar_to(-alpha%w, CNV);
    RI_ii += scalar_to( alpha%w, CNU);

    // --- collisional + frictional heating (writeup eq R_I energy rows) -
    RI_ii += scalar_to(alpha/(m_i+m_n) % (3.0*k_b*(Tn-Ti) + m_n*w%w) - alpha%vv%w, CEI);
    RI_ii += scalar_to(alpha/(m_i+m_n) % (3.0*k_b*(Ti-Tn) + m_i*w%w) + alpha%uu%w, CEN);

    // --- conservative field-aligned heat conduction (writeup eq 4.4):
    //   C = (B/ds) [ K_{i+1/2}/B_{i+1/2} (T_{i+1}-T_i)/ds
    //              - K_{i-1/2}/B_{i-1/2} (T_i-T_{i-1})/ds ]
    // Face conductivities by arithmetic average.
    const Vec Kt_iph = 0.5 * ((Ke + Ki) + (Ke_ip1 + Ki_ip1));
    const Vec Kt_imh = 0.5 * ((Ke + Ki) + (Ke_im1 + Ki_im1));
    const Vec Kn_iph = 0.5 * (Kn + Kn_ip1);
    const Vec Kn_imh = 0.5 * (Kn + Kn_im1);

    const Vec q_i_iph = Kt_iph / B_iph % (Ti_ip1 - Ti)     / ds_i;
    const Vec q_i_imh = Kt_imh / B_imh % (Ti     - Ti_im1) / ds_i;
    const Vec q_n_iph = Kn_iph / B_iph % (Tn_ip1 - Tn)     / ds_i;
    const Vec q_n_imh = Kn_imh / B_imh % (Tn     - Tn_im1) / ds_i;

    const Vec C_i = B_i % (q_i_iph - q_i_imh) / ds_i;
    const Vec C_n = B_i % (q_n_iph - q_n_imh) / ds_i;
    RI_ii += scalar_to(C_i, CEI);
    RI_ii += scalar_to(C_n, CEN);

    return RI_ii;
}


// calculate max(V+c_fast) according to eq(23) in the paper.
// Returns an array of ns elements
Vec get_max_v_i(const Vec& xn_ii) {
    Vec spc_r_ii = find_spectral_radius_ii(xn_ii);
    return get_scalar(spc_r_ii, 0);
}

Vec cal_dt_i(const Vec& xn_ii) {
    // return CFL*ds_i/get_max_v_i(xn_ii);
    Vec dt_i = CFL*ds_i/get_max_v_i(xn_ii);
    float dt_min_i = arma::min(dt_i);
    return dt_min_i*arma::ones<Vec>(ns);
}

Vec advance_Euler_ii(const Vec& xn_ii, const Vec& dt_i) {
    const float tol = 1e-6;
    // Initialize dt_ii before all functions.
    dt_ii.zeros(num_of_elem);
    for(int i=0; i<num_of_eq; i++) {
        dt_ii += scalar_to(dt_i, i);
    }
    Vec RE_n_ii = rhs_explicit_ii(xn_ii);
    Vec RI_np1_kp1_ii(arma::size(xn_ii), arma::fill::zeros);
    Vec RI_np1_k_ii(arma::size(xn_ii), arma::fill::zeros);
    Vec xnp1_ii;
    Vec error, norm;
    uword it_cnt = 0;

    // Solve implicit terms using fix-point iteration
    do {
        RI_np1_k_ii = RI_np1_kp1_ii;
        xnp1_ii = xn_ii + dt_ii%(RE_n_ii + RI_np1_k_ii);
        // Picard iteration on R_I (writeup §3.7) diverges whenever
        // dt * |dR_I/du| > 1 — i.e., when ion-neutral drag (alpha*dt ~ 60 at
        // Model C7 conditions) or heat conduction makes R_I stiff. A true
        // implicit solve (Newton with the R_I Jacobian) is required to enable
        // this branch. Until then, leave R_I at zero so the step is pure
        // explicit Euler.
        // RI_np1_kp1_ii = rhs_implicit_ii(xnp1_ii);
        RI_np1_kp1_ii.zeros(arma::size(xn_ii));
        error = (RI_np1_kp1_ii - RI_np1_k_ii)%(RI_np1_kp1_ii - RI_np1_k_ii);
        norm = RI_np1_k_ii%RI_np1_k_ii;
        it_cnt++;
    } while(arma::sum(error)/arma::sum(norm) > tol && it_cnt < 25);
    // cout << "----------------------------------  it_cnt = " << it_cnt << endl;

    return xnp1_ii;
}

Vec advance_RK4(const Vec& xn_ii, const Vec& dt_i) {
    dt_ii.zeros(num_of_elem);
    for(int i=0; i<num_of_eq; i++) {
        dt_ii += scalar_to(dt_i, i);
    }
    const Vec k1 = dt_ii%rhs_explicit_ii(xn_ii);
    const Vec k2 = dt_ii%rhs_explicit_ii(xn_ii + 0.5*k1);
    const Vec k3 = dt_ii%rhs_explicit_ii(xn_ii + 0.5*k2);
    const Vec k4 = dt_ii%rhs_explicit_ii(xn_ii + k3);
    return xn_ii + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
}

void print_xn(const Vec& xn) {
    for(int i=0; i<ns; i++) {
        for(int j=0; j<num_of_eq; j++) {
            cout << xn(sub2ind(size(ns, num_of_eq),i, j)) << " ";
        }
        cout << endl;
    }
}

} // end namespace chromosphere

