#include "../chromosphere.hpp"
#include "armadillo"

//--------------------------------------------------------------------
// C/Fortran Interface
// The interfaces take in float arrays, convert them to Armadillo vectors,
// call the corresponding C++ functions, and convert the results back to float arrays.

struct Chromo_params {
    // SOLUTION INDEX
    // Conserved Variable Indices
    int CNI;        // Ion density
    int CNN;        // Neutral density
    int CNV;        // Ion velocity
    int CNU;        // Neutral velocity
    int CEI;        // Ion energy
    int CEN;        // Neutral energy
    // Primitive Variable Indices
    int PNI;        // Ion density
    int PNN;        // Neutral density
    int PV ;        // Ion velocity
    int PU ;        // Neutral velosity
    int PPI;        // Ion pressure
    int PPN;        // Neutral pressure

    int ns;       // Number of grid points in the x-direction
    int num_of_eq; // Number of equations (or number of unknowns)
    int num_of_elem; // Total number of elements (nx * nz * num_of_eq)

    // constants
    float gammamono; // Ratio of specific heats (adiabatic index for a monatomic ideal gas)
    float alpha_p;   // Divergence error propagation parameter (alpha_p)
    float pi;        // The value of π (pi)

    // Physical constants
    float m_i;       // Mass of proton in kilograms (also m0, normalization constant for mass)
    float m_n;       // mass of neutrals
    float m_e;       // mass of electron
    float g;         // Gravitational acceleration in normalized units, default is Sun value of 0.27395 km/s^2
    float mu_0;      // Permeability of free space (magnetic constant)
    float k_b;       // Boltzmann constant in SI units
    float e_;        // Elementary charge in Coulombs

    float* ds;       // Grid spacing in the x-direction in normalized domain units (array)
    float CFL;       // Courant-Friedrichs-Lewy condition

    // stream-aligned quantities (arrays)
    int is_allocated;
    float* dt;
    float* xn;
    float* cellX;
    float* cellY;
    float* cellZ;
    float* B_imh;
    float* B_iph;
    float* B_i;
    float* dinvB_ds;
    float* gPotential_imh;
    float* gPotential_iph;
    float* outer_bound0;
    float* outer_bound1;
    float* inner_bound0;
    float* inner_bound1;
};

void params2chromo(const struct Chromo_params* params) {
    // Copy params from struct back to chromosphere variables
    chromosphere::ns = params->ns;
    chromosphere::num_of_elem = params->num_of_elem;

    chromosphere::gammamono = params->gammamono;
    chromosphere::alpha_p = params->alpha_p;
    chromosphere::pi = params->pi;

    chromosphere::m_i = params->m_i;
    chromosphere::m_n = params->m_n;
    chromosphere::m_e = params->m_e;
    chromosphere::g = params->g;
    chromosphere::mu_0 = params->mu_0;
    chromosphere::k_b = params->k_b;
    chromosphere::e_ = params->e_;

    chromosphere::CFL = params->CFL;

    // Allocate memory
    chromosphere::ds_i.zeros(params->ns);
    chromosphere::cellX_i.zeros(params->ns);
    chromosphere::cellY_i.zeros(params->ns);
    chromosphere::cellZ_i.zeros(params->ns);
    chromosphere::B_imh.zeros(params->ns);
    chromosphere::B_iph.zeros(params->ns);
    chromosphere::B_i.zeros(params->ns);
    chromosphere::dinvB_ds_i.zeros(params->ns);
    chromosphere::gPotential_imh.zeros(params->ns);
    chromosphere::gPotential_iph.zeros(params->ns);
    chromosphere::outer_boundary0_i.zeros(params->num_of_eq);
    chromosphere::outer_boundary1_i.zeros(params->num_of_eq);
    chromosphere::inner_boundary0_i.zeros(params->num_of_eq);
    chromosphere::inner_boundary1_i.zeros(params->num_of_eq);
    
    // Cannot update xn here. Should do it elsewhere.
    array2vec(params->ds, chromosphere::ds_i, params->ns);
    array2vec(params->cellX, chromosphere::cellX_i, params->ns);
    array2vec(params->cellY, chromosphere::cellY_i, params->ns);
    array2vec(params->cellZ, chromosphere::cellZ_i, params->ns);
    array2vec(params->B_imh, chromosphere::B_imh, params->ns);
    array2vec(params->B_iph, chromosphere::B_iph, params->ns);
    array2vec(params->B_i, chromosphere::B_i, params->ns);
    array2vec(params->dinvB_ds, chromosphere::dinvB_ds_i, params->ns);
    array2vec(params->gPotential_imh, chromosphere::gPotential_imh, params->ns);
    array2vec(params->gPotential_iph, chromosphere::gPotential_iph, params->ns);
    array2vec(params->outer_bound0, chromosphere::outer_boundary0_i, params->num_of_eq);
    array2vec(params->outer_bound1, chromosphere::outer_boundary1_i, params->num_of_eq);
    array2vec(params->inner_bound0, chromosphere::inner_boundary0_i, params->num_of_eq);
    array2vec(params->inner_bound1, chromosphere::inner_boundary1_i, params->num_of_eq);

    // Initialize 2D arrays
    chromosphere::B_iimh.zeros(chromosphere::num_of_elem);
    chromosphere::B_iiph.zeros(chromosphere::num_of_elem);
    chromosphere::B_ii.zeros(chromosphere::num_of_elem);
    chromosphere::ds_ii.zeros(chromosphere::num_of_elem);
    chromosphere::dinvB_ds_ii.zeros(chromosphere::num_of_elem);
    for(int i=0; i<chromosphere::num_of_eq; i++) {
        chromosphere::B_iimh += chromosphere::scalar_to(chromosphere::B_imh, i);
        chromosphere::B_iiph += chromosphere::scalar_to(chromosphere::B_iph, i);
        chromosphere::B_ii += chromosphere::scalar_to(chromosphere::B_i, i);
        chromosphere::ds_ii += chromosphere::scalar_to(chromosphere::ds_i, i);
        chromosphere::dinvB_ds_ii += chromosphere::scalar_to(chromosphere::dinvB_ds_i, i);
    }
}

void chromo2params(struct Chromo_params* params) {
    // copy params to struct
    params->CNI = chromosphere::CNI;   // Ion density
    params->CNN = chromosphere::CNN;   // Neutral density
    params->CNV = chromosphere::CNV;   // Ion velocity
    params->CNU = chromosphere::CNU;   // Neutral velocity
    params->CEI = chromosphere::CEI;   // Ion energy
    params->CEN = chromosphere::CEN;   // Neutral energy

    // Primitive Variable Indices
    params->PNI = chromosphere::PNI;   // Ion density
    params->PNN = chromosphere::PNN;   // Neutral density
    params->PV  = chromosphere::PV;    // Ion velocity
    params->PU  = chromosphere::PU;    // Neutral velocity
    params->PPI = chromosphere::PPI;   // Ion pressure
    params->PPN = chromosphere::PPN;   // Neutral pressure

    params->ns = chromosphere::ns;
    params->num_of_eq = chromosphere::num_of_eq;
    params->num_of_elem = chromosphere::num_of_elem;

    params->gammamono = chromosphere::gammamono;
    params->alpha_p = chromosphere::alpha_p;
    params->pi = chromosphere::pi;

    params->m_i = chromosphere::m_i;
    params->m_n = chromosphere::m_n;
    params->m_e = chromosphere::m_e;
    params->g = chromosphere::g;
    params->mu_0 = chromosphere::mu_0;
    params->k_b = chromosphere::k_b;
    params->e_ = chromosphere::e_;

    params->CFL = chromosphere::CFL;

    // Allocate arrays
    Chromo_allocate_params_C(params);

    // Cannot update xn here.
    vec2array(chromosphere::ds_i, params->ds, chromosphere::ns);
    vec2array(chromosphere::cellX_i, params->cellX, chromosphere::ns);
    vec2array(chromosphere::cellY_i, params->cellY, chromosphere::ns);
    vec2array(chromosphere::cellZ_i, params->cellZ, chromosphere::ns);
    vec2array(chromosphere::B_imh, params->B_imh, chromosphere::ns);
    vec2array(chromosphere::B_iph, params->B_iph, chromosphere::ns);
    vec2array(chromosphere::B_i, params->B_i, chromosphere::ns);
    vec2array(chromosphere::dinvB_ds_i, params->dinvB_ds, chromosphere::ns);
    vec2array(chromosphere::gPotential_iph, params->gPotential_iph, chromosphere::ns);
    vec2array(chromosphere::gPotential_imh, params->gPotential_imh, chromosphere::ns);
    vec2array(chromosphere::outer_boundary0_i, params->outer_bound0, chromosphere::num_of_eq);
    vec2array(chromosphere::outer_boundary1_i, params->outer_bound1, chromosphere::num_of_eq);
    vec2array(chromosphere::inner_boundary0_i, params->inner_bound0, chromosphere::num_of_eq);
    vec2array(chromosphere::inner_boundary1_i, params->inner_bound1, chromosphere::num_of_eq);
}

extern "C" int Chromo_sub2ind(const int* ns, const int* i, const int* k) {
    return (int)arma::sub2ind(
            arma::size(*ns, chromosphere::num_of_eq), *i, *k);
}

extern "C" void Chromo_init_CPP(struct Chromo_params* params) {
    cout << "----------Initializing chromosphere.----------" << endl;
    params2chromo(params);
    Vec xn_vec(chromosphere::num_of_elem);
    // copy xn_vec to xn
    // vec2array(xn_vec, xn, chromosphere::num_of_elem);
    array2vec(params->xn, xn_vec, chromosphere::num_of_elem);
    chromo2params(params);

    cout << "------end chromo_init_cpp----------" << endl;
}

extern "C" void Chromo_cal_dt(const struct Chromo_params* params) {
    // xn intent: in
    // dt intent: out
    Vec xn_vec(chromosphere::num_of_elem), dt_vec(chromosphere::ns);
    array2vec(params->xn, xn_vec, chromosphere::num_of_elem);
    chromosphere::CFL = params->CFL;

    dt_vec = chromosphere::cal_dt_i(xn_vec);
    vec2array(dt_vec, params->dt, chromosphere::ns);
}

// By Keheng: use Euler method to advance the state.
// xn is an array of size nx*nz*num_of_eq, dt is the time step.
extern "C" void Chromo_advance_Euler(struct Chromo_params* params) {
    // xn intent: inout
    // dt intent: in
    Vec xn_vec(chromosphere::num_of_elem);
    Vec dt_vec(chromosphere::ns);

    params2chromo(params);
    array2vec(params->xn, xn_vec, chromosphere::num_of_elem);
    array2vec(params->dt, dt_vec, chromosphere::ns);

    // advance the state
    xn_vec = chromosphere::advance_Euler_ii(xn_vec, dt_vec);

    vec2array(xn_vec, params->xn, chromosphere::num_of_elem);
    chromo2params(params);
}

extern "C" void Chromo_prim2cons(float* p_xn) {
    // p_xn intent: inout
    Vec p_xn_vec(chromosphere::num_of_elem), c_xn_vec(chromosphere::num_of_elem);
    array2vec(p_xn, p_xn_vec, chromosphere::num_of_elem);
    c_xn_vec = chromosphere::prim2cons(p_xn_vec);
    vec2array(c_xn_vec, p_xn, chromosphere::num_of_elem);
}

extern "C" void Chromo_cons2prim(float* c_xn) {
    // c_xn intent: inout
    Vec c_xn_vec(chromosphere::num_of_elem), p_xn_vec(chromosphere::num_of_elem);
    array2vec(c_xn, c_xn_vec, chromosphere::num_of_elem);
    p_xn_vec = chromosphere::cons2prim(c_xn_vec);
    vec2array(p_xn_vec, c_xn, chromosphere::num_of_elem);
}

extern "C" void Chromo_allocate_params_C(struct Chromo_params* pC) {
    if(!pC->is_allocated) {
        // cout << "allocated memory in C++" << endl;
        pC->ds = new float[pC->ns];
        pC->dt = new float[pC->ns];
        pC->xn = new float[pC->num_of_elem];
        pC->cellX = new float[pC->ns];
        pC->cellY = new float[pC->ns];
        pC->cellZ = new float[pC->ns];
        pC->B_imh = new float[pC->ns];
        pC->B_iph = new float[pC->ns];
        pC->B_i = new float[pC->ns];
        pC->dinvB_ds = new float[pC->ns];
        pC->gPotential_imh = new float[pC->ns];
        pC->gPotential_iph = new float[pC->ns];
        pC->outer_bound0 = new float[pC->num_of_eq];
        pC->outer_bound1 = new float[pC->num_of_eq];
        pC->inner_bound0 = new float[pC->num_of_eq];
        pC->inner_bound1 = new float[pC->num_of_eq];
        pC->is_allocated = 1;
    }
}

extern "C" void Chromo_deallocate_params_C(struct Chromo_params* pC) {
    if(pC->is_allocated) {
        // cout << "deallocated memory in C++" << endl;
        delete[] pC->ds;
        delete[] pC->dt;
        delete[] pC->xn;
        delete[] pC->cellX;
        delete[] pC->cellY;
        delete[] pC->cellZ;
        delete[] pC->B_imh;
        delete[] pC->B_iph;
        delete[] pC->B_i;
        delete[] pC->dinvB_ds;
        delete[] pC->gPotential_imh;
        delete[] pC->gPotential_iph;
        delete[] pC->outer_bound0;
        delete[] pC->outer_bound1;
        delete[] pC->inner_bound0;
        delete[] pC->inner_bound1;
        pC->is_allocated = 0;
    }
}
