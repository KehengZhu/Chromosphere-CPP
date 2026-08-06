#include "chromosphere.hpp"
#include "physics.hpp"
#include "profiling.hpp"
#include "parallel.hpp"
#include "scenarios/scenario.hpp"

#include <armadillo>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

bool env_enabled(const char* name, bool default_value = false) {
    const char* value = std::getenv(name);
    if (!value) return default_value;
    return std::string(value) != "0" && std::string(value) != "false"
        && std::string(value) != "off";
}

uint32_t rotr(uint32_t x, unsigned n) { return (x >> n) | (x << (32-n)); }

std::string sha256_file(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot hash Gamma1 table: "+path);
    std::vector<uint8_t> bytes((std::istreambuf_iterator<char>(input)), {});
    const uint64_t bit_length = static_cast<uint64_t>(bytes.size())*8;
    bytes.push_back(0x80);
    while (bytes.size()%64 != 56) bytes.push_back(0);
    for (int shift = 56; shift >= 0; shift -= 8)
        bytes.push_back(static_cast<uint8_t>(bit_length >> shift));
    static constexpr uint32_t k[64] = {
        0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
        0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
        0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
        0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
        0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
        0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
        0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
        0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
    std::array<uint32_t,8> h = {0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
                                0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    for (std::size_t off = 0; off < bytes.size(); off += 64) {
        uint32_t w[64];
        for (int i = 0; i < 16; ++i)
            w[i] = (uint32_t(bytes[off+4*i])<<24)|(uint32_t(bytes[off+4*i+1])<<16)
                 | (uint32_t(bytes[off+4*i+2])<<8)|uint32_t(bytes[off+4*i+3]);
        for (int i = 16; i < 64; ++i) {
            const uint32_t s0 = rotr(w[i-15],7)^rotr(w[i-15],18)^(w[i-15]>>3);
            const uint32_t s1 = rotr(w[i-2],17)^rotr(w[i-2],19)^(w[i-2]>>10);
            w[i] = w[i-16]+s0+w[i-7]+s1;
        }
        uint32_t a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],hh=h[7];
        for (int i = 0; i < 64; ++i) {
            const uint32_t s1=rotr(e,6)^rotr(e,11)^rotr(e,25);
            const uint32_t ch=(e&f)^((~e)&g);
            const uint32_t t1=hh+s1+ch+k[i]+w[i];
            const uint32_t s0=rotr(a,2)^rotr(a,13)^rotr(a,22);
            const uint32_t maj=(a&b)^(a&c)^(b&c);
            const uint32_t t2=s0+maj;
            hh=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
        }
        h[0]+=a;h[1]+=b;h[2]+=c;h[3]+=d;h[4]+=e;h[5]+=f;h[6]+=g;h[7]+=hh;
    }
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (uint32_t word : h) out << std::setw(8) << word;
    return out.str();
}

double gamma_solver_conductivity(const chromosphere::Grid& grid,
                                 const chromosphere::MixtureThermo& th,
                                 double local_spacing) {
    return chromosphere::solver_effective_conductivity(
        grid, th.n_e, th.n_HI, th.T,
        chromosphere::equilibrium_heat_capacity(th.rho, th.T), local_spacing);
}

} // namespace

int main(int argc, char** argv) {
    using namespace chromosphere;

    // Usage: chromo_main [output_path] [mode] [ionization] [scenario] [data_path] [time_mult] [cooling]
    //   output_path:  defaults to "output.txt"
    //   mode:         "full" (semi-implicit, default) or "explicit" (R_I ≡ 0)
    //   ionization:   "ionization" (default — Stage E enabled) or "no-ionization"
    //   scenario:     "model_column" (default — the unified chromosphere→corona column;
    //                 accepts "model_isentropic" / "model_gentle" as aliases) |
    //                 "model_c7" | "model_flare" | "pfss_field_line" | "analytic_canopy"
    //   data_path:    required for tabulated scenarios (e.g. pfss_field_line);
    //                 ignored otherwise (pass "-" or "" for scenarios that don't use it)
    //   time_mult:    multiplier on the default total_time (default 1.0). Step
    //                 cap scales accordingly so a longer run is not truncated.
    //   cooling:      "cooling" (Stage R: CL2012 optically-thick + optically-thin
    //                 radiative cooling, enabled by DEFAULT) or "no-cooling"
    const std::string out_path      = (argc > 1) ? argv[1] : "outputs/output.txt";
    const std::string mode          = (argc > 2) ? argv[2] : "full";
    const std::string ioniz_arg     = (argc > 3) ? argv[3] : "ionization";
    const std::string scenario_name = (argc > 4) ? argv[4] : "model_column";
    const std::string data_path     = (argc > 5) ? argv[5] : "";
    const float       time_mult     = (argc > 6) ? std::stof(argv[6]) : 1.0f;
    const std::string cool_arg      = (argc > 7) ? argv[7] : "cooling";
    const bool explicit_only        = (mode == "explicit");
    const bool ionization_on        = (ioniz_arg != "no-ionization");
    const bool cooling_on           = (cool_arg != "no-cooling");
    const bool profiling_on         = env_enabled("CHROMO_PROFILE");
    const bool eos_counting_on       = env_enabled("CHROMO_EOS_COUNTS");
    const bool write_output         = env_enabled("CHROMO_OUTPUT", true);
    const bool write_gamma_diag     = env_enabled("CHROMO_GAMMA_DIAG", true);
    // Read-only face-flux diagnostic (default OFF). CHROMO_FACE_FLUX_DIAG=1 writes
    // a <out>.faceflux sidecar holding the PRODUCTION Rusanov total-mass face flux
    // and its central/diffusive split, taken straight out of rhs_explicit_mixture.
    // CHROMO_FACE_FLUX_TOP limits each record to the top N cells (0 = all cells);
    // CHROMO_FACE_FLUX_STRIDE sets the capture cadence in steps (default: the
    // snapshot stride). Enabling it changes no numerical result.
    const bool write_face_flux      = env_enabled("CHROMO_FACE_FLUX_DIAG", false);
    // Read-only outer-face conduction diagnostic (default OFF). =1 writes a
    // <out>.outercond sidecar with the OUTER-face quantities of the final converged
    // conduction Newton iteration: T_top, T_wall, the physical and numerical face
    // conductivities, and the corresponding q_phys / q_num / q_total. One row per
    // captured step. Enabling it changes no numerical result.
    const bool write_outer_cond     = env_enabled("CHROMO_OUTER_COND_DIAG", false);
    set_runtime_profiling(profiling_on);
    reset_runtime_profile();
    set_eos_operation_counting(eos_counting_on);
    reset_eos_operation_counts();
    require_supported_parallel_runtime();
    std::cerr << "[openmp] compiled=" << (openmp_compiled() ? 1 : 0)
              << " max_threads=" << parallel_max_threads()
              << " conduction_team="
              << std::min(kMaximumConductionThreads, parallel_max_threads());
    if (const char* requested = std::getenv("OMP_NUM_THREADS"))
        std::cerr << " requested=" << requested;
    std::cerr << '\n';

    float cfl = 0.25f;
    if (const char* e = std::getenv("CHROMO_CFL")) { try { cfl = std::stof(e); } catch (...) {} }

    Scenario sc;
    try {
        sc = make_scenario(scenario_name, data_path);
    } catch (const std::exception& e) {
        std::cerr << "scenario error: " << e.what() << std::endl;
        return 1;
    }

    Grid grid;
    grid.init(sc.peek_ns(), cfl);
    std::string gamma_table_path;
    std::string gamma_table_sha256;
    if (const char* gamma_path = std::getenv("GAMMA_TABLE")) {
        try {
            if (std::getenv("ISO_GAMMA"))
                throw std::logic_error("GAMMA_TABLE and explicitly supplied ISO_GAMMA are mutually exclusive");
            gamma_table_path = gamma_path;
            gamma_table_sha256 = sha256_file(gamma_table_path);
            grid.eos_gamma_table = EosGammaTable::load(gamma_table_path);
            if (explicit_only)
                throw std::logic_error("GAMMA_TABLE requires the full Euler integrator, not explicit mode");
            if (scenario_name != "model_column" && scenario_name != "model_isentropic")
                throw std::logic_error("GAMMA_TABLE currently supports only model_column/model_isentropic");
        } catch (const std::exception& e) {
            std::cerr << "gamma-table mode error: " << e.what() << std::endl;
            return 1;
        }
    }
    grid.enable_ionization = ionization_on;
    grid.enable_radiative_cooling = cooling_on;
    // "neutrals off" experiment: SINGLE_FLUID=1 slaves neutrals to the ion fluid
    // (single-fluid limit). Default (unset/0) is the full two-fluid model.
    if (const char* e = std::getenv("SINGLE_FLUID")) {
        try { grid.single_fluid = (std::stof(e) != 0.0f); } catch (...) {}
    }
    // Separate electron temperature T_e ≠ T_i (docs/electron_temperature_plan.md).
    // ENABLE_TE=1 switches on the three-temperature model; default (unset/0) is
    // the single-temperature baseline that reproduces the pre-T_e physics.
    if (const char* e = std::getenv("ENABLE_TE")) {
        try { grid.enable_Te = (std::stof(e) != 0.0f); } catch (...) {}
    }
    if (!grid.eos_gamma_table.empty()) {
        if (ionization_on) {
            std::cerr << "gamma-table mode error: finite-rate ionization was requested; "
                         "pass no-ionization" << std::endl;
            return 1;
        }
        if (grid.enable_Te) {
            std::cerr << "gamma-table mode error: ENABLE_TE=1 is incompatible with "
                         "the common-temperature equilibrium closure" << std::endl;
            return 1;
        }
        if (std::getenv("SINGLE_FLUID") && !grid.single_fluid) {
            std::cerr << "gamma-table mode error: explicit SINGLE_FLUID=0 is incompatible "
                         "with equilibrium projection" << std::endl;
            return 1;
        }
    }
    Vec xn = sc.ic(grid);
    if (!grid.eos_gamma_table.empty()
        && (!grid.single_fluid || grid.enable_Te || grid.enable_ionization)) {
        std::cerr << "gamma-table mode error: scenario toggles must leave single_fluid=1, "
                     "ENABLE_TE=0, and ionization disabled" << std::endl;
        return 1;
    }

    const float c_s_target = 2.0e4f; // m/s, ion sound speed scale (writeup §4)
    float total_time = time_mult * 10.0f * arma::sum(grid.ds_i) / c_s_target;
    // Optional absolute stop time for reproducible long diagnostics. Unset keeps
    // the historical time_mult-derived behavior byte-for-byte unchanged.
    if (const char* end = std::getenv("CHROMO_T_END")) {
        try {
            const float requested = std::stof(end);
            if (requested > 0.0f && std::isfinite(requested)) total_time = requested;
        } catch (...) {}
    }
    const int   step_cap   = static_cast<int>(std::max(10000.0f, 10000.0f * time_mult));

    std::ofstream logf("outputs/output.log", std::ios::app);
    logf << "[" << out_path << " mode=" << mode << "] Total time = " << total_time << std::endl;
    std::cout << "[" << mode << "] Total time = " << total_time << std::endl;

    // Multi-snapshot output. Format:
    //   line 1: "ns num_of_eq"
    //   line 2: ns cumulative cell heights in km (offset by grid.out_base_km, the
    //           domain base — 1003 km for C7-based scenarios, 0 km for the
    //           photosphere-anchored model_isentropic C7 IC; ds_i is in m, hence
    //           the 1e-3 conversion)
    //   then, repeated: a "# t = T step = S" marker followed by ns lines of
    //   num_of_eq space-separated conserved-variable values.
    std::ofstream fout;
    if (write_output) fout.open(out_path);
    if (fout) {
        fout << grid.ns << " " << num_of_eq << '\n';
        float cum_km = 0.0f;
        for (arma::uword i = 0; i < grid.ns; ++i) {
            cum_km += grid.ds_i(i) * 1.0e-3f;
            fout << "  " << (cum_km + grid.out_base_km);
        }
        fout << '\n';
    }
    if (fout && !grid.eos_gamma_table.empty() && write_gamma_diag)
        fout << "# EOS_MODE=gamma_table diagnostics=" << out_path << ".gamma_diag\n";
    std::ofstream gamma_diag;
    if (!grid.eos_gamma_table.empty() && write_gamma_diag) {
        gamma_diag.open(out_path+".gamma_diag");
        if (!gamma_diag) throw std::runtime_error("cannot open gamma diagnostic sidecar");
        gamma_diag << grid.ns << " 10\n";
        float cum_km = 0.0f;
        for (arma::uword i = 0; i < grid.ns; ++i) {
            cum_km += grid.ds_i(i)*1.0e-3f;
            gamma_diag << "  " << (cum_km+grid.out_base_km);
        }
        gamma_diag << "\n# EOS_MODE=gamma_table\n"
                   << "# GAMMA_TABLE=" << gamma_table_path << "\n"
                   << "# GAMMA_TABLE_SHA256=" << gamma_table_sha256 << "\n"
                   << "# potential=cell_center_mean_of_phi_g_imh_phi_g_iph\n"
                   << "# columns=rho_total v_cm T x_eq n_e n_HI p_total Gamma1 kappa_physical kappa_solver\n";
    }

    auto write_frame = [&](float t_now, int step_now) {
        ProfileScope output_timer(ProfileRegion::Output);
        if (fout) {
            fout << "# t = " << t_now << " step = " << step_now << '\n';
            for (arma::uword i = 0; i < grid.ns; ++i) {
                for (arma::uword k = 0; k < num_of_eq; ++k) {
                    fout << "  " << xn(arma::sub2ind(
                        arma::size(grid.ns,num_of_eq),i,k));
                }
                fout << '\n';
            }
        }
        if (gamma_diag) {
            gamma_diag << "# t = " << t_now << " step = " << step_now << '\n';
            const auto sz = arma::size(grid.ns, num_of_eq);
            for (arma::uword i = 0; i < grid.ns; ++i) {
                auto at = [&](arma::uword k) -> double {
                    return static_cast<double>(xn(arma::sub2ind(sz, i, k)));
                };
                const double phi = 0.5*static_cast<double>(grid.phi_g_imh(i)+grid.phi_g_iph(i));
                const MixtureThermo th = decode_equilibrium_mixture(
                    grid.eos_gamma_table, at(cons::RHO_I), at(cons::RHO_N),
                    at(cons::MOM_I), at(cons::MOM_N), at(cons::E_I), at(cons::E_N), phi,
                    std::numeric_limits<double>::quiet_NaN(), grid.eos_gamma_debug_clamp);
                const double velocity = (at(cons::MOM_I)+at(cons::MOM_N))/th.rho;
                gamma_diag << "  " << th.rho << "  " << velocity << "  " << th.T
                           << "  " << th.x_eq << "  " << th.n_e << "  " << th.n_HI
                           << "  " << (th.p_i+th.p_n) << "  " << th.gamma1
                           << "  " << physical_conductivity(th.n_e,th.n_HI,th.T)
                           << "  " << gamma_solver_conductivity(
                                  grid, th, grid.ds_i(i)) << '\n';
            }
        }
    };

    // Scale snapshot cadence with run length so that frame count stays
    // bounded (~500 frames per run) regardless of time_mult. Keeps render
    // time roughly constant when sweeping time_mult — animation length at
    // 30 fps is ~17 s independent of physical simulation duration.
    int frame_stride = std::max(10,
        static_cast<int>(10.0f * std::max(1.0f, time_mult / 5.0f)));
    if (const char* stride = std::getenv("CHROMO_FRAME_STRIDE")) {
        try { frame_stride = std::max(1, std::stoi(stride)); } catch (...) {}
    }
    // Face-flux sidecar. One record per captured step, written AFTER the step but
    // holding the RHS evaluation of the state at the START of that step (the flux
    // that produced the step). Every column comes from the one production RHS call,
    // so a record is self-consistent and needs no cross-matching with the snapshot
    // file. Face index i is the upper face i+1/2 of cell i; the last row is the
    // outer boundary face.
    std::ofstream face_flux;
    int  face_flux_stride = frame_stride;
    long face_flux_top    = 0;
    arma::uword face_flux_lo = 0;
    Vec face_cell_km, face_face_km;
    if (const char* e = std::getenv("CHROMO_FACE_FLUX_STRIDE")) {
        try { face_flux_stride = std::max(1, std::stoi(e)); } catch (...) {}
    }
    if (const char* e = std::getenv("CHROMO_FACE_FLUX_TOP")) {
        try { face_flux_top = std::max(0L, std::stol(e)); } catch (...) {}
    }
    if (write_face_flux) {
        if (grid.eos_gamma_table.empty())
            throw std::runtime_error(
                "CHROMO_FACE_FLUX_DIAG requires gamma-table mode (the capture lives "
                "in rhs_explicit_mixture)");
        face_flux.open(out_path+".faceflux");
        if (!face_flux) throw std::runtime_error("cannot open face-flux sidecar");
        float cum_km = 0.0f;
        face_cell_km.set_size(grid.ns);
        face_face_km.set_size(grid.ns);
        for (arma::uword i = 0; i < grid.ns; ++i) {
            face_cell_km(i) = cum_km + 0.5f*grid.ds_i(i)*1.0e-3f + grid.out_base_km;
            cum_km += grid.ds_i(i)*1.0e-3f;
            face_face_km(i) = cum_km + grid.out_base_km;
        }
        face_flux_lo = (face_flux_top > 0
            && static_cast<arma::uword>(face_flux_top) < grid.ns)
            ? grid.ns - static_cast<arma::uword>(face_flux_top) : 0;
        face_flux << "# face-flux diagnostic (read-only capture of rhs_explicit_mixture)\n"
                  << "# ns=" << grid.ns << " first_cell=" << face_flux_lo
                  << " uniform_mesh=" << grid.uniform_mesh
                  << " mc3=" << grid.mc3_limiter << " beta=" << grid.limiter_beta
                  << " eq_wb=" << grid.eq_wb << " ds_km=" << (grid.ds_i(0)*1.0e-3f)
                  << "\n# mass flux = RHO_I row + RHO_N row (the conserved total)\n"
                  << "# columns=cell_km face_km rho_cell v_cell T_cell "
                     "rho_L rho_R v_L v_R T_L T_R cs_L cs_R a_face "
                     "f_central f_diff f_total eq_residual_mass "
                     "r_rho phi_plus_rho r_ip1_rho phi_minus_rho "
                     "r_v phi_plus_v r_T phi_plus_T\n";
        face_flux.precision(10);
    }

    // Outer-face conduction sidecar. One row per captured step, holding the outer-face
    // quantities of the FINAL CONVERGED Newton iteration of that step's conduction
    // solve. t/step label the START of the step, matching the .faceflux convention.
    std::ofstream outer_cond;
    int outer_cond_stride = face_flux_stride;
    if (const char* e = std::getenv("CHROMO_OUTER_COND_STRIDE")) {
        try { outer_cond_stride = std::max(1, std::stoi(e)); } catch (...) {}
    }
    if (write_outer_cond) {
        if (grid.eos_gamma_table.empty())
            throw std::runtime_error(
                "CHROMO_OUTER_COND_DIAG requires gamma-table mode (the capture lives "
                "in apply_gamma_conduction_stage)");
        outer_cond.open(out_path+".outercond");
        if (!outer_cond) throw std::runtime_error("cannot open outer-conduction sidecar");
        outer_cond << "# outer-face conduction diagnostic (read-only capture of the "
                      "final converged apply_gamma_conduction_stage iteration)\n"
                   << "# ns=" << grid.ns << " ds_km=" << (grid.ds_i(grid.ns-1)*1.0e-3f)
                   << " numerical_diffusivity_per_length="
                   << grid.numerical_diffusivity_per_length
                   << " uniform_mesh=" << grid.uniform_mesh
                   << " impose_outer_heat_flux=" << grid.impose_outer_heat_flux << '\n'
                   << "# q = (kappa_phys_face + kappa_num_face)*(T_wall - T_top)/ds_face"
                      "  [W m^-2, positive = into the top cell]\n"
                   << "# columns=t step T_top T_wall kappa_phys_face chi_num_face "
                      "kappa_num_face "
                      "ds_face area_ratio q_phys q_num q_total\n";
        outer_cond.precision(10);
    }

    auto write_outer_cond_record = [&](float t_now, int step_now) {
        const OuterConductionCapture& oc = grid.outer_conduction_capture;
        if (!outer_cond || !oc.valid) return;
        outer_cond << t_now << ' ' << step_now << ' '
                   << oc.T_top << ' ' << oc.T_wall << ' '
                   << oc.kappa_phys_face << ' ' << oc.chi_num_face << ' '
                   << oc.kappa_num_face << ' '
                   << oc.ds_face << ' ' << oc.area_ratio << ' '
                   << oc.q_phys << ' ' << oc.q_num << ' ' << oc.q_total << '\n';
    };

    auto write_face_record = [&](float t_now, int step_now) {
        const GammaFaceFluxCapture& c = grid.face_flux_capture;
        if (!face_flux || !c.valid) return;
        {
            face_flux << "# t = " << t_now << " step = " << step_now << '\n';
            for (arma::uword i = face_flux_lo; i < grid.ns; ++i) {
                face_flux << face_cell_km(i) << ' ' << face_face_km(i) << ' '
                    << c.rho_cell[i] << ' ' << c.v_cell[i] << ' ' << c.T_cell[i] << ' '
                    << c.rho_L[i] << ' ' << c.rho_R[i] << ' '
                    << c.v_L[i] << ' ' << c.v_R[i] << ' '
                    << c.T_L[i] << ' ' << c.T_R[i] << ' '
                    << c.cs_L[i] << ' ' << c.cs_R[i] << ' ' << c.a_face[i] << ' '
                    << c.f_central[i] << ' ' << c.f_diff[i] << ' ' << c.f_total[i] << ' '
                    << c.eq_residual_mass[i] << ' '
                    << c.r_rho[i] << ' ' << c.phi_plus_rho[i] << ' '
                    << c.r_ip1_rho[i] << ' ' << c.phi_minus_rho[i] << ' '
                    << c.r_v[i] << ' ' << c.phi_plus_v[i] << ' '
                    << c.r_T[i] << ' ' << c.phi_plus_T[i] << '\n';
            }
        }
    };

    float     time         = 0.0f;
    int       step         = 0;
    std::uint64_t state_generation = 0;
    DecodedMixtureField decoded_storage[2];
    DecodedMixtureField* decoded=&decoded_storage[0];
    const DecodedMixtureField* previous_decoded=nullptr;
    write_frame(time, step);

    while (time < total_time && step < step_cap) {
        grid.sim_time = time;   // expose current time to time-dependent terms (beam window)
        {
            ProfileScope boundary_timer(ProfileRegion::Boundary);
            sc.update_bc(grid, xn);
        }
        Vec dt;
        if (!grid.eos_gamma_table.empty()) {
            decode_mixture_field_into(
                grid,xn,*decoded,state_generation,previous_decoded);
            dt = cal_dt_i(grid,xn,*decoded);
        } else {
            dt = cal_dt_i(grid,xn);
        }
        const float dt_avg = arma::mean(dt);
        // Arm the read-only face-flux capture for THIS step's single RHS call. dt is
        // already known here, so the final step of the run is always captured.
        if (write_face_flux)
            grid.capture_face_flux = (step % face_flux_stride == 0)
                || (time + dt_avg >= total_time) || (step + 1 >= step_cap);
        const bool capture_cond_now = write_outer_cond
            && ((step % outer_cond_stride == 0)
                || (time + dt_avg >= total_time) || (step + 1 >= step_cap));
        grid.capture_outer_conduction = capture_cond_now;

        if (step % 100 == 0) {
            std::cout << "[" << mode << "] step = " << step
                      << "  dt = " << dt_avg
                      << "  time = " << time
                      << "  T_c = " << grid.trac_cutoff_T << std::endl;
        }

        if (!grid.eos_gamma_table.empty())
            xn = advance_Euler_state(grid,xn,dt,*decoded);
        else
            xn = explicit_only ? advance_Euler_explicit_state(grid,xn,dt)
                               : advance_Euler_state(grid,xn,dt);
        if (grid.capture_face_flux) {
            write_face_record(time, step);   // t/step of the state the RHS was taken at
            grid.capture_face_flux = false;
        }
        if (capture_cond_now) {
            write_outer_cond_record(time, step);
            grid.capture_outer_conduction = false;
        }
        time += dt_avg;
        ++step;
        ++state_generation;
        if (!grid.eos_gamma_table.empty()) {
            previous_decoded=decoded;
            decoded = decoded == &decoded_storage[0]
                ? &decoded_storage[1] : &decoded_storage[0];
        }

        if (step % frame_stride == 0) write_frame(time, step);
    }
    if (step % frame_stride != 0) write_frame(time, step);

    std::cout << "[" << mode << "] END! step=" << step << " time=" << time << std::endl;
    print_runtime_profile(std::cout);
    if (eos_counting_on) {
        const EosOperationCounts counts = eos_operation_counts();
        std::cout << "eos.gamma1_queries=" << counts.gamma1_queries << '\n'
                  << "eos.temperature_logs=" << counts.temperature_logs << '\n'
                  << "eos.n_h_logs=" << counts.n_h_logs << std::endl;
    }
    return 0;
}
