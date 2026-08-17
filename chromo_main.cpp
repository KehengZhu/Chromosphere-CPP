/*!
 * @file chromo_main.cpp
 * @brief Command-line driver: parse arguments and environment, build the
 *        scenario, run the timestep loop, and write snapshots and sidecars.
 * @ingroup driver
 *
 * The driver owns no physics. It selects a solver, steps it, and writes files:
 *
 *  1. make_scenario() builds the requested IC/BC pair.
 *  2. Grid::init() allocates the mesh; the scenario IC fills the geometry,
 *     gravity, ghost buffers and the initial conserved state (and may resize
 *     the Grid when it generates a locally refined mesh).
 *  3. A loaded `Gamma1` table selects the RELEASE single-fluid solver
 *     (mixture_advance()); without one the driver runs the historical
 *     two-fluid solver (advance_Euler_state() / advance_Euler_explicit_state()).
 *     The release path rejects `SINGLE_FLUID`, `ENABLE_TE`, `ISO_TWO_FLUID`,
 *     `ISO_IONIZATION` and finite-rate ionization outright.
 *  4. Each step refreshes the boundary, computes the CFL timestep, advances,
 *     and writes a snapshot when the OutputSchedule says it is due.
 *  5. On exit it prints `termination=end_time|step_cap|other` together with the
 *     requested end time, final time, final step and effective cap.
 *
 * The full command-line and environment reference is on the
 * @ref configuration "Configuration reference" page; the file formats written
 * here are described on the @ref io_formats "Output formats" page.
 */
#include "chromosphere.hpp"
#include "single_fluid/mixture.hpp"
#include "two_fluid/two_fluid.hpp"
#include "physics.hpp"
#include "profiling.hpp"
#include "parallel.hpp"
#include "run_control.hpp"
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

} // namespace

/// Command-line entry point.
///
/// Reads the positional arguments and the `CHROMO_*` / `ISO_*` environment,
/// builds the requested Scenario, allocates and initialises the Grid, then runs
/// the timestep loop — boundary refresh, CFL timestep, solver advance, snapshot
/// when due — until the end time or the step cap is reached, writing the
/// snapshot, log and sidecar files on the way and the `termination=` summary at
/// the end. See the @ref configuration "Configuration reference" for the full
/// argument and environment list.
///
/// @param argc Argument count.
/// @param argv `[output_path] [mode] [ionization] [scenario] [data_path]
///             [time_mult] [cooling]`, each optional and positional; the
///             defaults are documented at the top of the function body.
/// @return 0 on a completed run; non-zero after a fatal configuration or
///         solver error, whose message is printed to stderr.
int main(int argc, char** argv) {
    using namespace chromosphere;

    // Usage: chromo_main [output_path] [mode] [ionization] [scenario] [data_path] [time_mult] [cooling]
    //   output_path:  defaults to "output.txt"
    //   mode:         "full" (semi-implicit, default) or "explicit" (R_I ≡ 0)
    //   ionization:   "ionization" or "no-ionization" (default: no-ionization for
    //                 the default model_column release, ionization otherwise)
    //   scenario:     "model_column" (default — the single-fluid RELEASE column) |
    //                 "model_gentle" | "model_c7" | "model_flare" |
    //                 "pfss_field_line" | "analytic_canopy"  (all two-fluid)
    //   data_path:    required for tabulated scenarios (e.g. pfss_field_line);
    //                 ignored otherwise (pass "-" or "" for scenarios that don't use it)
    //   time_mult:    multiplier on the default total_time (default 1.0). Step
    //                 cap scales accordingly so a longer run is not truncated.
    //   cooling:      "cooling" or "no-cooling" (default: no-cooling for the
    //                 default model_column release, cooling otherwise)
    const std::string out_path      = (argc > 1) ? argv[1] : "outputs/output.txt";
    const std::string mode          = (argc > 2) ? argv[2] : "full";
    const std::string scenario_name = (argc > 4) ? argv[4] : "model_column";
    const bool release_column = scenario_name == "model_column";
    const std::string ioniz_arg     = (argc > 3) ? argv[3]
                                                 : (release_column ? "no-ionization" : "ionization");
    const std::string data_path     = (argc > 5) ? argv[5] : "";
    const float       time_mult     = (argc > 6) ? std::stof(argv[6]) : 1.0f;
    const std::string cool_arg      = (argc > 7) ? argv[7]
                                                 : (release_column ? "no-cooling" : "cooling");
    const bool explicit_only        = (mode == "explicit");
    const bool ionization_on        = (ioniz_arg != "no-ionization");
    const bool cooling_on           = (cool_arg != "no-cooling");
    const bool profiling_on         = env_enabled("CHROMO_PROFILE");
    const bool eos_counting_on       = env_enabled("CHROMO_EOS_COUNTS");
    const bool write_output         = env_enabled("CHROMO_OUTPUT", true);
    const bool write_gamma_diag     = env_enabled("CHROMO_GAMMA_DIAG", true);
    // Read-only face-flux diagnostic (default OFF). CHROMO_FACE_FLUX_DIAG=1 writes
    // a <out>.faceflux sidecar holding the PRODUCTION total-mass face flux
    // and its central/diffusive split, taken straight out of mixture_rhs_explicit.
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
                throw std::logic_error("the release solver has no explicit-only mode");
            if (scenario_name != "model_column")
                throw std::logic_error("the release solver currently supports only model_column");
        } catch (const std::exception& e) {
            std::cerr << "release-solver error: " << e.what() << std::endl;
            return 1;
        }
    }
    // A loaded Gamma1 table selects the RELEASE single-fluid solver; without it
    // the driver runs the historical two-fluid research solver.
    const bool release_mode = !grid.eos_gamma_table.empty();
    grid.enable_ionization = ionization_on;
    grid.enable_radiative_cooling = cooling_on;
    // Legacy two-fluid research knobs. They have no meaning for the release
    // solver, which is a single-fluid common-temperature equilibrium mixture.
    if (!release_mode) {
        // "neutrals off" experiment: SINGLE_FLUID=1 slaves neutrals to the ion
        // fluid. Default (unset/0) is the full two-fluid model.
        if (const char* e = std::getenv("SINGLE_FLUID")) {
            try { grid.single_fluid = (std::stof(e) != 0.0f); } catch (...) {}
        }
        // Separate electron temperature T_e != T_i.
        if (const char* e = std::getenv("ENABLE_TE")) {
            try { grid.enable_Te = (std::stof(e) != 0.0f); } catch (...) {}
        }
    } else {
        if (ionization_on) {
            std::cerr << "release-solver error: finite-rate ionization was requested, "
                         "but the release closure is Saha equilibrium; pass "
                         "no-ionization" << std::endl;
            return 1;
        }
        for (const char* legacy : {"SINGLE_FLUID", "ENABLE_TE", "ISO_TWO_FLUID",
                                   "ISO_IONIZATION"}) {
            if (std::getenv(legacy)) {
                std::cerr << "release-solver error: " << legacy
                          << " is a legacy two-fluid setting and has no meaning for "
                             "the single-fluid release solver" << std::endl;
                return 1;
            }
        }
    }
    Vec xn = sc.ic(grid);
    if (release_mode && xn.n_elem != grid.n_mixture_state) {
        std::cerr << "release-solver error: the scenario returned a state of "
                  << xn.n_elem << " elements, expected " << grid.n_mixture_state
                  << " (rho, rho u, E per cell)" << std::endl;
        return 1;
    }

    const float c_s_target = 2.0e4f; // m/s, ion sound speed scale (writeup §4)
    float total_time = time_mult * 10.0f * arma::sum(grid.ds_i) / c_s_target;
    // Optional absolute stop time for reproducible long diagnostics. Unset keeps
    // the historical time_mult-derived behavior byte-for-byte unchanged.
    const char* end_time_env = std::getenv("CHROMO_T_END");
    bool end_time_requested = false;
    if (end_time_env) {
        try {
            const float requested = std::stof(end_time_env);
            if (requested > 0.0f && std::isfinite(requested)) {
                total_time = requested;
                end_time_requested = true;
            }
        } catch (...) {}
    }
    // Step cap. An explicit CHROMO_STEP_CAP always wins. Without it, a run that
    // sets CHROMO_T_END is governed by the end time (the legacy time_mult cap
    // would truncate it); a run without CHROMO_T_END keeps the legacy cap.
    const StepCapSelection cap_sel = select_step_cap(
        std::getenv("CHROMO_STEP_CAP"),
        end_time_requested ? end_time_env : nullptr, time_mult);
    if (!cap_sel.valid()) {
        std::cerr << "step-cap error: " << cap_sel.error << std::endl;
        return 1;
    }
    const long long step_cap = cap_sel.cap;

    std::ofstream logf("outputs/output.log", std::ios::app);
    logf << "[" << out_path << " mode=" << mode << "] Total time = " << total_time << std::endl;
    std::cout << "[" << mode << "] Total time = " << total_time << std::endl;

    // Multi-snapshot output. Format:
    //   line 1: "ns num_of_eq"
    //   line 2: ns cumulative cell heights in km (offset by grid.out_base_km, the
    //           domain base — 1003 km for C7-based scenarios, 0 km for the
    //           photosphere-anchored model_column C7 IC; ds_i is in m, hence
    //           the 1e-3 conversion)
    //   then, repeated: a "# t = T step = S" marker followed by ns lines of
    //   num_of_eq space-separated conserved-variable values.
    const arma::uword state_rows = release_mode ? num_of_mixture_eq : num_of_eq;
    std::ofstream fout;
    if (write_output) fout.open(out_path);
    // NOTE: a default-constructed std::ofstream has goodbit, so `if (stream)` is
    // TRUE even when nothing was ever opened. Always test is_open().
    if (fout.is_open()) {
        fout << grid.ns << " " << state_rows << '\n';
        float cum_km = 0.0f;
        for (arma::uword i = 0; i < grid.ns; ++i) {
            cum_km += grid.ds_i(i) * 1.0e-3f;
            fout << "  " << (cum_km + grid.out_base_km);
        }
        fout << '\n';
    }
    if (fout.is_open() && release_mode && write_gamma_diag)
        fout << "# EOS_MODE=gamma_table diagnostics=" << out_path << ".gamma_diag\n";
    std::ofstream gamma_diag;
    if (release_mode && write_gamma_diag) {
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
                   << "# release conduction is PHYSICAL only, so kappa_solver is\n"
                      "# identically kappa_physical; the column is kept for format\n"
                      "# stability with existing analysis scripts.\n"
                   << "# columns=rho_total v_cm T x_eq n_e n_HI p_total Gamma1 kappa_physical kappa_solver\n";
    }

    auto write_frame = [&](float t_now, long long step_now) {
        ProfileScope output_timer(ProfileRegion::Output);
        if (fout.is_open()) {
            fout << "# t = " << t_now << " step = " << step_now << '\n';
            for (arma::uword i = 0; i < grid.ns; ++i) {
                for (arma::uword k = 0; k < state_rows; ++k) {
                    fout << "  " << xn(arma::sub2ind(
                        arma::size(grid.ns,state_rows),i,k));
                }
                fout << '\n';
            }
        }
        if (gamma_diag.is_open()) {
            gamma_diag << "# t = " << t_now << " step = " << step_now << '\n';
            const auto sz = arma::size(grid.ns, num_of_mixture_eq);
            for (arma::uword i = 0; i < grid.ns; ++i) {
                auto at = [&](arma::uword k) -> double {
                    return static_cast<double>(xn(arma::sub2ind(sz, i, k)));
                };
                const MixtureThermo th = decode_equilibrium_mixture(
                    grid.eos_gamma_table, at(mix::RHO), at(mix::MOM),
                    at(mix::ENERGY), mixture_cell_phi(grid, i),
                    std::numeric_limits<double>::quiet_NaN(), grid.eos_gamma_debug_clamp);
                const double velocity = at(mix::MOM)/th.rho;
                gamma_diag << "  " << th.rho << "  " << velocity << "  " << th.T
                           << "  " << th.x << "  " << th.n_e << "  " << th.n_HI
                           << "  " << th.p << "  " << th.gamma1
                           << "  " << physical_conductivity(th.n_e,th.n_HI,th.T)
                           << "  " << physical_conductivity(th.n_e,th.n_HI,th.T)
                           << '\n';
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
    // CHROMO_FRAME_DT switches the snapshot cadence from a step stride to
    // physical seconds. Required for long CHROMO_T_END runs: the stride above is
    // derived from time_mult, so a 1000 s run at a 20 s stride would emit tens of
    // thousands of snapshots and spend most of its wall time formatting ASCII.
    // Unset -> the legacy stride cadence is reproduced exactly.
    OutputSchedule schedule = OutputSchedule::from_stride(frame_stride);
    int progress_stride = 100;
    if (const char* e = std::getenv("CHROMO_PROGRESS_STRIDE")) {
        try { progress_stride = std::max(1, std::stoi(e)); } catch (...) {}
    }
    if (const char* e = std::getenv("CHROMO_FRAME_DT")) {
        double frame_dt = 0.0;
        try { frame_dt = std::stod(e); } catch (...) {}
        if (!(frame_dt > 0.0) || !std::isfinite(frame_dt)) {
            std::cerr << "CHROMO_FRAME_DT must be a positive finite number of "
                         "physical seconds, got '" << e << "'" << std::endl;
            return 1;
        }
        schedule = OutputSchedule::from_frame_dt(frame_dt);
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
        if (!release_mode)
            throw std::runtime_error(
                "CHROMO_FACE_FLUX_DIAG requires the release solver (the capture "
                "lives in mixture_rhs_explicit)");
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
        face_flux << "# face-flux diagnostic (read-only capture of mixture_rhs_explicit)\n"
                  << "# ns=" << grid.ns << " first_cell=" << face_flux_lo
                  << " uniform_mesh=" << grid.uniform_mesh
                  << " mc3=" << grid.mc3_limiter << " beta=" << grid.limiter_beta
                  << " roe=" << grid.roe_characteristic_flux
                  << " swmf_godunov=" << grid.swmf_godunov_flux
                  << " pressure_reconstruct=" << grid.pressure_reconstruct
                  << " ds_km=" << (grid.ds_i(0)*1.0e-3f)
                  << "\n# mass flux = the mix::RHO continuity row (total mass "
                     "density; the release state has no carrier rows)\n"
                  << "# columns=cell_km face_km rho_cell v_cell T_cell "
                     "rho_L rho_R v_L v_R T_L T_R cs_L cs_R a_face "
                     "f_central f_diff f_total "
                     "r_rho phi_plus_rho r_ip1_rho phi_minus_rho "
                     "r_v phi_plus_v r_T phi_plus_T "
                     "p_cell p_L p_R\n";
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
        if (!release_mode)
            throw std::runtime_error(
                "CHROMO_OUTER_COND_DIAG requires the release solver (the capture "
                "lives in mixture_apply_conduction)");
        outer_cond.open(out_path+".outercond");
        if (!outer_cond) throw std::runtime_error("cannot open outer-conduction sidecar");
        outer_cond << "# outer-face conduction diagnostic (read-only capture of the "
                      "final converged mixture_apply_conduction iteration)\n"
                   << "# ns=" << grid.ns << " ds_km=" << (grid.ds_i(grid.ns-1)*1.0e-3f)
                   << " uniform_mesh=" << grid.uniform_mesh
                   << " impose_outer_heat_flux=" << grid.impose_outer_heat_flux << '\n'
                   << "# release conduction is PHYSICAL only (kappa_e + kappa_n; no "
                      "TRAC, no artificial diffusivity)\n"
                   << "# ds_face = half the top cell when an external face "
                      "temperature is imposed, else the centre-to-ghost-centre span\n"
                   << "# q_face = kappa_face*(T_wall - T_top)/ds_face"
                      "  [W m^-2, positive = into the top cell]\n"
                   << "# columns=t step T_top T_wall kappa_face "
                      "ds_face area_ratio q_face\n";
        outer_cond.precision(10);
    }

    auto write_outer_cond_record = [&](float t_now, long long step_now) {
        const OuterConductionCapture& oc = grid.outer_conduction_capture;
        if (!outer_cond.is_open() || !oc.valid) return;
        outer_cond << t_now << ' ' << step_now << ' '
                   << oc.T_top << ' ' << oc.T_wall << ' '
                   << oc.kappa_face << ' '
                   << oc.ds_face << ' ' << oc.area_ratio << ' '
                   << oc.q_face << '\n';
    };

    auto write_face_record = [&](float t_now, long long step_now) {
        const MixtureFaceFluxCapture& c = grid.face_flux_capture;
        if (!face_flux.is_open() || !c.valid) return;
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
                    << c.r_rho[i] << ' ' << c.phi_plus_rho[i] << ' '
                    << c.r_ip1_rho[i] << ' ' << c.phi_minus_rho[i] << ' '
                    << c.r_v[i] << ' ' << c.phi_plus_v[i] << ' '
                    << c.r_T[i] << ' ' << c.phi_plus_T[i] << ' '
                    << c.p_cell[i] << ' ' << c.p_L[i] << ' ' << c.p_R[i] << '\n';
            }
        }
    };

    float     time         = 0.0f;
    long long step         = 0;
    std::uint64_t state_generation = 0;
    MixtureField decoded_storage[2];
    MixtureField* decoded=&decoded_storage[0];
    const MixtureField* previous_decoded=nullptr;
    if (schedule.due(step, time)) {
        write_frame(time, step);
        schedule.note_written(step, time);
    }

    while (time < total_time && step < step_cap) {
        grid.sim_time = time;   // expose current time to time-dependent terms (beam window)
        {
            ProfileScope boundary_timer(ProfileRegion::Boundary);
            sc.update_bc(grid, xn);
        }
        Vec dt;
        if (release_mode) {
            mixture_decode_into(grid,xn,*decoded,state_generation,previous_decoded);
            dt = mixture_timestep(grid,xn,*decoded);
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

        if (step % progress_stride == 0) {
            std::cout << "[" << mode << "] step = " << step
                      << "  dt = " << dt_avg
                      << "  time = " << time
                      << "  T_c = " << grid.trac_cutoff_T << std::endl;
        }

        if (release_mode)
            xn = mixture_advance(grid,xn,dt,*decoded);
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
        if (release_mode) {
            previous_decoded=decoded;
            decoded = decoded == &decoded_storage[0]
                ? &decoded_storage[1] : &decoded_storage[0];
        }

        if (schedule.due(step, time)) {
            write_frame(time, step);
            schedule.note_written(step, time);
        }
    }
    // Always capture the final state, even when it does not land on the cadence.
    if (!schedule.already_written(step)) {
        write_frame(time, step);
        schedule.note_written(step, time);
    }

    std::cout << "[" << mode << "] END! step=" << step << " time=" << time << std::endl;
    const Termination termination =
        classify_termination(time, total_time, step, step_cap);
    std::cout << "termination=" << termination_name(termination) << '\n'
              << "requested_end_time=" << total_time << '\n'
              << "final_physical_time=" << time << '\n'
              << "final_step=" << step << '\n'
              << "effective_step_cap=" << step_cap
              << " (source=" << step_cap_source_name(cap_sel.source) << ")" << std::endl;
    if (termination == Termination::StepCap) {
        std::cerr << "WARNING: run stopped at the step cap (" << step_cap
                  << ") after " << time << " of the requested " << total_time
                  << " physical seconds. Raise or unset CHROMO_STEP_CAP." << std::endl;
    }
    // Release swmf_godunov flux: report what the exact Riemann solver did, so a
    // silent fallback to Rusanov can never be mistaken for a Godunov result.
    if (grid.swmf_godunov_flux) {
        const GodunovFluxStats& g = grid.godunov_stats;
        std::cout << "godunov.faces=" << g.faces << '\n'
                  << "godunov.exact=" << g.exact << '\n'
                  << "godunov.fallbacks=" << g.fallbacks()
                  << " (bad_input=" << g.fallback_bad_input
                  << " vacuum=" << g.fallback_vacuum
                  << " negative_p=" << g.fallback_negative_p
                  << " no_converge=" << g.fallback_no_converge
                  << " bad_sample=" << g.fallback_bad_sample << ")\n"
                  << "godunov.iterations_max=" << g.max_iterations
                  << " mean="
                  << (g.exact ? static_cast<double>(g.total_iterations)
                                / static_cast<double>(g.exact) : 0.0)
                  << std::endl;
        if (g.fallbacks() != 0)
            std::cerr << "WARNING: the swmf-godunov flux fell back to Rusanov on "
                      << g.fallbacks() << " of " << g.faces << " face solves."
                      << std::endl;
    }
    print_runtime_profile(std::cout);
    if (eos_counting_on) {
        const EosOperationCounts counts = eos_operation_counts();
        std::cout << "eos.gamma1_queries=" << counts.gamma1_queries << '\n'
                  << "eos.temperature_logs=" << counts.temperature_logs << '\n'
                  << "eos.n_h_logs=" << counts.n_h_logs << std::endl;
    }
    return 0;
}
