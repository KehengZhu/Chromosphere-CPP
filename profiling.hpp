#pragma once
/*!
 * @file profiling.hpp
 * @brief Opt-in runtime profiling: region timers, EOS-inversion and conduction
 *        iteration statistics, and the active timestep limiter.
 * @ingroup runtime
 *
 * Everything here is disabled unless the driver calls set_runtime_profiling(true)
 * (`CHROMO_PROFILE=1`), so production hot paths pay only a predictable false
 * branch. Profiling never alters solver state or the timestep.
 */

#include <chrono>
#include <cstdint>
#include <iosfwd>

namespace chromosphere {

/** @addtogroup runtime
 *  @{
 */

/// Wall-clock timer buckets. One ProfileScope in each stage attributes the
/// seconds it was alive to exactly one of these; the report prints them as
/// `time.<name>_s`. The buckets are stage labels, not a partition of the step:
/// they neither overlap nor necessarily cover the whole timestep.
enum class ProfileRegion : unsigned {
    Cfl,         ///< Timestep sizing (mixture_timestep / the legacy cal_dt_i).
    Decode,      ///< Conserved-to-thermodynamic decode (mixture_decode_into).
    Rhs,         ///< Explicit MUSCL-Hancock/Roe hydro RHS (mixture_rhs_explicit).
    Conduction,  ///< Implicit backward-Euler conduction (mixture_apply_conduction).
    Boundary,    ///< Scenario boundary-condition update at the top of each step.
    Output,      ///< Formatting and writing one ASCII snapshot frame.
    Count        ///< Number of regions; array bound, never a timed region.
};

/// Which physical constraint set the timestep actually taken. Recorded once per
/// timestep computation via profile_note_timestep_limiter() and reported as
/// `limiter.<name>` counts, so a run shows what was limiting it and for how many
/// steps. The release solver always reports Acoustic — with no volumetric source
/// and an unconditionally stable implicit conduction stage, the acoustic CFL
/// condition is its only constraint. The remaining values come from the legacy
/// two-fluid timestep, which lowers dt when a source term's heating time is
/// shorter than the acoustic step.
enum class TimestepLimiter : unsigned {
    Acoustic,        ///< Acoustic CFL, CFL*ds/(|u|+c_s). The only release limiter.
    BeamHeating,     ///< Flare beam energy-gain cap (two-fluid, enable_beam_heating).
    CoronalHeating,  ///< Coronal heating energy-gain cap (two-fluid, enable_coronal_heating).
    Cooling,         ///< Radiative-cooling category; no current rule selects it.
    ElectronEnergy,  ///< Electron-energy category; no current rule selects it.
    Other,           ///< Catch-all category; no current rule selects it.
    Count            ///< Number of limiters; array bound, never a recorded value.
};

/// Counters for the (rho,p) -> T face inversion of the pressure-based MUSCL
/// reconstruction ONLY. Deliberately separate from EosInversionProfile, which
/// counts the (rho,E) -> T energy inversion of the decode/projection stages.
struct PressureInversionProfile {
    /// Safeguarded Newton/bisection solves that returned a temperature. One is
    /// recorded per successful inversion; a solve that exhausts its 128-iteration
    /// budget throws instead and is never counted.
    std::uint64_t calls = 0;
    std::uint64_t hinted_calls = 0;  ///< a caller-supplied guess was usable
    std::uint64_t evaluations = 0;  ///< Saha/caloric evaluations consumed
    /// Largest number of Saha/caloric evaluations any single solve needed — a
    /// maximum over calls (and over threads), not a sum, so it exposes the worst
    /// face rather than the typical one.
    std::uint64_t maximum_evaluations = 0;
    /// Newton steps that were discarded because the candidate temperature fell
    /// outside the current bracket and had to be replaced by the bracket
    /// midpoint. Summed over steps, so it can exceed #calls.
    std::uint64_t bisection_fallbacks = 0;
};

/// Counters for the (rho,E) -> T ENERGY inversion, i.e. the safeguarded
/// Newton/bisection solve of e_int(rho,T) = e_int that the decode and projection
/// stages run per cell. Deliberately separate from PressureInversionProfile,
/// which counts the (rho,p) -> T face inversion of the MUSCL reconstruction.
struct EosInversionProfile {
    /// Energy inversions attempted: one per call that passes argument validation,
    /// counted before any fast path, so it is the denominator for the averages in
    /// the report.
    std::uint64_t calls = 0;
    /// Calls whose supplied temperature guess already satisfied the energy
    /// tolerance, so the root was returned with zero Newton updates.
    std::uint64_t initial_guess_accepts = 0;
    /// Calls that converged after exactly one Newton update — the measured common
    /// case on the production column, and what the fused fast path targets.
    std::uint64_t one_update_convergences = 0;
    /// Sum over calls of the Newton updates each one reported; the report divides
    /// it by #calls to print `inversion.average_iterations`.
    std::uint64_t total_iterations = 0;
    /// Largest per-call Newton-update count — a maximum over calls (and over
    /// threads at merge time), not a sum.
    std::uint64_t maximum_iterations = 0;
    /// Newton steps rejected by a safeguard (candidate outside the bracket, or
    /// progress slower than halving the previous bracket) and replaced by the
    /// bracket midpoint. Summed over steps, so it can exceed #calls.
    std::uint64_t bisection_fallbacks = 0;
    /// Caloric evaluations spent establishing the table-wide bracket: two per
    /// call (the minimum and maximum table temperatures), added only when both
    /// fast paths miss. An accepted fast path deliberately contributes zero, so
    /// this counter measures how often the full solver was actually entered.
    std::uint64_t bracket_evaluations = 0;
};

/// Turn profiling on or off globally. The driver enables it for `CHROMO_PROFILE=1`
/// and nothing else; with it off, every profile_* entry point below is a predicted
/// false branch. Profiling never alters solver state or the timestep.
void set_runtime_profiling(bool enabled);
/// Whether profiling is currently on (see set_runtime_profiling()).
bool runtime_profiling_enabled() noexcept;
/// Clear every per-thread counter slot back to zero. Also validates that the
/// OpenMP thread count fits the fixed slot capacity, so it must be called from
/// serial code before any profiled parallel region runs.
void reset_runtime_profile();
/// Write the merged report (region seconds, limiter counts, both inversion
/// profiles, conduction Newton statistics) as `key=value` lines. Does nothing
/// when profiling is off.
void print_runtime_profile(std::ostream& out);

/// Add elapsed wall-clock seconds to one region bucket of the calling thread.
/// Normally reached through ProfileScope rather than called directly.
void profile_add_seconds(ProfileRegion region, double seconds) noexcept;
/// Record which constraint set the timestep for one timestep computation.
void profile_note_timestep_limiter(TimestepLimiter limiter) noexcept;
/// Count one attempted (rho,E) -> T energy inversion
/// (EosInversionProfile::calls).
void profile_note_inversion_call() noexcept;
/// Record that an energy inversion was answered by its supplied temperature guess
/// with no Newton update (EosInversionProfile::initial_guess_accepts).
void profile_note_inversion_initial_guess() noexcept;
/// Record how many Newton updates one energy inversion consumed: accumulates the
/// total, updates the per-call maximum, and counts a value of exactly one as a
/// one-update convergence.
void profile_note_inversion_iterations(std::uint64_t iterations) noexcept;
/// Record one Newton step of the energy inversion that a safeguard replaced with
/// a bracket bisection.
void profile_note_inversion_bisection() noexcept;
/// Add table-endpoint caloric evaluations charged to the energy inversion's
/// bracket setup (the solver passes 2, one per endpoint).
void profile_note_inversion_bracket_evaluations(std::uint64_t count) noexcept;
/// Record one implicit conduction solve and the number of Newton updates it took
/// to converge; feeds `conduction.calls`, the average, and the maximum.
void profile_note_conduction_iterations(std::uint64_t iterations) noexcept;
/// Snapshot of the energy-inversion counters summed over all thread slots.
EosInversionProfile eos_inversion_profile();
/// Record one completed (rho,p) -> T face inversion: the Saha/caloric evaluations
/// it consumed, whether a caller-supplied temperature hint was usable as the
/// starting point, and how many of its Newton steps fell back to bisection.
void profile_note_pressure_inversion(std::uint64_t evaluations, bool hinted,
                                     std::uint64_t bisections) noexcept;
/// Snapshot of the pressure-inversion counters summed over all thread slots.
PressureInversionProfile pressure_inversion_profile();

/// RAII wall-clock timer for one ProfileRegion: it samples a steady_clock stamp
/// at construction and, at destruction, adds the elapsed seconds to that region's
/// bucket for the calling thread. Whether profiling is on is latched at
/// construction, so toggling it mid-scope cannot produce a one-sided measurement.
/// Non-copyable, and a no-op beyond one branch when profiling is off.
class ProfileScope {
public:
    /// Start timing @p region (no clock read at all when profiling is off).
    explicit ProfileScope(ProfileRegion region) noexcept;
    /// Stop timing and charge the elapsed seconds to the region.
    ~ProfileScope();
    ProfileScope(const ProfileScope&) = delete;
    ProfileScope& operator=(const ProfileScope&) = delete;
private:
    ProfileRegion region_;  ///< bucket the elapsed time is charged to
    bool enabled_;  ///< profiling state latched at construction
    std::chrono::steady_clock::time_point start_;  ///< set only when enabled_
};

/** @} */

} // namespace chromosphere
