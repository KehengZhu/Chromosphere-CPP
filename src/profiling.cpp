#include "profiling.hpp"

#include <algorithm>
#include <array>
#include <iomanip>
#include <ostream>

namespace chromosphere {
namespace {

bool enabled = false;
std::array<double, static_cast<unsigned>(ProfileRegion::Count)> seconds{};
std::array<std::uint64_t, static_cast<unsigned>(TimestepLimiter::Count)> limiters{};
EosInversionProfile inversions;
std::uint64_t conduction_calls=0, conduction_total_iterations=0;
std::uint64_t conduction_maximum_iterations=0;

const char* region_name(ProfileRegion region) {
    switch (region) {
    case ProfileRegion::Cfl: return "cfl";
    case ProfileRegion::Decode: return "decode";
    case ProfileRegion::Rhs: return "rhs";
    case ProfileRegion::Projection: return "projection";
    case ProfileRegion::Conduction: return "conduction";
    case ProfileRegion::Boundary: return "boundary";
    case ProfileRegion::Output: return "output";
    case ProfileRegion::Count: break;
    }
    return "unknown";
}

const char* limiter_name(TimestepLimiter limiter) {
    switch (limiter) {
    case TimestepLimiter::Acoustic: return "acoustic";
    case TimestepLimiter::BeamHeating: return "beam_heating";
    case TimestepLimiter::CoronalHeating: return "coronal_heating";
    case TimestepLimiter::Cooling: return "cooling";
    case TimestepLimiter::ElectronEnergy: return "electron_energy";
    case TimestepLimiter::Other: return "other";
    case TimestepLimiter::Count: break;
    }
    return "unknown";
}

} // namespace

void set_runtime_profiling(bool value) { enabled = value; }
bool runtime_profiling_enabled() noexcept { return enabled; }

void reset_runtime_profile() {
    seconds.fill(0.0);
    limiters.fill(0);
    inversions = EosInversionProfile{};
    conduction_calls=conduction_total_iterations=conduction_maximum_iterations=0;
}

void profile_add_seconds(ProfileRegion region, double value) noexcept {
    if (enabled) seconds[static_cast<unsigned>(region)] += value;
}

void profile_note_timestep_limiter(TimestepLimiter limiter) noexcept {
    if (enabled) ++limiters[static_cast<unsigned>(limiter)];
}

void profile_note_inversion_call() noexcept {
    if (enabled) ++inversions.calls;
}

void profile_note_inversion_initial_guess() noexcept {
    if (enabled) ++inversions.initial_guess_accepts;
}

void profile_note_inversion_iterations(std::uint64_t count) noexcept {
    if (!enabled) return;
    inversions.total_iterations += count;
    inversions.maximum_iterations = std::max(inversions.maximum_iterations, count);
    if (count == 1) ++inversions.one_update_convergences;
}

void profile_note_inversion_bisection() noexcept {
    if (enabled) ++inversions.bisection_fallbacks;
}

void profile_note_inversion_bracket_evaluations(std::uint64_t count) noexcept {
    if (enabled) inversions.bracket_evaluations += count;
}

EosInversionProfile eos_inversion_profile() { return inversions; }

void profile_note_conduction_iterations(std::uint64_t count) noexcept {
    if (!enabled) return;
    ++conduction_calls;
    conduction_total_iterations+=count;
    conduction_maximum_iterations=std::max(conduction_maximum_iterations,count);
}

ProfileScope::ProfileScope(ProfileRegion region) noexcept
    : region_(region), enabled_(enabled) {
    if (enabled_) start_ = std::chrono::steady_clock::now();
}

ProfileScope::~ProfileScope() {
    if (!enabled_) return;
    const auto stop = std::chrono::steady_clock::now();
    profile_add_seconds(region_, std::chrono::duration<double>(stop-start_).count());
}

void print_runtime_profile(std::ostream& out) {
    if (!enabled) return;
    out << "\n===== CHROMO runtime profile =====\n";
    out << std::fixed << std::setprecision(6);
    for (unsigned i = 0; i < static_cast<unsigned>(ProfileRegion::Count); ++i)
        out << "time." << region_name(static_cast<ProfileRegion>(i)) << "_s="
            << seconds[i] << '\n';
    for (unsigned i = 0; i < static_cast<unsigned>(TimestepLimiter::Count); ++i)
        out << "limiter." << limiter_name(static_cast<TimestepLimiter>(i)) << '='
            << limiters[i] << '\n';
    const double average = inversions.calls
        ? static_cast<double>(inversions.total_iterations)/inversions.calls : 0.0;
    out << "inversion.calls=" << inversions.calls << '\n'
        << "inversion.initial_guess_accepts=" << inversions.initial_guess_accepts << '\n'
        << "inversion.one_update_convergences=" << inversions.one_update_convergences << '\n'
        << "inversion.average_iterations=" << average << '\n'
        << "inversion.maximum_iterations=" << inversions.maximum_iterations << '\n'
        << "inversion.bisection_fallbacks=" << inversions.bisection_fallbacks << '\n'
        << "inversion.bracket_evaluations=" << inversions.bracket_evaluations << '\n'
        << "conduction.calls=" << conduction_calls << '\n'
        << "conduction.average_newton_updates="
        << (conduction_calls ? static_cast<double>(conduction_total_iterations)/conduction_calls : 0.0) << '\n'
        << "conduction.maximum_newton_updates=" << conduction_maximum_iterations << '\n'
        << "==================================\n";
}

} // namespace chromosphere
