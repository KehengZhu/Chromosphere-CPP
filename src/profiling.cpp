#include "profiling.hpp"
#include "parallel.hpp"

#include <algorithm>
#include <array>
#include <iomanip>
#include <ostream>

namespace chromosphere {
namespace {

bool enabled = false;

struct alignas(64) ThreadProfileSlot {
    std::array<double, static_cast<unsigned>(ProfileRegion::Count)> seconds{};
    std::array<std::uint64_t, static_cast<unsigned>(TimestepLimiter::Count)> limiters{};
    EosInversionProfile inversions;
    std::uint64_t conduction_calls = 0;
    std::uint64_t conduction_total_iterations = 0;
    std::uint64_t conduction_maximum_iterations = 0;
};

std::array<ThreadProfileSlot, kMaximumParallelThreads> slots;

ThreadProfileSlot& current_slot() noexcept {
    return slots[static_cast<std::size_t>(parallel_thread_index())];
}

ThreadProfileSlot merged_profile() {
    ThreadProfileSlot merged;
    for (int tid = 0; tid < kMaximumParallelThreads; ++tid) {
        const ThreadProfileSlot& slot = slots[static_cast<std::size_t>(tid)];
        for (unsigned i = 0; i < static_cast<unsigned>(ProfileRegion::Count); ++i)
            merged.seconds[i] += slot.seconds[i];
        for (unsigned i = 0; i < static_cast<unsigned>(TimestepLimiter::Count); ++i)
            merged.limiters[i] += slot.limiters[i];
        merged.inversions.calls += slot.inversions.calls;
        merged.inversions.initial_guess_accepts += slot.inversions.initial_guess_accepts;
        merged.inversions.one_update_convergences += slot.inversions.one_update_convergences;
        merged.inversions.total_iterations += slot.inversions.total_iterations;
        merged.inversions.maximum_iterations = std::max(
            merged.inversions.maximum_iterations, slot.inversions.maximum_iterations);
        merged.inversions.bisection_fallbacks += slot.inversions.bisection_fallbacks;
        merged.inversions.bracket_evaluations += slot.inversions.bracket_evaluations;
        merged.conduction_calls += slot.conduction_calls;
        merged.conduction_total_iterations += slot.conduction_total_iterations;
        merged.conduction_maximum_iterations = std::max(
            merged.conduction_maximum_iterations, slot.conduction_maximum_iterations);
    }
    return merged;
}

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
    require_supported_parallel_runtime();
    for (auto& slot : slots) slot = ThreadProfileSlot{};
}

void profile_add_seconds(ProfileRegion region, double value) noexcept {
    if (enabled) current_slot().seconds[static_cast<unsigned>(region)] += value;
}

void profile_note_timestep_limiter(TimestepLimiter limiter) noexcept {
    if (enabled) ++current_slot().limiters[static_cast<unsigned>(limiter)];
}

void profile_note_inversion_call() noexcept {
    if (enabled) ++current_slot().inversions.calls;
}

void profile_note_inversion_initial_guess() noexcept {
    if (enabled) ++current_slot().inversions.initial_guess_accepts;
}

void profile_note_inversion_iterations(std::uint64_t count) noexcept {
    if (!enabled) return;
    EosInversionProfile& inversions = current_slot().inversions;
    inversions.total_iterations += count;
    inversions.maximum_iterations = std::max(inversions.maximum_iterations, count);
    if (count == 1) ++inversions.one_update_convergences;
}

void profile_note_inversion_bisection() noexcept {
    if (enabled) ++current_slot().inversions.bisection_fallbacks;
}

void profile_note_inversion_bracket_evaluations(std::uint64_t count) noexcept {
    if (enabled) current_slot().inversions.bracket_evaluations += count;
}

EosInversionProfile eos_inversion_profile() { return merged_profile().inversions; }

void profile_note_conduction_iterations(std::uint64_t count) noexcept {
    if (!enabled) return;
    ThreadProfileSlot& slot = current_slot();
    ++slot.conduction_calls;
    slot.conduction_total_iterations += count;
    slot.conduction_maximum_iterations = std::max(slot.conduction_maximum_iterations, count);
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
    const ThreadProfileSlot merged = merged_profile();
    out << "\n===== CHROMO runtime profile =====\n";
    out << std::fixed << std::setprecision(6);
    for (unsigned i = 0; i < static_cast<unsigned>(ProfileRegion::Count); ++i)
        out << "time." << region_name(static_cast<ProfileRegion>(i)) << "_s="
            << merged.seconds[i] << '\n';
    for (unsigned i = 0; i < static_cast<unsigned>(TimestepLimiter::Count); ++i)
        out << "limiter." << limiter_name(static_cast<TimestepLimiter>(i)) << '='
            << merged.limiters[i] << '\n';
    const EosInversionProfile& inversions = merged.inversions;
    const double average = inversions.calls
        ? static_cast<double>(inversions.total_iterations)/inversions.calls : 0.0;
    out << "inversion.calls=" << inversions.calls << '\n'
        << "inversion.initial_guess_accepts=" << inversions.initial_guess_accepts << '\n'
        << "inversion.one_update_convergences=" << inversions.one_update_convergences << '\n'
        << "inversion.average_iterations=" << average << '\n'
        << "inversion.maximum_iterations=" << inversions.maximum_iterations << '\n'
        << "inversion.bisection_fallbacks=" << inversions.bisection_fallbacks << '\n'
        << "inversion.bracket_evaluations=" << inversions.bracket_evaluations << '\n'
        << "conduction.calls=" << merged.conduction_calls << '\n'
        << "conduction.average_newton_updates="
        << (merged.conduction_calls
            ? static_cast<double>(merged.conduction_total_iterations)/merged.conduction_calls
            : 0.0) << '\n'
        << "conduction.maximum_newton_updates=" << merged.conduction_maximum_iterations << '\n'
        << "==================================\n";
}

} // namespace chromosphere
