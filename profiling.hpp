#pragma once

#include <chrono>
#include <cstdint>
#include <iosfwd>

namespace chromosphere {

enum class ProfileRegion : unsigned {
    Cfl,
    Decode,
    Rhs,
    Projection,
    Conduction,
    Boundary,
    Output,
    Count
};

enum class TimestepLimiter : unsigned {
    Acoustic,
    BeamHeating,
    CoronalHeating,
    Cooling,
    ElectronEnergy,
    Other,
    Count
};

struct EosInversionProfile {
    std::uint64_t calls = 0;
    std::uint64_t initial_guess_accepts = 0;
    std::uint64_t one_update_convergences = 0;
    std::uint64_t total_iterations = 0;
    std::uint64_t maximum_iterations = 0;
    std::uint64_t bisection_fallbacks = 0;
    std::uint64_t bracket_evaluations = 0;
};

void set_runtime_profiling(bool enabled);
bool runtime_profiling_enabled() noexcept;
void reset_runtime_profile();
void print_runtime_profile(std::ostream& out);

void profile_add_seconds(ProfileRegion region, double seconds) noexcept;
void profile_note_timestep_limiter(TimestepLimiter limiter) noexcept;
void profile_note_inversion_call() noexcept;
void profile_note_inversion_initial_guess() noexcept;
void profile_note_inversion_iterations(std::uint64_t iterations) noexcept;
void profile_note_inversion_bisection() noexcept;
void profile_note_inversion_bracket_evaluations(std::uint64_t count) noexcept;
void profile_note_conduction_iterations(std::uint64_t iterations) noexcept;
EosInversionProfile eos_inversion_profile();

class ProfileScope {
public:
    explicit ProfileScope(ProfileRegion region) noexcept;
    ~ProfileScope();
    ProfileScope(const ProfileScope&) = delete;
    ProfileScope& operator=(const ProfileScope&) = delete;
private:
    ProfileRegion region_;
    bool enabled_;
    std::chrono::steady_clock::time_point start_;
};

} // namespace chromosphere
