#pragma once
// Run-control helpers for the driver loop: step-cap selection and snapshot
// scheduling. Both are pure, header-only, and free of OpenMP/Armadillo so they
// can be unit tested directly without running a simulation.
//
// Step cap (CHROMO_STEP_CAP / CHROMO_T_END)
//   Historically the cap was derived from the positional time_mult argument
//   (max(10000, 10000*time_mult)). That is correct for legacy time_mult-driven
//   runs but silently truncates a run whose end time comes from CHROMO_T_END:
//   CHROMO_T_END=1000 with time_mult=20 caps at 200000 steps, roughly 140 of
//   the 1000 requested physical seconds. select_step_cap() keeps the legacy
//   behavior only when CHROMO_T_END is absent.
//
// Snapshot cadence (CHROMO_FRAME_DT)
//   The legacy cadence is a step stride, also derived from time_mult, so a long
//   CHROMO_T_END run keeps a short-run stride and can emit tens of thousands of
//   snapshots. OutputSchedule adds a physical-time cadence; with CHROMO_FRAME_DT
//   unset it reproduces the stride behavior exactly.

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>

namespace chromosphere {

// Effectively unbounded, but far below LLONG_MAX so neither the cap nor the
// step counter can overflow. At 1 ms/step this is ~31 years of wall time.
constexpr long long kUnboundedStepCap = 1000000000000LL; // 1e12

enum class StepCapSource {
    Explicit,       // CHROMO_STEP_CAP
    LegacyTimeMult, // max(10000, 10000*time_mult), no CHROMO_T_END
    EndTime,        // CHROMO_T_END set, no explicit cap -> unbounded
    Invalid
};

struct StepCapSelection {
    long long     cap    = kUnboundedStepCap;
    StepCapSource source = StepCapSource::EndTime;
    std::string   error;                        // non-empty iff source == Invalid
    bool valid() const { return source != StepCapSource::Invalid; }
};

// step_cap_env / end_time_env are the raw getenv() results (nullptr when unset).
// An explicit CHROMO_STEP_CAP always wins; an unparseable or non-positive value
// is rejected rather than silently replaced by a dangerous cap.
inline StepCapSelection select_step_cap(const char* step_cap_env,
                                        const char* end_time_env,
                                        float time_mult) {
    StepCapSelection out;
    if (step_cap_env && *step_cap_env) {
        errno = 0;
        char* tail = nullptr;
        const long long requested = std::strtoll(step_cap_env, &tail, 10);
        while (tail && (*tail == ' ' || *tail == '\t' || *tail == '\n')) ++tail;
        if (errno == ERANGE || !tail || *tail != '\0') {
            out.source = StepCapSource::Invalid;
            out.error  = "CHROMO_STEP_CAP is not a valid integer: '"
                       + std::string(step_cap_env) + "'";
            return out;
        }
        if (requested <= 0) {
            out.source = StepCapSource::Invalid;
            out.error  = "CHROMO_STEP_CAP must be a positive integer, got '"
                       + std::string(step_cap_env) + "'";
            return out;
        }
        out.cap    = requested < kUnboundedStepCap ? requested : kUnboundedStepCap;
        out.source = StepCapSource::Explicit;
        return out;
    }
    if (end_time_env && *end_time_env) {
        // Explicit end time, no explicit cap: the end time governs.
        out.cap    = kUnboundedStepCap;
        out.source = StepCapSource::EndTime;
        return out;
    }
    // Legacy path, bit-for-bit the historical expression, guarded against a
    // float overflow / out-of-range cast for absurd time_mult values.
    const float legacy = std::max(10000.0f, 10000.0f * time_mult);
    out.cap = (!std::isfinite(legacy) || legacy >= 1.0e12f)
        ? kUnboundedStepCap : static_cast<long long>(legacy);
    out.source = StepCapSource::LegacyTimeMult;
    return out;
}

inline const char* step_cap_source_name(StepCapSource s) {
    switch (s) {
        case StepCapSource::Explicit:       return "explicit_step_cap";
        case StepCapSource::LegacyTimeMult: return "legacy_time_mult";
        case StepCapSource::EndTime:        return "end_time";
        default:                            return "invalid";
    }
}

enum class Termination { EndTime, StepCap, Other };

inline const char* termination_name(Termination t) {
    switch (t) {
        case Termination::EndTime: return "end_time";
        case Termination::StepCap: return "step_cap";
        default:                   return "other";
    }
}

// Reached end time wins over the cap when both fire on the same step.
inline Termination classify_termination(double final_time, double requested_end_time,
                                        long long final_step, long long cap) {
    if (final_time >= requested_end_time) return Termination::EndTime;
    if (final_step >= cap)                return Termination::StepCap;
    return Termination::Other;
}

// Snapshot scheduler. Two modes:
//   stride mode (legacy):  write when step % stride == 0
//   time mode:             write when the simulation time reaches
//                          frame_index * frame_dt
// The scheduled time is always recomputed from an integer frame index, so the
// cadence cannot drift no matter how many intervals elapse. A step that crosses
// several nominal output times still writes the state only once.
struct OutputSchedule {
    bool      time_based       = false;
    double    frame_dt         = 0.0;
    long long frame_index      = 0;   // index of the NEXT scheduled output time
    long long stride           = 1;
    long long last_written_step = -1;

    static OutputSchedule from_stride(long long s) {
        OutputSchedule o;
        o.stride = s > 0 ? s : 1;
        return o;
    }
    static OutputSchedule from_frame_dt(double dt) {
        OutputSchedule o;
        o.time_based = true;
        o.frame_dt   = dt;
        return o;
    }

    double next_output_time() const {
        return static_cast<double>(frame_index) * frame_dt;
    }
    double tolerance() const { return 1.0e-9 * frame_dt; }

    // Should a snapshot be written for the state at (step, t)? Never true twice
    // for the same step.
    bool due(long long step, double t) const {
        if (step == last_written_step) return false;
        if (!time_based) return step % stride == 0;
        return t + tolerance() >= next_output_time();
    }

    // Record that (step, t) was written and advance past every scheduled time
    // this state already covers.
    void note_written(long long step, double t) {
        last_written_step = step;
        if (!time_based) return;
        const long long crossed =
            static_cast<long long>(std::floor((t + tolerance()) / frame_dt));
        if (crossed + 1 > frame_index) frame_index = crossed + 1;
    }

    bool already_written(long long step) const { return step == last_written_step; }
};

} // namespace chromosphere
