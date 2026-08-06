#pragma once

#include <array>
#include <cstddef>
#include <exception>
#include <limits>
#include <stdexcept>

#ifdef CHROMO_USE_OPENMP
#include <omp.h>
#endif

namespace chromosphere {

constexpr int kMaximumParallelThreads = 256;
constexpr int kMaximumConductionThreads = 4;

inline bool openmp_compiled() noexcept {
#ifdef CHROMO_USE_OPENMP
    return true;
#else
    return false;
#endif
}

inline int parallel_thread_index() noexcept {
#ifdef CHROMO_USE_OPENMP
    return omp_in_parallel() ? omp_get_thread_num() : 0;
#else
    return 0;
#endif
}

inline int parallel_max_threads() noexcept {
#ifdef CHROMO_USE_OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

inline int parallel_team_size() noexcept {
#ifdef CHROMO_USE_OPENMP
    return omp_in_parallel() ? omp_get_num_threads() : 1;
#else
    return 1;
#endif
}

inline void require_supported_parallel_runtime() {
    const int count = parallel_max_threads();
    if (count < 1 || count > kMaximumParallelThreads)
        throw std::runtime_error("OpenMP thread count exceeds fixed solver slot capacity");
}

/// Deterministic exception collector. Each worker writes only its own slot; after
/// the region, the caller rethrows the exception from the lowest physical cell.
class ParallelFailure {
public:
    ParallelFailure() { reset(); }

    void reset() noexcept {
        indices_.fill(std::numeric_limits<std::size_t>::max());
        for (auto& error : errors_) error = std::exception_ptr{};
    }

    void capture(std::size_t cell, std::exception_ptr error) noexcept {
        const int tid = parallel_thread_index();
        if (tid < 0 || tid >= kMaximumParallelThreads) return;
        const std::size_t slot = static_cast<std::size_t>(tid);
        if (cell < indices_[slot]) {
            indices_[slot] = cell;
            errors_[slot] = error;
        }
    }

    bool failed() const noexcept {
        for (const auto& error : errors_)
            if (error) return true;
        return false;
    }

    void rethrow_lowest() const {
        std::size_t best = std::numeric_limits<std::size_t>::max();
        std::exception_ptr chosen;
        for (int tid = 0; tid < kMaximumParallelThreads; ++tid) {
            const std::size_t slot = static_cast<std::size_t>(tid);
            if (errors_[slot] && indices_[slot] < best) {
                best = indices_[slot];
                chosen = errors_[slot];
            }
        }
        if (chosen) std::rethrow_exception(chosen);
    }

private:
    alignas(64) std::array<std::size_t, kMaximumParallelThreads> indices_;
    alignas(64) std::array<std::exception_ptr, kMaximumParallelThreads> errors_;
};

template <class Function>
void parallel_for_cells(std::size_t count, Function function) {
#ifdef CHROMO_USE_OPENMP
    require_supported_parallel_runtime();
    if (parallel_max_threads() == 1) {
        for (std::size_t i = 0; i < count; ++i) function(i);
        return;
    }
    ParallelFailure failure;
#pragma omp parallel for schedule(static)
    for (long long raw = 0; raw < static_cast<long long>(count); ++raw) {
        const std::size_t i = static_cast<std::size_t>(raw);
        try {
            function(i);
        } catch (...) {
            failure.capture(i, std::current_exception());
        }
    }
    failure.rethrow_lowest();
#else
    for (std::size_t i = 0; i < count; ++i) function(i);
#endif
}

} // namespace chromosphere
