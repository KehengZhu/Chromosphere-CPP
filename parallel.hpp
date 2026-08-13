#pragma once
/*!
 * @file parallel.hpp
 * @brief Thin OpenMP abstraction: thread queries, a deterministic per-thread
 *        exception collector, and the cell-parallel loop the solvers use.
 * @ingroup runtime
 *
 * The whole file degrades to plain serial loops when the build is configured
 * without `CHROMO_ENABLE_OPENMP`, so solver code never has to branch on it.
 * parallel_for_cells() makes failure deterministic: a throwing cell is recorded
 * per thread and the exception from the LOWEST physical cell index is rethrown
 * after the region, so a parallel run reports the same error as a serial one.
 */

#include <array>
#include <cstddef>
#include <exception>
#include <limits>
#include <stdexcept>

#ifdef CHROMO_USE_OPENMP
#include <omp.h>
#endif

namespace chromosphere {

/** @addtogroup runtime
 *  @{
 */

/// Fixed capacity of every per-thread slot array in the solver (ParallelFailure's
/// index/exception slots, the profiling counter slots, the EOS operation counters,
/// the conduction residual and step slots). Those arrays are statically sized, so
/// a thread index must never reach this value; require_supported_parallel_runtime()
/// enforces the bound instead of letting a large OMP_NUM_THREADS overrun them.
constexpr int kMaximumParallelThreads = 256;

/// Upper bound on the team size used for the implicit conduction stage, which
/// takes min(this, parallel_max_threads()) threads. The stage is a barrier-heavy
/// nonlinear Newton sequence around a serial Thomas solve, so on production-shaped
/// meshes it stops scaling and then regresses: measured per 0.1 physical-second
/// probe at 2638 cells, 2 threads ~0.092 s and 4 threads ~0.087 s, versus ~0.132 s
/// at 8 and ~0.204 s at 12 (docs/openmp_parallelization_recap.md). Every other
/// stage still uses the full requested team.
constexpr int kMaximumConductionThreads = 4;

/// True when the build was configured with OpenMP (`CHROMO_ENABLE_OPENMP`, which
/// defines CHROMO_USE_OPENMP); false for the serial build. Reported in the runtime
/// banner so a run records which build produced it.
inline bool openmp_compiled() noexcept {
#ifdef CHROMO_USE_OPENMP
    return true;
#else
    return false;
#endif
}

/// Index of the calling thread within the current parallel region, and 0 whenever
/// the caller is outside one (or in a serial build). This is the slot selector for
/// every per-thread array, so serial code and worker 0 share one slot.
inline int parallel_thread_index() noexcept {
#ifdef CHROMO_USE_OPENMP
    return omp_in_parallel() ? omp_get_thread_num() : 0;
#else
    return 0;
#endif
}

/// Largest team the OpenMP runtime would give the next parallel region (typically
/// `OMP_NUM_THREADS`); always 1 in a serial build. This is the count checked
/// against kMaximumParallelThreads and used to size the conduction team.
inline int parallel_max_threads() noexcept {
#ifdef CHROMO_USE_OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

/// Number of threads in the parallel region the caller is actually executing in,
/// and 1 outside any region. Unlike parallel_max_threads() this reports the team
/// that was granted — the conduction stage records it for its runtime banner.
inline int parallel_team_size() noexcept {
#ifdef CHROMO_USE_OPENMP
    return omp_in_parallel() ? omp_get_num_threads() : 1;
#else
    return 1;
#endif
}

/// Fail closed if the runtime would hand out more threads than the fixed
/// per-thread slot arrays can hold. Call it from serial code before entering a
/// region that writes those slots; it throws std::runtime_error rather than
/// letting a thread index run past kMaximumParallelThreads.
inline void require_supported_parallel_runtime() {
    const int count = parallel_max_threads();
    if (count < 1 || count > kMaximumParallelThreads)
        throw std::runtime_error("OpenMP thread count exceeds fixed solver slot capacity");
}

/// Deterministic exception collector. Each worker writes only its own slot; after
/// the region, the caller rethrows the exception from the lowest physical cell.
class ParallelFailure {
public:
    /// Start empty: no slot holds an exception.
    ParallelFailure() { reset(); }

    /// Drop every recorded exception and reset all cell indices to "none", so the
    /// collector can be reused for another region.
    void reset() noexcept {
        indices_.fill(std::numeric_limits<std::size_t>::max());
        for (auto& error : errors_) error = std::exception_ptr{};
    }

    /// Record that physical cell @p cell threw @p error. The calling worker writes
    /// only its own slot — no lock, no atomic — and keeps just the LOWEST cell
    /// index it has seen, which is what makes rethrow_lowest() deterministic. A
    /// thread index outside the slot capacity is dropped rather than written out
    /// of bounds. Safe to call from a `catch (...)` block: it cannot itself throw.
    void capture(std::size_t cell, std::exception_ptr error) noexcept {
        const int tid = parallel_thread_index();
        if (tid < 0 || tid >= kMaximumParallelThreads) return;
        const std::size_t slot = static_cast<std::size_t>(tid);
        if (cell < indices_[slot]) {
            indices_[slot] = cell;
            errors_[slot] = error;
        }
    }

    /// True if any slot holds an exception. Used inside a region (from an `omp
    /// single` block) to decide whether the remaining work should be abandoned.
    bool failed() const noexcept {
        for (const auto& error : errors_)
            if (error) return true;
        return false;
    }

    /// Rethrow the exception belonging to the lowest physical cell index recorded
    /// by any thread; do nothing if none was recorded. Call this after the region
    /// has joined. Because the winner is chosen by cell index and not by which
    /// thread happened to fail first, the reported error is independent of the
    /// thread count and matches what a serial run would have thrown.
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
    /// Lowest failing physical cell index seen by each thread; size_t max means
    /// "that thread has not failed". Cache-line aligned to avoid false sharing.
    alignas(64) std::array<std::size_t, kMaximumParallelThreads> indices_;
    /// The exception each thread captured for its indices_ entry, one slot per
    /// thread, empty when that thread has not failed.
    alignas(64) std::array<std::exception_ptr, kMaximumParallelThreads> errors_;
};

/// Apply @p function to every cell index in [0, @p count), in parallel when the
/// build and runtime allow it. A single-thread runtime and the serial build both
/// take a plain loop, so the callable sees identical arguments either way; the
/// caller must therefore keep the per-index work independent.
///
/// Exceptions are deterministic rather than fatal: OpenMP forbids escaping a
/// parallel region, so a throwing index is captured per thread and, after the
/// region joins, the exception from the lowest cell index is rethrown. A parallel
/// run reports exactly the error a serial run would.
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

/** @} */

} // namespace chromosphere
