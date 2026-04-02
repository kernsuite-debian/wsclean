#ifndef AOCOMMON_STATIC_FOR_H_
#define AOCOMMON_STATIC_FOR_H_

#include "threadpool.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <condition_variable>
#include <cstring>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

#include <sched.h>

namespace aocommon {

/**
 * The StaticFor class implements a parallel for loop that
 * is statically distributed over all threads. It is suitable
 * for large loops that approximately take an equal amount of
 * time per iteration. If one thread finishes its chunk
 * earlier, it won't be rescheduled to do more.
 *
 * The advantage of this is that it doesn't require communication
 * between iterations, and thus the loop is faster, as long as
 * iterations take similar time.
 *
 * StaticFor makes use of the static ThreadPool instance OR if
 * a RecursiveFor instance is alive, it will nest itself inside
 * the RecursiveFor. To use multiple threads, call
 * ThreadPool::GetInstance().SetNThreads() beforehand.
 *
 * An example (unnested) for-loop to count to 1000 in parallel:
 * StaticFor<size_t> loop;
 * loop.Run(0, 1000, [&](size_t start, size_t end) {
 *   for(size_t i=start; i!=end; ++i) {
 *     std::cout << i << '\n';
 *   }
 * }
 *
 * It is also possible to acquire the thread index, by providing
 * a function with 3 parameters:
 * StaticFor<size_t> loop;
 * loop.Run(0, 1000, [&](size_t start, size_t end, size_t thread) {
 *   for(size_t i=start; i!=end; ++i) {
 *     std::cout << i << " from thread " << thread << '\n';
 *   }
 * }
 *
 * Exceptions are not handled, which implies that an uncaught exception
 * thrown while iterating might get thrown in the main thread or might cause
 * immediate termination when it occurred in a separate thread.
 */
template <typename Iter>
class StaticFor {
 public:
  StaticFor() = default;

  /**
   * Iteratively call a function in parallel. Requires start <= end.
   *
   * The provided function is expected to accept two parameters, the
   * start and end indices of this thread, e.g.:
   *   void loopFunction(size_t chunk_start, size_t chunk_end);
   */
  void Run(Iter start, Iter end, std::function<void(Iter, Iter)> function) {
    loop_function_without_id_ = std::move(function);
    run(start, end);
    loop_function_without_id_ = nullptr;
  }

  /**
   * Iteratively call a function in parallel with thread id. Requires start <=
   * end.
   *
   * The provided function is expected to accept three parameters, the start
   * and end indices of this thread and the thread index, e.g.:
   *   void loopFunction(size_t chunk_start, size_t chunk_end, size_t
   * thread_index);
   */
  void Run(Iter start, Iter end,
           std::function<void(Iter, Iter, size_t)> function) {
    loop_function_with_id_ = std::move(function);
    run(start, end);
    loop_function_with_id_ = nullptr;
  }

  /**
   * Like @ref Run(), but with limited number of threads. This
   * may be useful when: i) each thread takes memory and the total
   * memory needs to be constrained -- in combination with
   * a nested for (using the @ref RecursiveFor class), it may
   * be possible to reduce the memory without affecting performance;
   * ii) if the run is part of a nested run, and synchronization
   * of all threads causes too much overhead.
   *
   * If @p max_threads is set to 1, the provided @p function is
   * directly called without any thread synchronization. This also
   * implies that a nested run does not require a RecursiveFor to
   * exist. This allows conditional nested parallelization without
   * having to duplicate the function for the case that doing a
   * nested run is not efficient (e.g. because the outer loop
   * dimension is large enough to occupy all threads).
   */
  void ConstrainedRun(Iter start, Iter end, size_t max_threads,
                      std::function<void(Iter, Iter)> function) {
    loop_function_without_id_ = std::move(function);
    run(start, end, max_threads);
    loop_function_without_id_ = nullptr;
  }

  /**
   * Same as the other ConstrainedRun() overload, but for
   * an iterating function that includes a thread index. Unlike
   * the corresponding Run() call, this function makes sure that the
   * @c thread_index passed to the iteration function is always
   * lower than @p max_threads. This comes at a slight overhead of
   * one extra std::function indirection.
   */
  void ConstrainedRun(Iter start, Iter end, size_t max_threads,
                      std::function<void(Iter, Iter, size_t)> function) {
    std::atomic<size_t> thread_counter = 0;
    loop_function_without_id_ = [&](Iter start, Iter end) {
      const size_t thread_index = thread_counter.fetch_add(1);
      function(start, end, thread_index);
    };
    run(start, end, max_threads);
    loop_function_without_id_ = nullptr;
  }

 private:
  StaticFor(const StaticFor&) = delete;

  void run(Iter start, Iter end,
           size_t max_chunks = std::numeric_limits<size_t>::max());

  void callFunction(Iter start, Iter end, size_t threadId) const {
    if (loop_function_without_id_)
      loop_function_without_id_(start, end);
    else
      loop_function_with_id_(start, end, threadId);
  }

  void loop(size_t thread_index) {
    if (thread_index < n_chunks_) {
      Iter chunk_start =
          iter_start_ + (iter_end_ - iter_start_) * thread_index / n_chunks_;
      Iter chunk_end = iter_start_ + (iter_end_ - iter_start_) *
                                         (thread_index + 1) / n_chunks_;

      callFunction(chunk_start, chunk_end, thread_index);
    }
  }

  /// The number of chunks (parts) that the full range is divided into.
  size_t n_chunks_;
  Iter iter_start_;
  Iter iter_end_;
  std::function<void(Iter, Iter)> loop_function_without_id_;
  std::function<void(Iter, Iter, size_t)> loop_function_with_id_;
};

}  // namespace aocommon

#include "recursivefor.h"

namespace aocommon {

template <typename Iter>
inline void StaticFor<Iter>::run(Iter start, Iter end, size_t max_chunks) {
  assert(start <= end);
  const size_t n_threads = ThreadPool::GetInstance().NThreads();
  n_chunks_ = std::min({n_threads, end - start, max_chunks});
  if (n_chunks_ <= 1) {
    callFunction(start, end, 0);
  } else {
    iter_start_ = start;
    iter_end_ = end;
    if (RecursiveFor* recursive_for = RecursiveFor::GetInstance();
        recursive_for) {
      recursive_for->Run(0, n_chunks_,
                         [&](size_t thread_index) { loop(thread_index); });
    } else {
      ThreadPool::GetInstance().ExecuteInParallel(
          [&](size_t thread_index) { loop(thread_index); });
    }
  }
}

template <typename Iter>
inline void RunStaticFor(Iter start, Iter end,
                         std::function<void(Iter, Iter)> function) {
  StaticFor<Iter>().Run(start, end, std::move(function));
}

template <typename Iter>
inline void RunStaticFor(Iter start, Iter end,
                         std::function<void(Iter, Iter, size_t)> function) {
  StaticFor<Iter>().Run(start, end, std::move(function));
}

template <typename Iter>
inline void RunConstrainedStaticFor(Iter start, Iter end, size_t max_threads,
                                    std::function<void(Iter, Iter)> function) {
  StaticFor<Iter>().ConstrainedRun(start, end, max_threads,
                                   std::move(function));
}

template <typename Iter>
inline void RunConstrainedStaticFor(
    Iter start, Iter end, size_t max_threads,
    std::function<void(Iter, Iter, size_t)> function) {
  StaticFor<Iter>().ConstrainedRun(start, end, max_threads,
                                   std::move(function));
}

}  // namespace aocommon

#endif
