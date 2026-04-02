#ifndef AOCOMMON_DYNAMIC_FOR_H_
#define AOCOMMON_DYNAMIC_FOR_H_

#include "recursivefor.h"
#include "threadpool.h"

#include <cstring>
#include <mutex>
#include <stdexcept>
#include <vector>

namespace aocommon {

/**
 * Run a loop in parallel. In this class, loops are "load balanced", i.e.,
 * if the iteration of one thread takes more time, other threads will perform
 * more iterations (this is sometimes called "dynamic" timing).
 *
 * The downside of load balancing is that every iteration involves a virtual
 * call, which may be relatively expensive when iterations themselves involve
 * very little work. In those cases, either use the @ref StaticFor class
 * or increase the load per iteration (thereby decreasing the nr of iterations).
 *
 * Once started, the threads are reused in second calls of Run(), and are kept
 * alive until the DynamicFor is destructed.
 */
template <typename IterType>
class DynamicFor {
 public:
  DynamicFor() = default;

  /**
   * Iteratively call a function in parallel (with thread id).
   *
   * The provided function is expected to accept two parameters, the loop
   * index and the thread id, e.g.:
   *   void loopFunction(size_t iteration, size_t threadID);
   * It is called (end-start) times unless an exception occurs.
   *
   * This function is very similar to RecursiveFor::For(), but does not
   * support recursion. For non-recursive loop, this function will be
   * faster. The function will block until all iterations have been
   * performed.
   *
   * If exceptions occur, the latest occurring exception will be
   * rethrown in the calling thread. In such cases, not all iterations
   * might be performed.
   */
  void Run(IterType start, IterType end,
           std::function<void(IterType, size_t)> function) {
    if (RecursiveFor* recursive_for = RecursiveFor::GetInstance();
        recursive_for) {
      // There's a recursive for alive, meaning this call is a nested loop.
      // Dispatch this Run to RecursiveFor:
      recursive_for->Run(start, end, function);
    } else {
      const size_t n_threads = ThreadPool::GetInstance().NThreads();
      if (end == start + 1 || n_threads == 1) {
        for (IterType iter = start; iter != end; ++iter) function(iter, 0);
      } else {
        current_ = start;
        end_ = end;
        loop_function_1_parameter_ = {};
        loop_function_2_parameters_ = std::move(function);
        ThreadPool::GetInstance().ExecuteInParallel(
            [&](size_t thread_index) { Loop(thread_index); });
        CheckForException();
      }
    }
  }

  /**
   * Iteratively call a function in parallel (without thread id).
   *
   * The provided function is expected to take only the loop index
   * as parameter. If the thread ID is required, use the other overload.
   * This function behaves otherwise equal to the other overload.
   *
   * For further info including exception behaviour, see the other overload:
   * @ref Run(IterType, IterType, std::function<void(IterType, size_t)>)
   *
   */
  void Run(IterType start, IterType end,
           std::function<void(IterType)> function) {
    if (RecursiveFor* recursive_for = RecursiveFor::GetInstance();
        recursive_for) {
      // There's a recursive for alive, meaning this call is a nested loop.
      // Dispatch this Run to RecursiveFor:
      recursive_for->Run(start, end, function);
    } else {
      const size_t n_threads = ThreadPool::GetInstance().NThreads();
      if (end == start + 1 || n_threads == 1) {
        for (IterType iter = start; iter != end; ++iter) function(iter);
      } else {
        current_ = start;
        end_ = end;
        loop_function_1_parameter_ = std::move(function);
        loop_function_2_parameters_ = {};
        ThreadPool::GetInstance().ExecuteInParallel(
            [&](size_t thread_index) { Loop(thread_index); });
        CheckForException();
      }
    }
  }

 private:
  DynamicFor(const DynamicFor&) = delete;

  /**
   * Throw if an exception occurred and reset exception state.
   */
  void CheckForException() {
    if (most_recent_exception_) {
      std::exception_ptr to_throw = std::move(most_recent_exception_);
      most_recent_exception_ = std::exception_ptr();
      std::rethrow_exception(to_throw);
    }
  }

  /**
   * Keep doing iterations until there are no more iterations necessary.
   */
  void Loop(size_t thread) {
    try {
      IterType iter;
      while (Next(iter)) {
        if (loop_function_2_parameters_) {
          loop_function_2_parameters_(iter, thread);
        } else {
          loop_function_1_parameter_(iter);
        }
      }
    } catch (std::exception&) {
      std::lock_guard<std::mutex> lock(mutex_);
      most_recent_exception_ = std::current_exception();
    }
  }

  /**
   * Obtain the next iteration number. Method is safe to call from multiple
   * threads.
   * @returns false if there are no more iterations necessary.
   */
  bool Next(IterType& iter) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (current_ == end_) {
      return false;
    } else {
      iter = current_;
      ++current_;
      return true;
    }
  }

  IterType current_;
  IterType end_;
  std::mutex mutex_;
  std::function<void(size_t, size_t)> loop_function_2_parameters_;
  std::function<void(size_t)> loop_function_1_parameter_;
  std::exception_ptr most_recent_exception_;
};

}  // namespace aocommon

#endif
