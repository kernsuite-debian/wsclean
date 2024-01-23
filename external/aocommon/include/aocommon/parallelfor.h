#ifndef AOCOMMON_PARALLEL_FOR_H_
#define AOCOMMON_PARALLEL_FOR_H_

#include "barrier.h"

#include <atomic>
#include <condition_variable>
#include <cstring>
#include <mutex>
#include <stdexcept>
#include <thread>
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
 * alive until the ParallelFor is destructed.
 */
template <typename IterType>
class ParallelFor {
 public:
  ParallelFor(size_t n_threads)
      : n_threads_(n_threads),
        barrier_(n_threads, [&]() { has_tasks_ = false; }),
        stop_(false),
        has_tasks_(false) {}

  ~ParallelFor() {
    std::unique_lock<std::mutex> lock(mutex_);
    if (!threads_.empty()) {
      stop_ = true;
      has_tasks_ = true;
      condition_changed_.notify_all();
      lock.unlock();
      for (std::thread& thr : threads_) thr.join();
    }
  }

  /**
   * Iteratively call a function in parallel (with thread id).
   *
   * The provided function is expected to accept two parameters, the loop
   * index and the thread id, e.g.:
   *   void loopFunction(size_t iteration, size_t threadID);
   * It is called (end-start) times unless an exception occurs.
   *
   * This function is very similar to ThreadPool::For(), but does not
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
    if (end == start + 1 || n_threads_ == 1) {
      for (IterType iter = start; iter != end; ++iter) function(iter, 0);
    } else {
      std::unique_lock<std::mutex> lock(mutex_);
      current_ = start;
      end_ = end;
      loop_function_1_parameter_ = {};
      loop_function_2_parameters_ = std::move(function);
      has_tasks_ = true;
      if (threads_.empty()) StartThreads();
      condition_changed_.notify_all();
      lock.unlock();
      Loop(0);
      barrier_.wait();
      CheckForException();
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
    if (end == start + 1 || n_threads_ == 1) {
      for (IterType iter = start; iter != end; ++iter) function(iter);
    } else {
      std::unique_lock<std::mutex> lock(mutex_);
      has_tasks_ = true;
      current_ = start;
      end_ = end;
      loop_function_1_parameter_ = std::move(function);
      loop_function_2_parameters_ = {};
      if (threads_.empty()) StartThreads();
      condition_changed_.notify_all();
      lock.unlock();
      Loop(0);
      barrier_.wait();
      CheckForException();
    }
  }

  size_t NThreads() const { return n_threads_; }

  /**
   * This method is only allowed to be called before Run() is
   * called.
   */
  void SetNThreads(size_t n_threads) {
    if (threads_.empty()) {
      n_threads_ = n_threads;
      barrier_ = Barrier(n_threads, [&]() { has_tasks_ = false; });
    } else {
      throw std::runtime_error("Can not set NThreads after calling Run()");
    }
  }

 private:
  ParallelFor(const ParallelFor&) = delete;

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
   * Keep running loops until the class is destructed.
   */
  void RunLoops(size_t thread) {
    while (!stop_) {
      Loop(thread);
      barrier_.wait();
      WaitForTasks();
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

  void WaitForTasks() {
    std::unique_lock<std::mutex> lock(mutex_);
    while (!has_tasks_) condition_changed_.wait(lock);
  }

  void StartThreads() {
    if (n_threads_ > 1) {
      threads_.reserve(n_threads_ - 1);
      for (size_t t = 1; t != n_threads_; ++t)
        threads_.emplace_back(&ParallelFor::RunLoops, this, t);
    }
  }

  IterType current_;
  IterType end_;
  std::mutex mutex_;
  size_t n_threads_;
  Barrier barrier_;
  std::atomic<bool> stop_;
  bool has_tasks_;
  std::condition_variable condition_changed_;
  std::vector<std::thread> threads_;
  std::function<void(size_t, size_t)> loop_function_2_parameters_;
  std::function<void(size_t)> loop_function_1_parameter_;
  std::exception_ptr most_recent_exception_;
};
}  // namespace aocommon

#endif
