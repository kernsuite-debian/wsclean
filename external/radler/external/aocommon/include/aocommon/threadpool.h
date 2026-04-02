#ifndef AOCOMMON_THREAD_POOL_H_
#define AOCOMMON_THREAD_POOL_H_

#include <cassert>
#include <condition_variable>
#include <mutex>
#include <thread>

#include "barrier.h"

namespace aocommon {

/**
 * This class holds a running list of threads that can be used to
 * run multiple functions in parallel without having to stop and
 * recreate the threads in between.
 *
 * There's a static instance of the pool that can be used
 * to avoid passing the class between many functions.
 */
class ThreadPool {
 public:
  ThreadPool() = default;
  ~ThreadPool() { StopThreads(); }

  /**
   * Returns a static instance of this class. It is by default initialized
   * with only one thread, but may be changed globally using SetNThreads().
   */
  static ThreadPool& GetInstance() {
    static ThreadPool instance_;
    return instance_;
  }

  size_t NThreads() const { return threads_.size() + 1; }

  /**
   * This makes the class ready to execute with n_threads parallel
   * threads. It will spawn n-1 new threads,
   * because @ref ExecuteInParallel() also uses the calling thread.
   * The running threads can be terminated by setting n_threads to 1.
   * This function should not be called from parallel executed
   * functions themselves.
   */
  void SetNThreads(size_t n_threads) {
    assert(execute_function_ == nullptr);
    assert(n_threads != 0);
    StopThreads();
    StartThreads(n_threads);
  }

  /**
   * Stops all threads. Afterwards, threads can be restarted with a call to
   * SetNThreads().
   */
  void Stop() { SetNThreads(1); }

  /**
   * Run the specified function in parallel. The function will be called
   * n_threads times. The called function should not throw exceptions.
   * @param execute_function Function with the thread_id as a parameter.
   */
  void ExecuteInParallel(std::function<void(size_t)> execute_function) {
    assert(execute_function_ == nullptr);
    if (threads_.empty()) {
      execute_function(0);
    } else {
      std::unique_lock lock(mutex_);
      execute_function_ = execute_function;
      condition_.notify_all();
      lock.unlock();

      execute_function(0);

      lock.lock();
      BarrierWait(lock);
      execute_function_ = nullptr;
    }
  }

  /**
   * Similar to @ref ExecuteInParallel(), but non-blocking. This function
   * calls the execute_function "n_threads-1" times. The calling thread
   * is not used to call the function, thereby allowing the function to
   * return immediately. In this case, the thread index provided to
   * execute_function is always >= 1. A call to ExecuteInParallel() must
   * be followed by a call to @ref FinishParallelExecution().
   */
  void StartParallelExecution(std::function<void(size_t)> execute_function) {
    assert(execute_function_ == nullptr);
    if (!threads_.empty()) {
      std::lock_guard lock(mutex_);
      execute_function_ = execute_function;
      condition_.notify_all();
    }
  }

  /**
   * Finish a call to StartParallelExecution() and wait for all executing
   * functions to finish. This call blocks until that condition is satisfied.
   * This function should only be called once after a call to
   * StartParallelExecution(). It should not be called to finish a call to
   * ExecuteInParallel().
   */
  void FinishParallelExecution() {
    std::unique_lock lock(mutex_);
    BarrierWait(lock);
    execute_function_ = nullptr;
  }

 private:
  void StartThreads(size_t n_threads) {
    assert(n_threads != 0);
    n_executing_ = n_threads;
    for (size_t i = 1; i != n_threads; ++i) {
      threads_.emplace_back([i, this]() { Run(i); });
    }
  }

  void StopThreads() {
    std::unique_lock lock(mutex_);
    stop_ = true;
    condition_.notify_all();
    lock.unlock();

    for (std::thread& t : threads_) t.join();
    threads_.clear();
    stop_ = false;
  }

  void Run(size_t thread_index) {
    std::unique_lock lock(mutex_);
    while (true) {
      while (execute_function_ == nullptr && !stop_) condition_.wait(lock);
      if (stop_) break;

      lock.unlock();
      execute_function_(thread_index);
      lock.lock();
      BarrierWait(lock);
    }
  }

  void BarrierWait(std::unique_lock<std::mutex>& lock) {
    --n_executing_;
    if (n_executing_ == 0) {
      ++barrier_cycle_;
      n_executing_ = NThreads();
      // All threads have finished executing their function and have arrived at
      // the barrier: make the threads wait until the next call to
      // ExecuteInParallel().
      execute_function_ = nullptr;
      condition_.notify_all();
    } else {
      const size_t cycle = barrier_cycle_;
      while (cycle == barrier_cycle_) condition_.wait(lock);
    }
  }

  std::vector<std::thread> threads_;
  bool stop_ = false;
  size_t n_executing_ = 0;
  size_t barrier_cycle_ = 0;
  /// The function that is executed in parallel. When it is empty (nullptr),
  /// there's no task for the threads.
  std::function<void(size_t)> execute_function_;
  std::mutex mutex_;
  std::condition_variable condition_;
};

}  // namespace aocommon

#endif
