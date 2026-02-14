#ifndef AOCOMMON_RECURSIVE_FOR_H_
#define AOCOMMON_RECURSIVE_FOR_H_

#include "system.h"
#include "threadpool.h"

#include <sched.h>

#include <condition_variable>
#include <functional>
#include <map>
#include <mutex>
#include <thread>
#include <vector>

namespace aocommon {

/**
 * Parallellizes a recursive loop. The iterations are executed in a
 * first in last out progression, meaning that the iterations corresponding with
 * the latest call to @ref RecursiveFor::Run() will have highest priority.
 * Nevertheless, at least one thread (the calling thread) will stay assigned to
 * perform iterations. Iterations are scheduled dynamically.
 *
 * The global thread pool remains occupied as long as the RecursiveFor is alive.
 * Other uses of the thread pool (e.g. with @ref DynamicFor or @ref StaticFor)
 * are therefore not allowed during the lifetime of the @ref RecursiveFor.
 *
 * A reason to use this class over for example only a @ref DynamicFor for the
 * outer loop, is if the outer loop may not have enough iterations to keep all
 * threads busy. An example is to loop over number of directions and number of
 * sources in Dp3: the number of directions may be one, which would cause serial
 * processing. Similarly, the number of sources per direction may be low, and
 * only parallelizing over sources may thus also not be efficient.
 */
class RecursiveFor {
 public:
  RecursiveFor() {
    // There should never be more than one instance
    assert(instance_ == nullptr);
    instance_ = this;
    ThreadPool::GetInstance().StartParallelExecution(
        [&](size_t thread_index) { RunIterations(thread_index); });
  }

  RecursiveFor(const RecursiveFor&) = delete;
  RecursiveFor& operator=(const RecursiveFor&) = delete;

  ~RecursiveFor() {
    std::unique_lock<std::mutex> lock(mutex_);
    stopped_ = true;
    new_task_condition_.notify_all();
    lock.unlock();
    ThreadPool::GetInstance().FinishParallelExecution();
    instance_ = nullptr;
  }

  /**
   * If an instance is currently alive, return that instance. Otherwise,
   * return nullptr.
   */
  static RecursiveFor* GetInstance() { return instance_; }

  /**
   * Iteratively call a function in parallel.
   *
   * The provided function is expected to accept one or two size_t parameters,
   * the loop index and optionally the thread id, e.g.:
   *   void loopFunction(size_t iteration, size_t thread_id);
   * It is called (end-start) times. The function blocks until all iterations
   * have been done.
   *
   * Exceptions are caught: if an exception occurs in an iteration, the
   * iteration is skipped and the for loop continues. Once all iterations are
   * done with one or more exceptions, the @ref Run() call will rethrow the last
   * occured exception into the calling thread.
   */
  template <typename Function>
  void Run(size_t start, size_t end, Function func) {
    std::unique_lock<std::mutex> lock(mutex_);
    const size_t this_priority = priority_;
    ++priority_;

    Task& task =
        tasks_.emplace(this_priority, Task(func, start, end)).first->second;
    new_task_condition_.notify_all();
    FinishSpecificPriority(this_priority, task, lock);
  }

  /**
   * Like @ref Run(), but reuses the global instance if it is alive.
   * If no global instance is alive, it creates a temporary
   * RecursiveFor.
   */
  template <typename Function>
  static void NestedRun(size_t start, size_t end, Function func) {
    if (instance_)
      instance_->Run(start, end, func);
    else {
      RecursiveFor().Run(start, end, func);
    }
  }

  /**
   * Performs a constrained statically scheduled loop for a function without
   * thread id. It calls @ref StaticFor::ConstrainedRun();
   */
  void ConstrainedRun(size_t start, size_t end, size_t max_threads,
                      std::function<void(size_t, size_t)> function);

  /**
   * Performs a constrained statically scheduled loop for a function with
   * thread id. It calls @ref StaticFor::ConstrainedRun();
   */
  void ConstrainedRun(size_t start, size_t end, size_t max_threads,
                      std::function<void(size_t, size_t, size_t)> function);

 private:
  struct Task {
    Task(const std::function<void(size_t, size_t)>& with, size_t start,
         size_t _end)
        : function_with_thread_index(with),
          iterator(start),
          end(_end),
          n_unfinished(_end - start) {}
    Task(const std::function<void(size_t)>& without, size_t start, size_t _end)
        : function_without_thread_index(without),
          iterator(start),
          end(_end),
          n_unfinished(_end - start) {}

    void Call(size_t iterator, size_t thread_index) const {
      if (function_with_thread_index)
        function_with_thread_index(iterator, thread_index);
      else
        function_without_thread_index(iterator);
    }

    std::function<void(size_t, size_t)> function_with_thread_index;
    std::function<void(size_t)> function_without_thread_index;
    /// The iterator value for the next thread (updated at the start of an
    /// iteration)
    size_t iterator;
    size_t end;
    /// The number of iterator values left (updated after finishing an
    /// iteration)
    size_t n_unfinished;
    std::exception_ptr most_recent_exception_;
  };

  void RunIterations(size_t thread_index) {
    std::unique_lock<std::mutex> lock(mutex_);
    while (!stopped_) {
      Task* task = HighestPriorityTask();
      while (!stopped_ && !task) {
        new_task_condition_.wait(lock);
        task = HighestPriorityTask();
      }
      if (task) {
        const size_t iterator = task->iterator;
        ++task->iterator;
        lock.unlock();

        try {
          task->Call(iterator, thread_index);
          lock.lock();
        } catch (std::exception& e) {
          lock.lock();
          task->most_recent_exception_ = std::current_exception();
        }

        --task->n_unfinished;
        // Use notify_all(), since multiple FinishSpecificPriority calls
        // can be waiting for finished_task_condition_.
        if (task->n_unfinished == 0) finished_task_condition_.notify_all();
      }
    }
  }

  void FinishSpecificPriority(size_t priority, Task& task,
                              std::unique_lock<std::mutex>& lock) {
    // When we are inside the function, an unfinished Run() call is
    // ongoing. This implies that the RecursiveFor cannot be destructed
    // at that point. Therefore it is unnecessary to check stopped_.
    while (task.iterator != task.end) {
      const size_t iterator = task.iterator;
      ++task.iterator;
      lock.unlock();

      try {
        // A thread index of zero can be used because this index is only used
        // inside the FinishSpecificPriority function, and because each Run
        // call will call FinishSpecificPriority exactly once and with a unique
        // priority. Therefore, no for loop will be given thread index 0 twice.
        constexpr size_t kThreadIndex = 0;
        task.Call(iterator, kThreadIndex);
        lock.lock();
      } catch (std::exception& e) {
        lock.lock();
        task.most_recent_exception_ = std::current_exception();
      }

      --task.n_unfinished;
    }
    // All iterator values have been taken by threads. However, some threads
    // might still be running an iteration, so wait for them to finish
    while (task.n_unfinished != 0) {
      finished_task_condition_.wait(lock);
    }
    std::exception_ptr to_throw = std::move(task.most_recent_exception_);
    tasks_.erase(priority);
    if (to_throw) {
      std::rethrow_exception(to_throw);
    }
  }

  /**
   * This function must be called with a lock on @c mutex_.
   */
  Task* HighestPriorityTask() {
    // tasks_ is ordered using std::greater, thus with highest priority first.
    for (std::pair<const size_t, Task>& task : tasks_) {
      if (task.second.iterator != task.second.end) {
        return &task.second;
      }
    }
    return nullptr;
  }

  bool stopped_ = false;
  size_t priority_ = 0;
  std::map<size_t, Task, std::greater<size_t>> tasks_;
  std::mutex mutex_;
  std::condition_variable new_task_condition_;
  std::condition_variable finished_task_condition_;
  inline static RecursiveFor* instance_ = nullptr;
};

}  // namespace aocommon

#include "staticfor.h"

namespace aocommon {

inline void RecursiveFor::ConstrainedRun(
    size_t start, size_t end, size_t max_threads,
    std::function<void(size_t, size_t)> function) {
  StaticFor<size_t> loop;
  loop.ConstrainedRun(start, end, max_threads, std::move(function));
}

inline void RecursiveFor::ConstrainedRun(
    size_t start, size_t end, size_t max_threads,
    std::function<void(size_t, size_t, size_t)> function) {
  StaticFor<size_t> loop;
  loop.ConstrainedRun(start, end, max_threads, std::move(function));
}

}  // namespace aocommon

#endif
