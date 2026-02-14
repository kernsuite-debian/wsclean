#ifndef AOCOMMON_TASK_QUEUE_H_
#define AOCOMMON_TASK_QUEUE_H_

#include <cassert>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <utility>

namespace aocommon {

/**
 * Defines a thread-safe FIFO queue, aimed at distributing tasks to
 * parallel worker processes.
 * Similarly to aocommon::Lane, it supports communicating that all tasks are
 * done / no task will enter the queue anymore.
 * The main difference with aocommon::Lane is that aocommon::Lane has a fixed
 * size buffer while a TaskQueue grows dynamically.
 *
 * Typical usage:
 * 1. Create the queue.
 * 2. Start <n_threads> worker threads, which do:
 *    while (queue.Pop(task)) { <execute 'task'> }
 * 3. Use queue.Emplace() to add a first group of tasks to the queue.
 * 4. Call queue.Wait(<n_threads>) to wait until all tasks are finished.
 * 5. Repeat 3. and 4. for further groups of tasks.
 * 6. Call queue.Finish(). The worker threads will exit.
 * 7. Call thread.join() for the worker threads.
 */
template <typename T>
class TaskQueue {
 public:
  /**
   * Creates a queue that stores an unlimited number of tasks.
   * Push() always returns immediately.
   */
  TaskQueue() noexcept {}

  /**
   * Creates a queue that limits the number of tasks in the queue.
   * Push() will wait until there is room in the queue.
   * @param limit Maximum number of tasks in the queue.
   */
  explicit TaskQueue(size_t limit) noexcept : limit_{limit} {
    assert(limit != 0);
  }

  /**
   * Creates a task at the end of the queue, using perfect forwarding.
   * @param args Arguments for constructing the task object.
   */
  template <class... Args>
  void Emplace(Args&&... args) {
    std::unique_lock<std::mutex> lock(mutex_);
    assert(!finish_);
    if (limit_ > 0) {
      while (tasks_.size() >= limit_) remove_notifier_.wait(lock);
    }
    tasks_.emplace(std::forward<Args>(args)...);
    insert_notifier_.notify_one();
  }

  /**
   * Waits until all work is done, which means the queue is empty and
   * all threads are idle.
   * @param n_threads The number of threads pulling work from the queue.
   */
  void WaitForIdle(size_t n_threads) {
    std::unique_lock<std::mutex> lock(mutex_);
    while (!tasks_.empty() || wait_count_ < n_threads) {
      wait_notifier_.wait(lock);
    }
  }

  /**
   * Tells the queue that no more tasks will enter the queue.
   * After returning the last task, Pop() will only return false.
   */
  void Finish() {
    std::lock_guard<std::mutex> lock(mutex_);
    finish_ = true;
    insert_notifier_.notify_all();
  }

  /**
   * Gets a task from the front of the queue.
   * If no task is available, waits until there is one.
   * @param destination [out] Reference for storing the task.
   * @return true: A task was succesfully returned.
   *        false: Processing has stopped since Finish() was called.
   * @see Finish()
   */
  bool Pop(T& destination) {
    std::unique_lock<std::mutex> lock(mutex_);
    while (!finish_ && tasks_.empty()) {
      ++wait_count_;
      wait_notifier_.notify_one();
      insert_notifier_.wait(lock);
      --wait_count_;
    }
    if (finish_) {
      return false;
    } else {
      assert(!tasks_.empty());
      destination = std::move(tasks_.front());
      tasks_.pop();
      if (limit_ > 0) remove_notifier_.notify_one();
      return true;
    }
  }

 private:
  std::mutex mutex_;
  std::condition_variable insert_notifier_;
  std::condition_variable remove_notifier_;
  std::condition_variable wait_notifier_;
  std::queue<T> tasks_;
  bool finish_ = false;
  size_t wait_count_ = 0;
  size_t limit_ = 0;  // 0 means there is no queue limit.
};

}  // namespace aocommon

#endif