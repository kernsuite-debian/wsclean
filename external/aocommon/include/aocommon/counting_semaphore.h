#ifndef AOCOMMON_COUNTING_SEMAPHORE_H_
#define AOCOMMON_COUNTING_SEMAPHORE_H_

#include <cassert>
#include <condition_variable>
#include <limits>
#include <mutex>

namespace aocommon {

/**
 * This replicates std::counting_semaphore from C++20.
 */
class CountingSemaphore {
 public:
  CountingSemaphore(std::size_t initial_value) : value_(initial_value) {}

  void Acquire() {
    std::unique_lock lock(mutex_);
    while (value_ == 0) {
      condition_.wait(lock);
    }
    --value_;
  }

  void Release() {
    std::scoped_lock lock(mutex_);
    ++value_;
    condition_.notify_one();
  }

 private:
  std::mutex mutex_;
  std::condition_variable condition_;
  std::size_t value_;
};

class ScopedCountingSemaphoreLock {
 public:
  ScopedCountingSemaphoreLock(CountingSemaphore& semaphore)
      : semaphore_(semaphore) {
    semaphore_.Acquire();
  }

  ~ScopedCountingSemaphoreLock() { semaphore_.Release(); }

 private:
  CountingSemaphore& semaphore_;
};

}  // namespace aocommon

#endif  // AOCOMMON_COUNTING_SEMAPHORE_H_
