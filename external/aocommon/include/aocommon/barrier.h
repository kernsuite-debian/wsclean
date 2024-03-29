#ifndef AOCOMMON_BARRIER_H_
#define AOCOMMON_BARRIER_H_

#include <condition_variable>
#include <functional>
#include <mutex>

namespace aocommon {

/**
 * This class is unfortunately necessary because boost::barrier had a
 * bug in completion functions, and std::barrier is still experimental.
 */
class Barrier {
 public:
  /**
   * Construct barrier for n threads without a completion function.
   * @param n Number of threads to wait for
   * @param completionFunction void function that is called when all threads
   * have arrived, just before the threads are released.
   */
  Barrier(size_t n)
      : _n(n), _count(_n), _cycle(0), _completionFunction(&Barrier::nop) {}

  /**
   * Construct barrier for n threads with the given completion function.
   * @param n Number of threads to wait for
   * @param completionFunction void function that is called when all threads
   * have arrived, just before the threads are released.
   */
  Barrier(size_t n, std::function<void()> completionFunction)
      : _n(n), _count(_n), _cycle(0), _completionFunction(completionFunction) {}

  /**
   * @brief Move constructor
   *
   * count could be checked for != _n to assert no threads
   * are waiting, but I leave this the responsibility of the caller.
   */
  Barrier& operator=(Barrier&& rhs) {
    _n = rhs._n;
    _count = _n;
    _cycle = 0;
    _completionFunction = rhs._completionFunction;
    return *this;
  }

  /**
   * Wait until all threads are waiting for the barrier.
   */
  void wait() {
    std::unique_lock<std::mutex> lock(_mutex);
    --_count;
    if (_count == 0) {
      ++_cycle;
      _count = _n;
      _completionFunction();
      _condition.notify_all();
    } else {
      size_t cycle = _cycle;
      while (cycle == _cycle) _condition.wait(lock);
    }
  }

 private:
  static void nop() {}

  std::mutex _mutex;
  std::condition_variable _condition;
  size_t _n, _count, _cycle;
  std::function<void()> _completionFunction;
};
}  // namespace aocommon

#endif
