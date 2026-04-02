#include <atomic>
#include <iostream>
#include <cmath>

#include <unistd.h>  // for usleep

#include <aocommon/recursivefor.h>

#include <boost/test/unit_test.hpp>

using namespace aocommon;

namespace {
void RunSingleLoop(bool nested = false) {
  std::mutex mutex;
  std::vector<size_t> counts(10, 0);
  auto function = [&](size_t iter, size_t) {
    std::unique_lock<std::mutex> lock(mutex);
    counts[iter]++;
  };
  if (nested) {
    RecursiveFor::NestedRun(0, 10, function);
  } else {
    RecursiveFor loop;
    loop.Run(0, 10, function);
  }
  std::vector<size_t> ref(10, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}
}  // namespace

BOOST_AUTO_TEST_SUITE(recursive_for)

BOOST_AUTO_TEST_CASE(construction_single_threaded) {
  ThreadPool::GetInstance().SetNThreads(1);
  BOOST_CHECK(RecursiveFor::GetInstance() == nullptr);
  RecursiveFor loop;
  BOOST_CHECK(RecursiveFor::GetInstance() == &loop);
}

BOOST_AUTO_TEST_CASE(construction_multi_threaded) {
  ThreadPool::GetInstance().SetNThreads(4);
  BOOST_CHECK(RecursiveFor::GetInstance() == nullptr);
  RecursiveFor loop;
  BOOST_CHECK(RecursiveFor::GetInstance() == &loop);
}

BOOST_AUTO_TEST_CASE(serial_single_loop) {
  ThreadPool::GetInstance().SetNThreads(1);
  RunSingleLoop();
}

BOOST_AUTO_TEST_CASE(parallel_single_loop) {
  ThreadPool::GetInstance().SetNThreads(4);
  RunSingleLoop();
}

BOOST_AUTO_TEST_CASE(recursive) {
  ThreadPool::GetInstance().SetNThreads(4);
  RecursiveFor loop;
  std::mutex mutex;
  std::vector<size_t> counts(800, 0);
  loop.Run(0, 100, [&](size_t iter1, size_t) {
    loop.Run(0, 8, [&](size_t iter2, size_t) {
      std::unique_lock<std::mutex> lock(mutex);
      counts[iter1 + iter2 * 100]++;
    });
  });
  std::vector<size_t> ref(800, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}

BOOST_AUTO_TEST_CASE(slow_task) {
  constexpr size_t kNThread = 4;
  ThreadPool::GetInstance().SetNThreads(kNThread);
  std::mutex mutex;
  RecursiveFor loop;
  std::vector<size_t> counts(5 * kNThread, 0);
  loop.Run(0, 5 * kNThread, [&](size_t iter, size_t) {
    usleep(1000);
    std::unique_lock<std::mutex> lock(mutex);
    counts[iter]++;
  });
  std::vector<size_t> ref(5 * kNThread, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}

BOOST_AUTO_TEST_CASE(nested_run) {
  ThreadPool::GetInstance().SetNThreads(4);
  // Call NestedRun() without an active recursive for:
  RunSingleLoop(true);
  // Same, but now with an active recursive for:
  const RecursiveFor outer_recursive_for;
  RunSingleLoop(true);
}

BOOST_AUTO_TEST_CASE(reuse) {
  ThreadPool::GetInstance().SetNThreads(4);
  RecursiveFor recursive_for;
  for (size_t i = 0; i != 100; ++i) {
    volatile double x = 0.1;
    recursive_for.Run(0, 100, [&](size_t, size_t) { (void)std::sin(x); });
  }
}

BOOST_AUTO_TEST_CASE(exception_single_threaded) {
  ThreadPool::GetInstance().SetNThreads(1);
  RecursiveFor loop;
  auto test_call = [&loop]() {
    loop.Run(0, 1000, [](size_t) { throw std::exception(); });
  };
  BOOST_CHECK_THROW(test_call(), std::exception);
}

BOOST_AUTO_TEST_CASE(exception_multi_threaded) {
  ThreadPool::GetInstance().SetNThreads(4);
  RecursiveFor loop;
  auto test_call = [&loop](size_t tested_thread) {
    std::atomic<bool> found(false);
    loop.Run(0, 1000, [tested_thread, &found](size_t, size_t thread_index) {
      if (thread_index == tested_thread) {
        found = true;
        throw std::exception();
      } else {
        // Non-throwing threads are delayed to make sure that the
        // requested thread gets to run an iteration. Otherwise, other
        // threads might process all iterations before the requested
        // thread starts, causing it to never throw.
        while (!found) {
        }
      }
    });
  };
  // Exceptions in thread 0 and other threads are handled differently,
  // so test both
  BOOST_CHECK_THROW(test_call(0), std::exception);
  BOOST_CHECK_THROW(test_call(3), std::exception);
}

BOOST_AUTO_TEST_SUITE_END()
