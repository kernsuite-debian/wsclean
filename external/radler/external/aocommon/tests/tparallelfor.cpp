#include <aocommon/parallelfor.h>

#include <iostream>

#include <cmath>
#include <exception>
#include <mutex>
#include <vector>

#include <unistd.h>  // for sleep

#include <boost/test/unit_test.hpp>

using aocommon::ParallelFor;

namespace {
void BasicTest(ParallelFor<size_t>& loop, bool with_thread_id) {
  std::mutex mutex;
  const size_t n_threads = loop.NThreads();
  std::vector<size_t> counts(10, 0);
  if (with_thread_id) {
    loop.Run(0, 10, [&](size_t iter, size_t thread_id) {
      // If ParallelFor behaves properly, this lock isn't required, so this is
      // just to make sure the test behaves as expected even when ParallelFor
      // doesn't. This also holds for some of the later tests.
      std::unique_lock<std::mutex> lock(mutex);
      BOOST_CHECK_LT(thread_id, n_threads);
      counts[iter]++;
    });
  } else {
    loop.Run(0, 10, [&](size_t iter) {
      std::unique_lock<std::mutex> lock(mutex);
      counts[iter]++;
    });
  }
  const std::vector<size_t> ref(10, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}
}  // namespace

BOOST_AUTO_TEST_SUITE(parallelfor)

BOOST_AUTO_TEST_CASE(unused) { BOOST_CHECK_NO_THROW(ParallelFor<size_t>(4);); }

BOOST_AUTO_TEST_CASE(single_thread_with_thread_id) {
  ParallelFor<size_t> loop(1);
  BasicTest(loop, true);
}

BOOST_AUTO_TEST_CASE(single_thread_without_thread_id) {
  ParallelFor<size_t> loop(1);
  BasicTest(loop, false);
}

BOOST_AUTO_TEST_CASE(parallel_with_thread_id) {
  ParallelFor<size_t> loop(4);
  BasicTest(loop, true);
}

BOOST_AUTO_TEST_CASE(parallel_without_thread_id) {
  ParallelFor<size_t> loop(4);
  BasicTest(loop, false);
}

BOOST_AUTO_TEST_CASE(reuse_with_thread_id) {
  std::vector<size_t> counts(20, 0);
  std::mutex mutex;
  ParallelFor<size_t> loop(4);
  loop.Run(0, 10, [&](size_t iter, size_t thread_id) {
    // See BasicTest() for an explanation
    std::unique_lock<std::mutex> lock(mutex);
    BOOST_CHECK_LT(thread_id, 4);
    counts[iter]++;
  });
  std::vector<size_t> ref(20, 0);
  std::fill(ref.begin(), ref.begin() + 10, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());

  loop.Run(10, 20, [&](size_t iter, size_t thread_id) {
    // See BasicTest() for an explanation
    std::unique_lock<std::mutex> lock(mutex);
    BOOST_CHECK_LT(thread_id, 4);
    counts[iter]++;
  });
  ref = std::vector<size_t>(20, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}

BOOST_AUTO_TEST_CASE(reuse_without_thread_id) {
  std::vector<size_t> counts(20, 0);
  std::mutex mutex;
  ParallelFor<size_t> loop(4);
  loop.Run(0, 10, [&](size_t iter) {
    // See BasicTest() for an explanation
    std::unique_lock<std::mutex> lock(mutex);
    counts[iter]++;
  });
  std::vector<size_t> ref(20, 0);
  std::fill(ref.begin(), ref.begin() + 10, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());

  loop.Run(10, 20, [&](size_t iter) {
    std::unique_lock<std::mutex> lock(mutex);
    counts[iter]++;
  });
  ref = std::vector<size_t>(20, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}

BOOST_AUTO_TEST_CASE(exception_single_threaded) {
  ParallelFor<size_t> loop(1);
  auto test_call = [&loop]() {
    loop.Run(0, 1000, [](size_t) { throw std::exception(); });
  };
  BOOST_CHECK_THROW(test_call(), std::exception);

  // Object still behaving properly?
  BasicTest(loop, true);
}

BOOST_AUTO_TEST_CASE(exception_multi_threaded) {
  ParallelFor<size_t> loop(4);
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

  BasicTest(loop, true);

  BOOST_CHECK_THROW(test_call(3), std::exception);

  BasicTest(loop, true);
}

BOOST_AUTO_TEST_CASE(multiple_exceptions) {
  const size_t kNThreads = 32;
  ParallelFor<size_t> loop(kNThreads);
  auto test_call = [&loop]() {
    aocommon::Barrier barrier(kNThreads);
    loop.Run(0, 1000, [&barrier](size_t) {
      barrier.wait();
      throw std::exception();
    });
  };

  BOOST_CHECK_THROW(test_call(), std::exception);
  BasicTest(loop, true);
}

BOOST_AUTO_TEST_SUITE_END()
