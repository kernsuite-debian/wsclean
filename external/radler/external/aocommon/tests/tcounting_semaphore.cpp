#include <atomic>

#include <aocommon/counting_semaphore.h>
#include <aocommon/dynamicfor.h>

#include <boost/test/unit_test.hpp>

using aocommon::CountingSemaphore;
using aocommon::DynamicFor;
using aocommon::ScopedCountingSemaphoreLock;

BOOST_AUTO_TEST_SUITE(counting_semaphore)

BOOST_AUTO_TEST_CASE(acquire_and_release) {
  aocommon::ThreadPool::GetInstance().SetNThreads(20);
  const size_t kLimit = 4;
  CountingSemaphore semaphore(kLimit);
  DynamicFor<size_t> loop;
  std::atomic<size_t> max_count = 0;
  std::atomic<size_t> run_count = 0;
  // A value of 2000 is high enough for this test to (sometimes) produce errors
  // when leaving out the semaphore
  loop.Run(0, 2000, [&](size_t) {
    ScopedCountingSemaphoreLock lock(semaphore);
    run_count.fetch_add(1);
    max_count = std::max(max_count.load(), run_count.load());
    run_count.fetch_sub(1);
  });
  BOOST_CHECK_LE(max_count, kLimit);
  BOOST_CHECK_GE(max_count, 0);
}

BOOST_AUTO_TEST_SUITE_END()
