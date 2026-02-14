#include <aocommon/threadpool.h>

#include <boost/test/unit_test.hpp>

#include <set>
#include <vector>

using aocommon::ThreadPool;

BOOST_AUTO_TEST_SUITE(thread_pool)

BOOST_AUTO_TEST_CASE(restart_threads) {
  ThreadPool pool;
  BOOST_CHECK_EQUAL(pool.NThreads(), 1);

  pool.SetNThreads(4);
  BOOST_CHECK_EQUAL(pool.NThreads(), 4);

  pool.SetNThreads(1);
  BOOST_CHECK_EQUAL(pool.NThreads(), 1);
}

BOOST_AUTO_TEST_CASE(single_threaded) {
  ThreadPool pool;
  pool.SetNThreads(1);
  size_t run_count = 0;
  std::thread::id executing_thread_id;
  auto run = [&](size_t thread_index) {
    BOOST_CHECK_EQUAL(thread_index, 0);
    ++run_count;
    executing_thread_id = std::this_thread::get_id();
  };
  pool.ExecuteInParallel(run);
  BOOST_CHECK_EQUAL(run_count, 1);
  BOOST_CHECK(executing_thread_id == std::this_thread::get_id());

  executing_thread_id = std::thread::id();
  pool.ExecuteInParallel(run);
  BOOST_CHECK_EQUAL(run_count, 2);
  BOOST_CHECK(executing_thread_id == std::this_thread::get_id());
}

BOOST_AUTO_TEST_CASE(multi_threaded) {
  ThreadPool pool;
  constexpr size_t kNThreads = 10;
  pool.SetNThreads(kNThreads);
  std::vector<int> run_counts(kNThreads, 0);
  std::set<std::thread::id> thread_ids;
  std::mutex mutex;
  auto run = [&](size_t thread_index) {
    std::lock_guard lock(mutex);
    thread_ids.insert(std::this_thread::get_id());
    ++run_counts[thread_index];
  };
  pool.ExecuteInParallel(run);
  const std::vector<int> ref_a(kNThreads, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(run_counts.begin(), run_counts.end(),
                                ref_a.begin(), ref_a.end());
  BOOST_CHECK_EQUAL(thread_ids.size(), kNThreads);

  thread_ids.clear();
  pool.ExecuteInParallel(run);
  const std::vector<int> ref_b(kNThreads, 2);
  BOOST_CHECK_EQUAL_COLLECTIONS(run_counts.begin(), run_counts.end(),
                                ref_b.begin(), ref_b.end());
  BOOST_CHECK_EQUAL(thread_ids.size(), kNThreads);
}

BOOST_AUTO_TEST_CASE(non_blocking_single_threaded) {
  ThreadPool pool;
  pool.SetNThreads(1);
  bool is_executed = false;
  // With only one thread, the function should never be executed
  pool.StartParallelExecution([&](size_t) {
    is_executed = true;
    while (true)
      ;
  });
  BOOST_CHECK(!is_executed);
  pool.FinishParallelExecution();
  BOOST_CHECK(!is_executed);
}

BOOST_AUTO_TEST_CASE(non_blocking_multi_threaded) {
  ThreadPool pool;
  constexpr size_t kNThreads = 10;
  pool.SetNThreads(kNThreads);
  std::vector<int> run_counts(kNThreads, 0);
  std::set<std::thread::id> thread_ids;
  std::mutex mutex;
  auto run = [&](size_t thread_index) {
    std::lock_guard lock(mutex);
    thread_ids.insert(std::this_thread::get_id());
    ++run_counts[thread_index];
  };
  pool.StartParallelExecution(run);
  pool.FinishParallelExecution();
  std::vector<int> ref_a(kNThreads, 1);
  ref_a[0] = 0;
  BOOST_CHECK_EQUAL_COLLECTIONS(run_counts.begin(), run_counts.end(),
                                ref_a.begin(), ref_a.end());
  BOOST_CHECK_EQUAL(thread_ids.size(), kNThreads - 1);

  thread_ids.clear();
  pool.StartParallelExecution(run);
  pool.FinishParallelExecution();
  std::vector<int> ref_b(kNThreads, 2);
  ref_b[0] = 0;
  BOOST_CHECK_EQUAL_COLLECTIONS(run_counts.begin(), run_counts.end(),
                                ref_b.begin(), ref_b.end());
  BOOST_CHECK_EQUAL(thread_ids.size(), kNThreads - 1);
}

BOOST_AUTO_TEST_SUITE_END()
