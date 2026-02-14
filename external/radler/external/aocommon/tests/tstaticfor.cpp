#include <aocommon/recursivefor.h>
#include <aocommon/staticfor.h>

#include <atomic>
#include <set>
#include <mutex>

#include <unistd.h>  // for sleep

#include <boost/test/unit_test.hpp>

using aocommon::RecursiveFor;
using aocommon::StaticFor;
using aocommon::ThreadPool;

BOOST_AUTO_TEST_SUITE(staticfor)

BOOST_AUTO_TEST_CASE(run) {
  ThreadPool::GetInstance().SetNThreads(4);
  StaticFor<size_t> loop;
  std::mutex mutex;
  std::vector<size_t> counts(10, 0);
  loop.Run(0, 10, [&](size_t a, size_t b) {
    for (size_t iter = a; iter != b; ++iter) {
      std::unique_lock<std::mutex> lock(mutex);
      counts[iter]++;
    }
  });

  std::vector<size_t> ref(10, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}

BOOST_AUTO_TEST_CASE(empty_loop) {
  ThreadPool::GetInstance().SetNThreads(4);
  StaticFor<size_t> loop;
  std::atomic<bool> fail = false;
  loop.Run(10, 10, [&](size_t a, size_t b) {
    if (a != b) fail = true;
  });
  BOOST_CHECK(!fail);
}

BOOST_AUTO_TEST_CASE(single_threaded) {
  ThreadPool::GetInstance().SetNThreads(1);
  StaticFor<size_t> loop;
  std::vector<size_t> counts(10, 0);
  loop.Run(0, 10, [&](size_t a, size_t b) {
    for (size_t iter = a; iter != b; ++iter) {
      counts[iter]++;
    }
  });

  std::vector<size_t> ref(10, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}

BOOST_AUTO_TEST_CASE(resume_run) {
  std::vector<size_t> counts(20, 0);
  std::mutex mutex;
  ThreadPool::GetInstance().SetNThreads(40);
  StaticFor<size_t> loop;
  loop.Run(0, 10, [&](size_t a, size_t b) {
    for (size_t iter = a; iter != b; ++iter) {
      std::unique_lock<std::mutex> lock(mutex);
      counts[iter]++;
    }
  });
  std::vector<size_t> ref(20, 0);
  std::fill(ref.begin(), ref.begin() + 10, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());

  loop.Run(10, 20, [&](size_t a, size_t b) {
    for (size_t iter = a; iter != b; ++iter) {
      std::unique_lock<std::mutex> lock(mutex);
      counts[iter]++;
    }
  });
  ref = std::vector<size_t>(20, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}

BOOST_AUTO_TEST_CASE(run_with_thread_id) {
  ThreadPool::GetInstance().SetNThreads(4);
  StaticFor<size_t> loop;
  std::mutex mutex;
  std::vector<size_t> counts(10, 0);
  std::vector<size_t> threads(4, 0);
  loop.Run(0, 10, [&](size_t a, size_t b, size_t t) {
    for (size_t iter = a; iter != b; ++iter) {
      std::unique_lock<std::mutex> lock(mutex);
      counts[iter]++;
    }
    std::unique_lock<std::mutex> lock(mutex);
    threads[t]++;
  });

  std::vector<size_t> ref(10, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());

  // Not all threads might actually be used, because if one thread
  // finishes before a second thread is starting, the first
  // thread is used to perform the next loop.
  // Therefore all we can check is whether there weren't more than
  // 4 blocks:
  for (size_t i = 0; i != threads.size(); ++i) BOOST_CHECK_LT(threads[i], 5);
}

BOOST_AUTO_TEST_CASE(inside_recursive_for) {
  const size_t kNThreads = 7;
  ThreadPool::GetInstance().SetNThreads(kNThreads);
  RecursiveFor recursive_for;
  std::mutex mutex;
  std::map<std::pair<size_t, size_t>, size_t> counter;
  recursive_for.Run(0, 10, [&](size_t outer, size_t) {
    StaticFor<size_t> nested_for;
    nested_for.Run(0, 10, [&, outer](size_t start, size_t end) {
      for (size_t inner = start; inner != end; ++inner) {
        std::lock_guard lock(mutex);
        counter[std::make_pair(outer, inner)]++;
      }
    });
  });
  BOOST_CHECK_EQUAL(counter.size(), 100);
  for (const std::pair<const std::pair<size_t, size_t>, size_t>& count :
       counter) {
    BOOST_CHECK_LT(count.first.first, 10);   // outer
    BOOST_CHECK_LT(count.first.second, 10);  // inner
    BOOST_CHECK_EQUAL(count.second, 1);      // count
  }
}

BOOST_AUTO_TEST_CASE(constrained_run_without_id) {
  ThreadPool::GetInstance().SetNThreads(4);
  StaticFor<size_t> loop;
  std::mutex mutex;
  std::vector<size_t> counts(10, 0);
  std::set<std::thread::id> thread_ids;
  const size_t kMaxThreads = 3;
  loop.ConstrainedRun(0, 10, kMaxThreads, [&](size_t a, size_t b) {
    for (size_t iter = a; iter != b; ++iter) {
      std::scoped_lock<std::mutex> lock(mutex);
      counts[iter]++;
    }
    std::scoped_lock<std::mutex> lock(mutex);
    thread_ids.emplace(std::this_thread::get_id());
  });

  std::vector<size_t> ref(10, 1);
  BOOST_CHECK_EQUAL(thread_ids.size(), kMaxThreads);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}

BOOST_AUTO_TEST_CASE(constrained_run_with_id) {
  ThreadPool::GetInstance().SetNThreads(4);
  StaticFor<size_t> loop;
  std::mutex mutex;
  std::vector<size_t> counts(10, 0);
  std::set<std::thread::id> thread_ids;
  std::set<size_t> thread_indices;
  const size_t kMaxThreads = 3;
  loop.ConstrainedRun(0, 10, kMaxThreads,
                      [&](size_t a, size_t b, size_t thread_index) {
                        for (size_t iter = a; iter != b; ++iter) {
                          std::scoped_lock<std::mutex> lock(mutex);
                          counts[iter]++;
                        }
                        std::scoped_lock<std::mutex> lock(mutex);
                        thread_ids.emplace(std::this_thread::get_id());
                        thread_indices.emplace(thread_index);
                      });

  BOOST_CHECK_EQUAL(thread_ids.size(), kMaxThreads);
  BOOST_CHECK_EQUAL(thread_indices.size(), kMaxThreads);
  const std::set<size_t> expected_indices{0, 1, 2};
  BOOST_CHECK_EQUAL_COLLECTIONS(thread_indices.begin(), thread_indices.end(),
                                expected_indices.begin(),
                                expected_indices.end());
  std::vector<size_t> expected_counts(10, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(),
                                expected_counts.begin(), expected_counts.end());
}

BOOST_AUTO_TEST_SUITE_END()
