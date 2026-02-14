#include <aocommon/taskqueue.h>

#include <atomic>
#include <thread>
#include <vector>

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

using aocommon::TaskQueue;

BOOST_AUTO_TEST_SUITE(task_queue)

BOOST_DATA_TEST_CASE(single_thread, boost::unit_test::data::make({false, true}),
                     use_wait) {
  const std::vector<int> kValues{42, 43, 44, 45};
  constexpr int kDummyValue = 142;

  // Using a unique pointer ensures that TaskQueue cannot copy tasks.
  TaskQueue<std::unique_ptr<int>> queue;
  for (const int& value : kValues) {
    queue.Emplace(std::make_unique<int>(value));
  }

  for (const int& value : kValues) {
    std::unique_ptr<int> popped;
    BOOST_TEST(queue.Pop(popped));
    BOOST_REQUIRE(popped);
    BOOST_TEST(*popped == value);
  }

  if (use_wait) queue.WaitForIdle(0);

  queue.Finish();
  auto dummy = std::make_unique<int>(kDummyValue);
  BOOST_TEST(!queue.Pop(dummy));
  BOOST_REQUIRE(dummy);
  BOOST_TEST(*dummy == kDummyValue);
}

BOOST_DATA_TEST_CASE(multiple_threads_pop,
                     boost::unit_test::data::make({false, true}), use_wait) {
  const std::vector<int> kValues{42, 43, 44, 45};
  const size_t kLimit = 2;
  TaskQueue<int> queue{kLimit};
  std::mutex mutex;
  std::condition_variable notify;
  int popped_in_thread = 0;

  std::vector<std::thread> pop_threads;
  for (size_t i = 0; i < kValues.size(); ++i) {
    pop_threads.emplace_back([&] {
      int popped;
      const bool result = queue.Pop(popped);
      std::lock_guard<std::mutex> lock(mutex);
      BOOST_TEST_REQUIRE(result);
      popped_in_thread = popped;
      notify.notify_one();
    });
  }

  if (use_wait) queue.WaitForIdle(pop_threads.size());

  for (const int& value : kValues) {
    popped_in_thread = 0;
    queue.Emplace(value);
    std::unique_lock<std::mutex> lock(mutex);
    while (popped_in_thread == 0) notify.wait(lock);
    BOOST_TEST(popped_in_thread == value);
  }

  for (std::thread& thread : pop_threads) thread.join();

  if (use_wait) queue.WaitForIdle(0);
}

BOOST_AUTO_TEST_CASE(multiple_threads_done) {
  constexpr size_t kNThreads = 42;
  constexpr int kDummyValue = 142;
  TaskQueue<int> queue;

  std::mutex mutex;
  std::vector<std::thread> threads;
  for (size_t i = 0; i < kNThreads; ++i) {
    threads.emplace_back([&] {
      int dummy = kDummyValue;
      std::lock_guard<std::mutex> lock(mutex);
      BOOST_TEST(!queue.Pop(dummy));
      BOOST_TEST(dummy == kDummyValue);
    });
  }

  queue.Finish();

  // Joining the threads also tests that all threads are done.
  for (std::thread& thread : threads) thread.join();
}

BOOST_AUTO_TEST_CASE(wait_for_idle) {
  // Test that WaitForIdle really waits until kNThreads call Pop().
  const size_t kNThreads = 42;

  TaskQueue<int> queue;

  std::atomic<bool> waiting = false;
  std::atomic<bool> done_waiting = false;
  std::thread wait_thread([&] {
    waiting = true;
    queue.WaitForIdle(kNThreads);
    done_waiting = true;
  });
  // Wait until wait_thread starts waiting.
  while (!waiting) std::this_thread::yield();

  std::vector<std::thread> pop_threads;
  for (size_t i = 0; i < kNThreads; ++i) {
    BOOST_TEST(waiting);
    BOOST_TEST(!done_waiting);
    std::atomic<bool> popping = false;
    pop_threads.emplace_back([&] {
      popping = true;
      int dummy;
      queue.Pop(dummy);
    });
    // Wait until the thread starts popping.
    while (!popping) std::this_thread::yield();
  }

  // Wait until wait_thread stops waiting.
  while (!done_waiting) std::this_thread::yield();

  wait_thread.join();
  queue.Finish();
  for (std::thread& pop_thread : pop_threads) pop_thread.join();
}

BOOST_AUTO_TEST_SUITE_END()