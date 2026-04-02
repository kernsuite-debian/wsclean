#include <aocommon/taskqueue.h>

#include <atomic>
#include <numeric>
#include <random>
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
      int popped = 0;
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

// This test is not 100% deterministic, it is possible (though unlikely) that
// the threads of the 'Pop' for loop execute and finish before the threads of
// the 'TryPop' for loop ever get a chance to run. If this occurs the test may
// fail despite there being no actual code error. To avoid this each worker
// thread calls sleep with a small random interval, to allow the OS an
// opportunity to context switch to other threads. This should ensure that all
// threads in this test get a chance to execute and should avoid the failure
// described above. A random instead of fixed interval is used to better stress
// test the concurrency of the data structure by having different threads wake
// up/operate in different order.
//
// If you are reading this the test may have failed intermittently despite these
// precautions. If so it might become necessary to simplify or rework this test
// in other ways to try and avoid intermittent failures.
BOOST_AUTO_TEST_CASE(multiple_threads_try_pop) {
  std::random_device random_seed;
  std::mt19937 random_generator(random_seed());
  std::uniform_int_distribution<> random_number(1, 5);

  std::vector<int> values(10240);
  std::iota(values.begin(), values.end(), 0);

  TaskQueue<int> queue;
  for (const int& value : values) {
    queue.Emplace(value);
  }

  // Test that we can get tasks via 'Pop' and 'TryPop' concurrently without
  // issue.
  std::vector<std::thread> pop_threads;
  std::atomic<size_t> n_popped = 0;
  for (size_t i = 0; i < 3; ++i) {
    pop_threads.emplace_back([&] {
      int popped_value = 0;
      while (queue.Pop(popped_value)) {
        n_popped++;
        // Prevent thread starvation.
        // See comment at top of test for further explanation.
        std::this_thread::sleep_for(
            std::chrono::milliseconds(random_number(random_seed)));
      }
    });
  }
  std::atomic<size_t> n_try_popped = 0;
  for (size_t i = 0; i < 3; ++i) {
    pop_threads.emplace_back([&] {
      int popped_value;
      while (queue.TryPop(popped_value)) {
        n_try_popped++;
        // Prevent thread starvation.
        // See comment at top of test for further explanation.
        std::this_thread::sleep_for(
            std::chrono::milliseconds(random_number(random_seed)));
      }
      BOOST_CHECK(!queue.Pop(popped_value));
    });
  }
  queue.WaitForIdle(pop_threads.size());

  BOOST_CHECK_GT(n_popped, 0);
  BOOST_CHECK_GT(n_try_popped, 0);
  BOOST_CHECK_EQUAL(n_popped + n_try_popped, values.size());

  queue.Finish();
  for (std::thread& thread : pop_threads) thread.join();
}

// This test is not 100% deterministic, it is possible (though unlikely) that
// the threads of the 'Pop' for loop execute and finish before the threads of
// the 'TryPopN' for loop ever get a chance to run. See comment at top of
// 'multiple_threads_try_pop' test for further explanation.
BOOST_AUTO_TEST_CASE(multiple_threads_try_pop_n) {
  std::random_device random_seed;
  std::mt19937 random_generator(random_seed());
  std::uniform_int_distribution<> random_number(1, 5);

  std::vector<int> values(10240);
  std::iota(values.begin(), values.end(), 0);

  TaskQueue<int> queue;
  for (const int& value : values) {
    queue.Emplace(value);
  }

  // Test that we can get tasks via 'Pop' and 'TryPop' concurrently without
  // issue.
  // Also that TryPopN pops the correct amount of tasks and that they are
  // ordered.
  std::vector<std::thread> pop_threads;
  std::atomic<size_t> n_popped = 0;
  for (size_t i = 0; i < 4; ++i) {
    pop_threads.emplace_back([&] {
      int popped_value;
      while (queue.Pop(popped_value)) {
        n_popped++;
        // Prevent thread starvation.
        // See comment at top of test for further explanation.
        std::this_thread::sleep_for(
            std::chrono::milliseconds(random_number(random_seed)));
      }
    });
  }
  std::atomic<size_t> n_try_popped = 0;
  std::atomic<size_t> n_additional_popped = 0;
  for (size_t i = 0; i < 3; ++i) {
    pop_threads.emplace_back([&] {
      std::vector<int> popped_values;
      size_t read_size = random_number(random_generator);
      while (queue.TryPopN(popped_values, read_size)) {
        BOOST_CHECK_EQUAL(popped_values.size(), read_size);
        std::vector<int> expected_values(read_size);
        std::iota(expected_values.begin(), expected_values.end(),
                  popped_values[0]);
        BOOST_CHECK(popped_values == expected_values);
        n_try_popped += read_size;
        read_size = random_number(random_generator);
        // Prevent thread starvation.
        // See comment at top of test for further explanation.
        std::this_thread::sleep_for(
            std::chrono::milliseconds(random_number(random_seed)));
      }
      int popped_value = 0;
      // Note that while it might seem like this call can be removed to simplify
      // the test, it is mandatory that all worker threads end up blocking
      // inside a call to `Pop` otherwise `WaitForIdle` will not function
      // correctly.
      while (queue.Pop(popped_value)) {
        n_additional_popped++;
      }
    });
  }
  queue.WaitForIdle(pop_threads.size());

  BOOST_CHECK_GT(n_popped, 0);
  BOOST_CHECK_GT(n_try_popped, 0);
  BOOST_CHECK_EQUAL(n_popped + n_additional_popped + n_try_popped,
                    values.size());

  // Test that read past end of available tasks returns false.
  {
    int popped_value = 0;
    BOOST_CHECK(!queue.TryPop(popped_value));
  }
  {
    std::vector<int> popped_values;
    BOOST_CHECK(!queue.TryPopN(popped_values, 1));
    BOOST_CHECK(!queue.TryPopN(popped_values, 2));
    BOOST_CHECK(!queue.TryPopN(popped_values, 10));
  }

  queue.Finish();
  for (std::thread& thread : pop_threads) thread.join();
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
