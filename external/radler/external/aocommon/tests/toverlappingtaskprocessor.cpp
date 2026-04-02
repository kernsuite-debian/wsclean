#include <aocommon/overlappingtaskprocessor.h>

#include <mutex>
#include <numeric>
#include <random>
#include <thread>
#include <vector>

#include <boost/test/unit_test.hpp>

using aocommon::OverlappingTaskProcessor;

BOOST_AUTO_TEST_SUITE(overlappingtaskprocessor)

struct TestTask {
  size_t lane;
  size_t task;
};

struct TestWorkPhase {
  size_t lane;
  size_t task;
  size_t phase;
};

BOOST_AUTO_TEST_CASE(multiple_threads) {
  const size_t kNumSimulatedLanes = 100;
  const size_t kNumTasksPerLane = 10;
  const size_t kNumProcessingPhasesPerTask = 3;
  const size_t kNumThreads = 3;

  std::random_device random_seed;
  std::mt19937 random_generator(random_seed());
  std::uniform_int_distribution<> random_size_small(1, 50);
  std::uniform_int_distribution<> random_size_large(50, 150);

  // Single queue with threads for processing all work.
  aocommon::TaskQueue<std::function<void()>> run_task_queue;
  std::vector<std::thread> run_thread_pool;
  run_thread_pool.reserve(kNumThreads);
  for (size_t i = 0; i < kNumThreads; ++i) {
    run_thread_pool.emplace_back([&] {
      std::function<void()> operation;
      while (run_task_queue.Pop(operation)) {
        operation();
      }
    });
  }
  // Wrapped in an overlapped processor to allow overlaps.
  OverlappingTaskProcessor overlapping_processor(run_task_queue);

  std::mutex processing_order_mutex;
  std::vector<TestWorkPhase> processing_order;
  std::vector<std::mutex> phase_mutex(kNumProcessingPhasesPerTask);

  // Track concurrency.
  size_t max_concurrent_phases = 0;
  std::atomic<size_t> current_concurrent_phases = 0;

  // Simulate reading in multiple chunks for multiple tasks with a lane for each
  // task. Call Process() for each task. Limit Process() calls to maximum 2 in
  // parallel.
  std::counting_semaphore<2> task_semaphore(2);
  for (size_t lane = 0; lane < kNumSimulatedLanes; ++lane) {
    task_semaphore.acquire();
    aocommon::Lane<TestTask> task_lane(kNumTasksPerLane / 2);
    std::thread([&]() {
      for (size_t task = 0; task < kNumTasksPerLane; ++task) {
        TestTask chunk_data(lane, task);
        // Simulate data read time.
        std::this_thread::sleep_for(
            std::chrono::microseconds(random_size_small(random_seed)));
        task_lane.write(std::move(chunk_data));
      }
      task_lane.write_end();
    }).detach();
    overlapping_processor.Process<TestTask>(
        task_lane,
        [&](TestTask&& chunk_data, size_t chunk_index,
            std::binary_semaphore& processing_order_semaphore) mutable {
          size_t lane = chunk_data.lane;
          size_t task = chunk_data.task;
          current_concurrent_phases++;
          {
            std::lock_guard<std::mutex> lock_processing_order(
                processing_order_mutex);
            max_concurrent_phases = std::max(max_concurrent_phases,
                                             size_t{current_concurrent_phases});
            processing_order.push_back(TestWorkPhase(lane, task, 0));
          }
          {
            // Ensure lock transition between phases without a race in between.
            // When we assign to the unique_ptr the old lock will only be
            // destroyed after the new one is created so we are guaranteed to
            // gain the new lock before other threads can obtain the old one.
            std::unique_ptr<std::lock_guard<std::mutex>> phase_lock;
            for (size_t phase = 1; phase <= kNumProcessingPhasesPerTask;
                 ++phase) {
              phase_lock = std::make_unique<std::lock_guard<std::mutex>>(
                  phase_mutex[phase - 1]);
              if (phase == 1) {
                processing_order_semaphore.release();
              }
              // Simulate task processing time.
              std::this_thread::sleep_for(
                  std::chrono::microseconds(random_size_large(random_seed)));
              {
                std::lock_guard<std::mutex> lock_processing_order(
                    processing_order_mutex);
                processing_order.push_back(TestWorkPhase(lane, task, phase));
              }
            }
          }
          current_concurrent_phases--;
          task_semaphore.release();
        },
        "TestTask");
  }

  // Wait for all processing to complete and cleanup threads.
  run_task_queue.WaitForIdle(kNumThreads);
  run_task_queue.Finish();
  for (std::thread& thread : run_thread_pool) {
    thread.join();
  }
  // There should have been 2 tasks running concurrently throughout most of this
  // test.
  BOOST_CHECK_EQUAL(max_concurrent_phases, 2);
  // Ensure later tasks of later lanes never preceed earlier tasks of earlier
  // lanes.
  BOOST_CHECK(
      std::is_sorted(processing_order.begin(), processing_order.end(),
                     [](const TestWorkPhase& p1, const TestWorkPhase& p2) {
                       return p1.lane < p2.lane && p1.task < p2.task;
                     }));
  for (size_t lane = 0; lane < kNumSimulatedLanes; ++lane) {
    std::vector<TestWorkPhase> lane_processing_order;
    std::copy_if(processing_order.begin(), processing_order.end(),
                 std::back_inserter(lane_processing_order),
                 [&](TestWorkPhase item) { return item.lane == lane; });
    // Ensure phases within a single task are always processed in order.
    for (size_t task = 0; task < kNumTasksPerLane; ++task) {
      std::vector<TestWorkPhase> task_processing_order;
      std::copy_if(lane_processing_order.begin(), lane_processing_order.end(),
                   std::back_inserter(task_processing_order),
                   [&](TestWorkPhase item) { return item.task == task; });
      BOOST_CHECK(std::is_sorted(
          task_processing_order.begin(), task_processing_order.end(),
          [](const TestWorkPhase& p1, const TestWorkPhase& p2) {
            return p1.phase < p2.phase;
          }));
    }
    // Ensure later phases of later tasks never preceed earlier phases of
    // earlier tasks.
    BOOST_CHECK(std::is_sorted(
        lane_processing_order.begin(), lane_processing_order.end(),
        [](const TestWorkPhase& p1, const TestWorkPhase& p2) {
          return p1.task < p2.task && p1.phase < p2.phase;
        }));
    std::vector<TestWorkPhase> lane_start_processing_order;
    std::copy_if(lane_processing_order.begin(), lane_processing_order.end(),
                 std::back_inserter(lane_start_processing_order),
                 [&](TestWorkPhase item) { return item.phase == 0; });
    // Ensure tasks within a lane are always started in order.
    BOOST_CHECK(std::is_sorted(
        lane_start_processing_order.begin(), lane_start_processing_order.end(),
        [](const TestWorkPhase& p1, const TestWorkPhase& p2) {
          return p1.task < p2.task;
        }));
  }
}

BOOST_AUTO_TEST_SUITE_END()
