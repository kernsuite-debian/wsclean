#ifndef AOCOMMON_OVERLAPPING_TASK_PROCESSOR_H_
#define AOCOMMON_OVERLAPPING_TASK_PROCESSOR_H_

#include <functional>
#include <semaphore>

#include "lane.h"
#include "logger.h"
#include "taskqueue.h"

namespace aocommon {
/**
 * Read chunks from a lane and ensure that all chunks are processed in order
 * with at most 2 chunks processed simultaneously in an overlapped manner.
 * This enables overlapped processing across two chunks which can be used to
 * help ensure more continuous resource (CPU) usage.
 * Guarantee:
 *   1. Chunks are processed in order.
 *   2. Never more than two chunks are processed simultaneously.
 *
 * This class can read chunks from a second lane as well, in which case in
 * addition the the already described functionality it will also ensure that
 * tasks from lane1 are all processed before chunks from lane2 allowed to
 * overlap its first task.
 * Guarantee:
 *   1. Lanes are processed in order.
 *   2. The first chunk of lane2 can only begin processing when all but 1 chunk
 *   of lane1 is done processing.
 * NB! At most two concurrent @ref Process() calls are allowed.
 *
 * It is the callers responsibility to ensure that their @ref
 * processing_function passed into @ref Process() acquires an internal mutex on
 * its internal resources prior to calling release on the processing semaphore.
 * Failure to acquire an internal mutex will invalidate order guarantees.
 * Failure to release the semaphore will prevent processing overlap.
 * See @ref Process() for more details.
 *
 * It is the callers responsibility to ensure that at most two @ref Process()
 * calls are ever active at a time, otherwise correct ordering on the later @ref
 * Process() calls cannot be guaranteed.
 */
class OverlappingTaskProcessor {
 public:
  OverlappingTaskProcessor(
      aocommon::TaskQueue<std::function<void()>>& task_queue)
      : task_queue_(task_queue),
        processing_count_semaphore_(2),
        processing_order_semaphore_(1),
        lane_process_semaphore_(1){};
  OverlappingTaskProcessor(const OverlappingTaskProcessor&) = delete;
  OverlappingTaskProcessor& operator=(const OverlappingTaskProcessor&) = delete;
  OverlappingTaskProcessor(OverlappingTaskProcessor&&) = delete;
  OverlappingTaskProcessor& operator=(OverlappingTaskProcessor&&) = delete;

  /*
   * Read chunks from @ref data_lane and call @ref processing_function on them
   * with @ref processing_function() called on each chunk in order and at most
   * two chunks at any given time.
   *
   * @ref processing_function should take the chunk of data to process as its
   * first argument, an index as its second argument and a binary semaphore as
   * its third argument.
   * @ref processing_function must acquire its own internal resource mutex and
   * then call `release()` on the semaphore when it is ready to allow overlap.
   * If this is not done then no overlap will occur, if no mutex is acquired
   * then ordering cannot be guaranteed.
   */
  template <typename DataType>
  void Process(aocommon::Lane<DataType>& data_lane,
               std::function<void(DataType&&, size_t, std::binary_semaphore&)>&&
                   processing_function,
               const std::string& log_tag = "") {
    // Acquired before queuing each task and released after processing it.
    std::counting_semaphore<2> currently_processing_semaphore(2);

    ProcessAllChunksForLane(data_lane, std::move(processing_function),
                            currently_processing_semaphore, log_tag);

    // Acquire twice to ensure that all processing tasks from this call
    // to Process() have completed. Note that tasks from a previous call to
    // Process() can still be queued and/or processing.
    for (size_t i = 0; i < 2; ++i) {
      currently_processing_semaphore.acquire();
    }
    if (!log_tag.empty()) {
      aocommon::Logger::Debug << "All " << log_tag << " chunks processed.\n";
    }
  }

 private:
  template <typename DataType>
  void ProcessAllChunksForLane(
      aocommon::Lane<DataType>& data_lane,
      std::function<void(DataType&&, size_t, std::binary_semaphore&)>&&
          processing_function,
      std::counting_semaphore<2>& currently_processing_semaphore,
      const std::string& log_tag) {
    // Only allow one lane to queue tasks at a time.
    lane_process_semaphore_.acquire();
    DataType chunk_data;
    size_t chunk_index = 0;
    while (data_lane.read(chunk_data)) {
      if (!log_tag.empty()) {
        aocommon::Logger::Debug << "Queue " << log_tag << " chunk "
                                << chunk_index << ".\n";
      }
      QueueChunk(chunk_index, std::move(chunk_data),
                 std::move(processing_function),
                 currently_processing_semaphore);
      ++chunk_index;
    }
    // Allow a second lane to start queuing tasks.
    // Note that while we are done queuing tasks from this lane they have not
    // necessarily all completed processing yet. This is desired behaviour.
    lane_process_semaphore_.release();
    if (!log_tag.empty()) {
      aocommon::Logger::Debug << "All " << log_tag << " chunks queued.\n";
    }
  }

  template <typename DataType>
  void QueueChunk(
      size_t chunk_index, DataType&& chunk_data,
      std::function<void(DataType&&, size_t, std::binary_semaphore&)>&&
          processing_function,
      std::counting_semaphore<2>& currently_processing_semaphore) {
    // Acquire before queuing every task and release after processing every
    // task. This prevents more than 2 tasks from running at a time. This
    // applies also across multiple concurrent Process() calls.
    processing_count_semaphore_.acquire();
    currently_processing_semaphore.acquire();
    // Acquire before queuing a task and release after it starts processing.
    // As this is a binary semaphore it enforces ordering between tasks.
    processing_order_semaphore_.acquire();
    task_queue_.Emplace(
        [this, &processing_function, &currently_processing_semaphore,
         chunk_data = std::move(chunk_data), chunk_index]() mutable {
          // processing function must acquire its own internally owned mutex
          // and then call release() on processing_order_semaphore_ at the first
          // point after that where it would be desirable for an overlap to
          // commence.
          processing_function(std::move(chunk_data), chunk_index,
                              processing_order_semaphore_);
          currently_processing_semaphore.release();
          processing_count_semaphore_.release();
        });
  }

  // This task queue is used to process chunks as they become available.
  // The semaphores are to enforce concurrency/ordering guarantees.
  // See comments at top of header for more information.
  aocommon::TaskQueue<std::function<void()>>& task_queue_;
  std::counting_semaphore<2> processing_count_semaphore_;
  std::binary_semaphore processing_order_semaphore_;
  std::binary_semaphore lane_process_semaphore_;
};

}  // namespace aocommon

#endif  // AOCOMMON_OVERLAPPING_TASK_PROCESSOR_H_
