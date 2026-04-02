// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_FFT_RESAMPLER_H_
#define SCHAAPCOMMON_FFT_RESAMPLER_H_

#include <aocommon/image.h>
#include <aocommon/lane.h>
#include <aocommon/uvector.h>
#include <aocommon/windowfunction.h>

#include <vector>
#include <thread>

#include <fftw3.h>

namespace schaapcommon::math {

/**
 * Class that can resize an image using sinc interpolation. It does this by
 * calculating the FT of the image, followed by zero-padding, and
 * FT-ing the larger result back to image space. This class has also
 * functionality to perform several resampling operation in parallel.
 *
 * If parallel processing is not necessary, @ref Resample() can be used.
 * If only the real-to-complex FT is necessary with the given input dimensions,
 * @ref SingleFT() can be used.
 *
 * Any NaN values are replaced by zeros before FTing.
 *
 * @TODO The structure of this class is based around older behaviour. It would
 * make the functionality easier to by:
 * i) make stand-alone functions for SingleFT() and Resample();
 * ii) make use of the thread pool stucture in aocommon for parallelisation.
 */
class Resampler {
 private:
  struct Task {
    float *input, *output;
  };

 public:
  /**
   * Construct a resampler that can increase an image from the input
   * dimensions to the output dimensions.
   * @TODO The @p cpu_count parameter should be removed, as the resampler
   * should use the global thread pool.
   */
  Resampler(size_t input_width, size_t input_height, size_t output_width,
            size_t output_height, size_t cpu_count);

  ~Resampler();

  /**
   * Add a task to the queue. Before doing so, @ref Start() must have been
   * called. The @p output parameter can be assumed to be finished only
   * after @ref Finish() has been called and has returned. If the queue
   * is full (all cpus are working), the function will block until
   * one of the tasks is finished and a cpu is available.
   *
   * @param input an array of input_width x input_height values. The
   * input will be destroyed.
   * @param output array of output_width x output_height values in
   * which the resized image will be stored.
   */
  void AddTask(float* input, float* output) {
    Task task;
    task.input = input;
    task.output = output;
    tasks_.write(task);
  }

  /**
   * Claim threads and start processing any tasks that are added. This
   * function is non-blocking (returns immediately). After calling Start(),
   * tasks can be added using @ref AddTask(). Once all tasks are added,
   * @ref Finish() should be called.
   */
  void Start() {
    for (size_t i = 0; i != tasks_.capacity(); ++i) {
      threads_.emplace_back(&Resampler::RunThread, this);
    }
  }

  /**
   * Wait until all added tasks have been finished. @ref Start() should
   * have been called beforehand.
   * After calling @ref Finish(), no more tasks should be added until
   * @ref Start() is called.
   */
  void Finish() {
    tasks_.write_end();
    for (std::thread& t : threads_) t.join();
    threads_.clear();
    tasks_.clear();
  }

  /**
   * Directly perform a single resampling actions. This function does
   * not parallelise and blocks until finished.
   * @param input an array of input_width x input_height values. The
   * input will be destroyed.
   * @param output array of output_width x output_height values in
   * which the resized image will be stored.
   */
  void Resample(float* input, float* output) {
    Task task;
    task.input = input;
    task.output = output;
    RunSingle(task, false);
  }

  /**
   * Directly perform a single real-to-complex Fourier transform.
   * @param input Data, with size of input_width x input_height.
   * @param real_output Real values after FT-ing, of size input_width x
   * input_height.
   * @param imaginary_output Corresponding imaginary values.
   */
  void SingleFT(const float* input, float* real_output,
                float* imaginary_output);

  /**
   * Set the resampler to apply a Tukey window during the Fourier transform.
   * Only to be used with @ref Resample() (it makes resampling thread unsafe!)
   */
  void SetTukeyWindow(double inset_size, bool correct_window) {
    window_function_ = aocommon::WindowFunction::Tukey;
    tukey_inset_size_ = inset_size;
    correct_window_ = correct_window;
    window_row_in_.clear();
    window_col_in_.clear();
    window_out_.clear();
  }

  /**
   * Set the resampler to apply a window function during the Fourier transform.
   * Only to be used with @ref Resample() (it makes resampling thread unsafe!)
   */
  void SetWindowFunction(aocommon::WindowFunction::Type window,
                         bool correct_window) {
    window_function_ = window;
    correct_window_ = correct_window;
    window_row_in_.clear();
    window_col_in_.clear();
    window_out_.clear();
  }

 private:
  void RunThread();
  void RunSingle(const Task& task, bool skip_window) const;
  void ApplyWindow(float* data) const;
  void UnapplyWindow(float* data) const;
  void MakeWindow(aocommon::UVector<float>& data, size_t width) const;
  void MakeTukeyWindow(aocommon::UVector<float>& data, size_t width) const;

  size_t input_width_;
  size_t input_height_;
  size_t output_width_;
  size_t output_height_;
  size_t fft_width_;
  size_t fft_height_;
  aocommon::WindowFunction::Type window_function_;
  double tukey_inset_size_;
  mutable aocommon::UVector<float> window_row_in_;
  mutable aocommon::UVector<float> window_col_in_;
  mutable aocommon::UVector<float> window_out_;
  bool correct_window_;

  fftwf_plan in_to_f_plan_;
  fftwf_plan f_to_out_plan_;

  aocommon::Lane<Task> tasks_;
  std::vector<std::thread> threads_;
};

}  // namespace schaapcommon::math

#endif
