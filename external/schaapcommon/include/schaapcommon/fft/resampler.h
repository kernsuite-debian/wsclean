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

namespace schaapcommon {
namespace fft {

class Resampler {
 private:
  struct Task {
    float *input, *output;
  };

 public:
  Resampler(size_t input_width, size_t input_height, size_t output_width,
            size_t output_height, size_t cpu_count);

  ~Resampler();

  void AddTask(float* input, float* output) {
    Task task;
    task.input = input;
    task.output = output;
    tasks_.write(task);
  }

  void Start() {
    for (size_t i = 0; i != tasks_.capacity(); ++i) {
      threads_.emplace_back(&Resampler::RunThread, this);
    }
  }

  void Finish() {
    tasks_.write_end();
    for (std::thread& t : threads_) t.join();
    threads_.clear();
    tasks_.clear();
  }

  void Resample(float* input, float* output) {
    Task task;
    task.input = input;
    task.output = output;
    RunSingle(task, false);
  }

  void SingleFT(const float* input, float* real_output,
                float* imaginary_output);

  /**
   * Only to be used with SingleFT (it makes resampling thread unsafe!)
   */
  void SetTukeyWindow(double inset_size, bool correct_window) {
    window_function_ = aocommon::WindowFunction::Tukey;
    tukey_inset_size_ = inset_size;
    correct_window_ = correct_window;
    window_row_in_.clear();
    window_col_in_.clear();
    window_out_.clear();
  }

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
}  // namespace fft
}  // namespace schaapcommon

#endif
