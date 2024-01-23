#include "resampler.h"

#include <complex>

using aocommon::WindowFunction;
namespace schaapcommon {
namespace fft {

Resampler::Resampler(size_t input_width, size_t input_height,
                     size_t output_width, size_t output_height,
                     size_t cpu_count)
    : input_width_(input_width),
      input_height_(input_height),
      output_width_(output_width),
      output_height_(output_height),
      fft_width_(std::max(input_width, output_width)),
      fft_height_(std::max(input_height, output_height)),
      window_function_(WindowFunction::Rectangular),
      tukey_inset_size_(0.0),
      correct_window_(false),
      tasks_(cpu_count) {
  float* input_data = reinterpret_cast<float*>(
      fftwf_malloc(fft_width_ * fft_height_ * sizeof(float)));
  fftwf_complex* fft_data = reinterpret_cast<fftwf_complex*>(
      fftwf_malloc(fft_width_ * fft_height_ * sizeof(fftwf_complex)));
  in_to_f_plan_ = fftwf_plan_dft_r2c_2d(input_height_, input_width_, input_data,
                                        fft_data, FFTW_ESTIMATE);
  f_to_out_plan_ = fftwf_plan_dft_c2r_2d(output_height_, output_width_,
                                         fft_data, input_data, FFTW_ESTIMATE);
  fftwf_free(fft_data);
  fftwf_free(input_data);
}

Resampler::~Resampler() {
  Finish();
  fftwf_destroy_plan(in_to_f_plan_);
  fftwf_destroy_plan(f_to_out_plan_);
}

void Resampler::RunThread() {
  Task task;
  while (tasks_.read(task)) {
    RunSingle(task, false);
  }
}

void Resampler::RunSingle(const Task& task, bool skip_window) const {
  float* end_ptr = task.input + input_width_ * input_height_;
  for (float* i = task.input; i != end_ptr; ++i) {
    if (!std::isfinite(*i)) *i = 0.0;
  }

  if (window_function_ != WindowFunction::Rectangular && !skip_window) {
    ApplyWindow(task.input);
  }

  const size_t fft_in_width = input_width_ / 2 + 1;
  std::complex<float>* fft_data = reinterpret_cast<std::complex<float>*>(
      fftwf_malloc(fft_in_width * input_height_ * sizeof(std::complex<float>)));
  fftwf_execute_dft_r2c(in_to_f_plan_, task.input,
                        reinterpret_cast<fftwf_complex*>(fft_data));

  const size_t fft_out_width = output_width_ / 2 + 1;
  // TODO this can be done without allocating more mem!
  std::complex<float>* new_fft_data =
      reinterpret_cast<std::complex<float>*>(fftwf_malloc(
          fft_out_width * output_height_ * sizeof(std::complex<float>)));
  std::uninitialized_fill_n(new_fft_data, fft_out_width * output_height_,
                            std::complex<float>(0));

  const size_t old_mid_x = input_width_ / 2;
  const size_t new_mid_x = output_width_ / 2;

  const size_t min_width = std::min(input_width_, output_width_);
  const size_t min_height = std::min(input_height_, output_height_);

  const size_t min_mid_x = min_width / 2;
  const size_t min_mid_y = min_height / 2;

  const float factor = 1.0 / (min_width * min_height);

  for (size_t y = 0; y != min_height; ++y) {
    size_t old_y = y - min_mid_y + input_height_;
    size_t newY = y - min_mid_y + output_height_;
    if (old_y >= input_height_) old_y -= input_height_;
    if (newY >= output_height_) newY -= output_height_;

    // The last dimension is stored half
    for (size_t x = 0; x != min_mid_x; ++x) {
      const size_t old_x = x;
      const size_t new_x = x;
      const size_t old_index = old_x + old_y * (old_mid_x + 1);
      const size_t new_index = new_x + newY * (new_mid_x + 1);

      new_fft_data[new_index] = fft_data[old_index] * factor;

      // if((x == 0 && newY == 0) || (x==min_mid_x-1 && y==min_height-1))
      //	std::cout << new_fft_data[new_index] << " (" << old_x << " , "
      //<<
      // old_y << ") - (" << new_x << " , " << newY << ")\n";
    }
    if (input_width_ >= output_width_) {
      const size_t old_index = input_width_ / 2 + old_y * (old_mid_x + 1);
      const size_t new_index = output_width_ / 2 + newY * (new_mid_x + 1);
      new_fft_data[new_index] = fft_data[old_index] * factor;
    }
  }

  fftwf_free(fft_data);

  fftwf_execute_dft_c2r(f_to_out_plan_,
                        reinterpret_cast<fftwf_complex*>(new_fft_data),
                        task.output);

  fftwf_free(new_fft_data);

  if (correct_window_ && window_function_ != WindowFunction::Rectangular &&
      !skip_window) {
    UnapplyWindow(task.output);
  }
}

void Resampler::SingleFT(const float* input, float* real_output,
                         float* imaginary_output) {
  aocommon::UVector<float> data(input_width_ * input_height_);
  const size_t half_width = input_width_ / 2;
  const size_t half_height = input_height_ / 2;
  for (size_t y = 0; y != input_height_; ++y) {
    size_t y_in = y + half_height;
    if (y_in >= input_height_) y_in -= input_height_;
    float* row_out_ptr = &data[y * input_width_];
    const float* row_in_ptr = &input[y_in * input_width_];
    for (size_t x = 0; x != input_width_; ++x) {
      size_t x_in = x + half_width;
      if (x_in >= input_width_) x_in -= input_width_;
      if (std::isfinite(row_in_ptr[x_in])) {
        row_out_ptr[x] = row_in_ptr[x_in];
      } else {
        row_out_ptr[x] = 0.0;
      }
    }
  }

  const size_t fft_in_width = input_width_ / 2 + 1;
  std::complex<float>* fft_data = reinterpret_cast<std::complex<float>*>(
      fftwf_malloc(fft_in_width * input_height_ * sizeof(std::complex<float>)));

  fftwf_execute_dft_r2c(in_to_f_plan_, data.data(),
                        reinterpret_cast<fftwf_complex*>(fft_data));

  const size_t mid_x = input_width_ / 2;
  const size_t mid_y = input_height_ / 2;

  const float factor = 1.0 / sqrt(input_width_ * input_height_);

  for (size_t y = 0; y != input_height_; ++y) {
    size_t old_y = y + mid_y;
    if (old_y >= input_height_) old_y -= input_height_;

    // The last dimension is stored half
    for (size_t x = 0; x != mid_x + 1; ++x) {
      const size_t old_index = x + old_y * (mid_x + 1);
      const size_t newIndex1 = mid_x - x + y * input_width_;

      const std::complex<float>& val = fft_data[old_index] * factor;

      real_output[newIndex1] = val.real();
      imaginary_output[newIndex1] = val.imag();
      if (x != mid_x) {
        size_t yTo = input_height_ - y;
        if (yTo == input_height_) yTo = 0;
        size_t newIndex2 = mid_x + x + yTo * input_width_;
        real_output[newIndex2] = val.real();
        imaginary_output[newIndex2] = -val.imag();
      }
    }
  }

  fftwf_free(fft_data);
}

void Resampler::MakeWindow(aocommon::UVector<float>& data, size_t width) const {
  if (window_function_ == WindowFunction::Tukey) {
    MakeTukeyWindow(data, width);
  } else {
    data.resize(width);
    for (size_t x = 0; x != width; ++x) {
      data[x] = WindowFunction::Evaluate(window_function_, width, x) + 1e-5;
    }
  }
}

void Resampler::MakeTukeyWindow(aocommon::UVector<float>& data,
                                size_t width) const {
  // Make a Tukey window, which consists of
  // left: a cosine going from 0 to 1
  // mid: all 1
  // right: a cosine going from 1 to 0
  data.resize(width);
  for (size_t x = 0; x != width; ++x) {
    // left part of Tukey window
    const double x_sh = (0.5 + x) * 2;
    if (x_sh < width - tukey_inset_size_) {
      const double pos = x_sh / (width - tukey_inset_size_);
      data[x] = (std::cos((pos + 1.0) * M_PI) + 1.0) * 0.5;
    } else if (x_sh < width + tukey_inset_size_) {
      data[x] = 1.0;
    } else {
      const double pos =
          (x_sh - (width + tukey_inset_size_)) / (width - tukey_inset_size_);
      data[x] = (std::cos(pos * M_PI) + 1.0) * 0.5;
    }
  }
}

void Resampler::ApplyWindow(float* data) const {
  if (window_row_in_.empty()) {
    MakeWindow(window_row_in_, input_width_);
    MakeWindow(window_col_in_, input_height_);
    if (correct_window_) {
      aocommon::UVector<float> window_image_in(input_width_ * input_height_);
      float* in_ptr = window_image_in.data();
      for (size_t y = 0; y != input_height_; ++y) {
        for (size_t x = 0; x != input_width_; ++x) {
          *in_ptr = window_row_in_[x] * window_col_in_[y];
          ++in_ptr;
        }
      }

      window_out_.resize(output_width_ * output_height_);
      Task task;
      task.input = window_image_in.data();
      task.output = window_out_.data();
      RunSingle(task, true);
    }
  }
  for (size_t y = 0; y != input_height_; ++y) {
    for (size_t x = 0; x != input_width_; ++x) {
      *data *= window_row_in_[x] * window_col_in_[y];
      ++data;
    }
  }
}

void Resampler::UnapplyWindow(float* data) const {
  size_t n = output_width_ * output_height_;
  for (size_t i = 0; i != n; ++i) {
    data[i] /= window_out_[i];
  }
}
}  // namespace fft
}  // namespace schaapcommon
