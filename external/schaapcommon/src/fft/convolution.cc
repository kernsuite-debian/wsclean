// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "convolution.h"

#include <complex>
#include <stdexcept>

#include <aocommon/uvector.h>
#include <aocommon/staticfor.h>

#include <fftw3.h>

#include "compositefft.h"

#include <iostream>

namespace schaapcommon {
namespace fft {

void MakeFftwfPlannerThreadSafe() { fftwf_make_planner_thread_safe(); }

void ResizeAndConvolve(float* image, size_t image_width, size_t image_height,
                       const float* kernel, size_t kernel_size,
                       size_t thread_count) {
  aocommon::UVector<float> scaled_kernel(image_width * image_height, 0.0);
  PrepareSmallConvolutionKernel(scaled_kernel.data(), image_width, image_height,
                                kernel, kernel_size, thread_count);
  Convolve(image, scaled_kernel.data(), image_width, image_height,
           thread_count);
}

void PrepareSmallConvolutionKernel(float* dest, size_t image_width,
                                   size_t image_height, const float* kernel,
                                   size_t kernel_size, size_t thread_count) {
  if (kernel_size > image_width || kernel_size > image_height) {
    throw std::runtime_error("Kernel size > image dimension");
  }
  aocommon::StaticFor<size_t> loop(thread_count);
  loop.Run(0, kernel_size / 2, [&](size_t y_start, size_t y_end) {
    const float* kernel_iter = &kernel[y_start * kernel_size];
    for (size_t y = y_start; y != y_end; ++y) {
      const size_t y_dest = image_height - kernel_size / 2 + y;
      const size_t x_first = image_width - kernel_size / 2;
      float* dest_iter = &dest[y_dest * image_width + x_first];
      for (size_t x = 0; x != kernel_size / 2; ++x) {
        *dest_iter = *kernel_iter;
        ++kernel_iter;
        ++dest_iter;
      }
      dest_iter = &dest[y_dest * image_width];
      for (size_t x = kernel_size / 2; x != kernel_size; ++x) {
        *dest_iter = *kernel_iter;
        ++kernel_iter;
        ++dest_iter;
      }
    }
  });
  loop.Run(kernel_size / 2, kernel_size, [&](size_t y_start, size_t y_end) {
    const float* kernel_iter = &kernel[y_start * kernel_size];
    for (size_t y = y_start; y != y_end; ++y) {
      size_t x_first = image_width - kernel_size / 2;
      float* dest_iter = &dest[x_first + (y - kernel_size / 2) * image_width];
      for (size_t x = 0; x != kernel_size / 2; ++x) {
        *dest_iter = *kernel_iter;
        ++kernel_iter;
        ++dest_iter;
      }
      dest_iter = &dest[(y - kernel_size / 2) * image_width];
      for (size_t x = kernel_size / 2; x != kernel_size; ++x) {
        *dest_iter = *kernel_iter;
        ++kernel_iter;
        ++dest_iter;
      }
    }
  });
}

void PrepareConvolutionKernel(float* dest, const float* source,
                              size_t image_width, size_t image_height,
                              size_t thread_count) {
  aocommon::StaticFor<size_t> loop(thread_count);
  loop.Run(0, image_height / 2, [&](size_t y_start, size_t y_end) {
    const float* source_iter = &source[y_start * image_width];
    for (size_t y = y_start; y != y_end; ++y) {
      size_t y_dest = image_height - image_height / 2 + y;
      size_t x_first = image_width - image_width / 2;
      float* dest_iter = &dest[y_dest * image_width + x_first];
      for (size_t x = 0; x != image_width / 2; ++x) {
        *dest_iter = *source_iter;
        ++source_iter;
        ++dest_iter;
      }
      dest_iter = &dest[y_dest * image_width];
      for (size_t x = image_width / 2; x != image_width; ++x) {
        *dest_iter = *source_iter;
        ++source_iter;
        ++dest_iter;
      }
    }
  });
  loop.Run(image_height / 2, image_height, [&](size_t y_start, size_t y_end) {
    const float* source_iter = &source[y_start * image_width];
    for (size_t y = y_start; y != y_end; ++y) {
      size_t x_first = image_width - image_width / 2;
      float* dest_iter = &dest[x_first + (y - image_height / 2) * image_width];
      for (size_t x = 0; x != image_width / 2; ++x) {
        *dest_iter = *source_iter;
        ++source_iter;
        ++dest_iter;
      }
      dest_iter = &dest[(y - image_height / 2) * image_width];
      for (size_t x = image_width / 2; x != image_width; ++x) {
        *dest_iter = *source_iter;
        ++source_iter;
        ++dest_iter;
      }
    }
  });
}

void Convolve(float* image, const float* kernel, size_t image_width,
              size_t image_height, size_t thread_count) {
  const size_t image_size = image_width * image_height;
  const size_t complex_width = image_width / 2 + 1;
  const size_t complex_size = complex_width * image_height;
  float* temp_data = fftwf_alloc_real(image_size);
  fftwf_complex* fft_image_data = fftwf_alloc_complex(complex_size);
  fftwf_complex* fft_kernel_data = fftwf_alloc_complex(complex_size);

  fftwf_plan plan_r2c =
      fftwf_plan_dft_r2c_1d(image_width, nullptr, nullptr, FFTW_ESTIMATE);
  fftwf_plan plan_c2c_forward = fftwf_plan_dft_1d(
      image_height, nullptr, nullptr, FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_plan plan_c2c_backward = fftwf_plan_dft_1d(
      image_height, nullptr, nullptr, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwf_plan plan_c2r =
      fftwf_plan_dft_c2r_1d(image_width, nullptr, nullptr, FFTW_ESTIMATE);

  aocommon::StaticFor<size_t> loop(thread_count);

  FftR2CComposite(plan_r2c, plan_c2c_forward, image_height, image_width, image,
                  fft_image_data, loop);

  std::copy_n(kernel, image_size, temp_data);
  FftR2CComposite(plan_r2c, plan_c2c_forward, image_height, image_width,
                  temp_data, fft_kernel_data, loop);

  const float fact = 1.0 / image_size;
  loop.Run(0, image_height, [&](size_t y_start, size_t y_end) {
    for (size_t y = y_start; y != y_end; ++y) {
      for (size_t x = 0; x != complex_width; ++x) {
        const size_t i = y * complex_width + x;
        reinterpret_cast<std::complex<float>*>(fft_image_data)[i] *=
            fact * reinterpret_cast<std::complex<float>*>(fft_kernel_data)[i];
      }
    }
  });

  FftC2RComposite(plan_c2c_backward, plan_c2r, image_height, image_width,
                  fft_image_data, image, loop);

  fftwf_free(fft_image_data);
  fftwf_free(fft_kernel_data);
  fftwf_free(temp_data);

  fftwf_destroy_plan(plan_r2c);
  fftwf_destroy_plan(plan_c2c_forward);
  fftwf_destroy_plan(plan_c2c_backward);
  fftwf_destroy_plan(plan_c2r);
}

}  // namespace fft
}  // namespace schaapcommon