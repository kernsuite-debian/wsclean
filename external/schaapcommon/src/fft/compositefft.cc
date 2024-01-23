// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "compositefft.h"

#include <algorithm>
#include <cstring>

namespace {

// Partially unroll rows/columns with a factor of kUnroll
constexpr size_t kUnroll = 4;

// With kUnroll > 1, the temporary buffers need to be aligned
// for FFTW to work correctly.
constexpr size_t kAlignment = 64;

size_t RoundUp(size_t a, size_t b) { return ((a + b) / b) * b; }

}  // namespace

namespace schaapcommon {
namespace fft {
void FftR2CComposite(fftwf_plan plan_r2c, fftwf_plan plan_c2c,
                     size_t image_height, size_t image_width, const float* in,
                     fftwf_complex* out, aocommon::StaticFor<size_t>& loop) {
  const size_t complex_width = image_width / 2 + 1;
  const size_t complex_size = image_height * complex_width;

  fftwf_complex* temp1 = fftwf_alloc_complex(complex_size);

  loop.Run(0, image_height, [&](size_t y_start, size_t y_end) {
    fftwf_complex* temp2 = fftwf_alloc_complex(complex_width);
    float* temp2_ptr = reinterpret_cast<float*>(temp2);
    for (size_t y = y_start; y < y_end; y++) {
      float* temp1_ptr = reinterpret_cast<float*>(&temp1[y * complex_width]);
      std::copy_n(&in[y * image_width], image_width, temp2_ptr);
      fftwf_execute_dft_r2c(plan_r2c, temp2_ptr, temp2);
      std::copy_n(temp2_ptr, 2 * complex_width, temp1_ptr);
    }
    fftwf_free(temp2);
  });

  loop.Run(0, complex_width, [&](size_t x_start, size_t x_end) {
    // Partially kUnroll over columns
    size_t padded_height = RoundUp(image_height, kAlignment);
    fftwf_complex* temp2 = fftwf_alloc_complex(kUnroll * padded_height);

    for (size_t x = x_start; x < x_end; x += kUnroll) {
      // Copy input
      for (size_t y = 0; y < image_height; y++) {
        for (size_t i = 0; i < kUnroll; i++) {
          if ((x + i) < x_end) {
            float* temp1_ptr =
                reinterpret_cast<float*>(&temp1[y * complex_width + x + i]);
            float* temp2_ptr =
                reinterpret_cast<float*>(&temp2[i * padded_height + y]);
            std::copy_n(temp1_ptr, 2, temp2_ptr);
          }
        }
      }

      // Perform 1D FFT over columns
      for (size_t i = 0; i < kUnroll; i++) {
        fftwf_complex* temp2_ptr = &temp2[i * padded_height];
        fftwf_execute_dft(plan_c2c, temp2_ptr, temp2_ptr);
      }

      // Transpose output
      for (size_t y = 0; y < image_height; y++) {
        for (size_t i = 0; i < kUnroll; i++) {
          if ((x + i) < x_end) {
            float* temp2_ptr =
                reinterpret_cast<float*>(&temp2[i * padded_height + y]);
            float* out_ptr =
                reinterpret_cast<float*>(&out[y * complex_width + x + i]);
            std::copy_n(temp2_ptr, 2, out_ptr);
          }
        }
      }
    }

    fftwf_free(temp2);
  });

  fftwf_free(temp1);
}

void FftC2RComposite(fftwf_plan plan_c2c, fftwf_plan plan_c2r,
                     size_t image_height, size_t image_width,
                     const fftwf_complex* in, float* out,
                     aocommon::StaticFor<size_t>& loop) {
  const size_t complex_width = image_width / 2 + 1;

  size_t padded_height = RoundUp(image_height, kAlignment);
  size_t padded_size = padded_height * complex_width;
  fftwf_complex* temp1 = fftwf_alloc_complex(padded_size);

  loop.Run(0, complex_width, [&](size_t x_start, size_t x_end) {
    for (size_t x = x_start; x < x_end; x += kUnroll) {
      // Transpose input
      for (size_t y = 0; y < image_height; y++) {
        for (size_t i = 0; i < kUnroll; i++) {
          if ((x + i) < x_end) {
            const float* in_ptr =
                reinterpret_cast<const float*>(&in[y * complex_width + x + i]);
            float* temp1_ptr =
                reinterpret_cast<float*>(&temp1[(x + i) * padded_height + y]);
            std::copy_n(in_ptr, 2, temp1_ptr);
          }
        }
      }

      // Perform 1D C2C FFT over columns
      for (size_t i = 0; i < kUnroll; i++) {
        if ((x + i) < x_end) {
          fftwf_complex* temp1_ptr = &temp1[(x + i) * padded_height];
          fftwf_execute_dft(plan_c2c, temp1_ptr, temp1_ptr);
        }
      }
    }
  });

  loop.Run(0, image_height, [&](size_t y_start, size_t y_end) {
    size_t paddedWidth = RoundUp(complex_width, kAlignment);
    fftwf_complex* temp2 = fftwf_alloc_complex(kUnroll * paddedWidth);

    for (size_t y = y_start; y < y_end; y += kUnroll) {
      // Transpose input
      for (size_t x = 0; x < complex_width; x++) {
        for (size_t i = 0; i < kUnroll; i++) {
          if ((y + i) < y_end) {
            float* temp1_ptr =
                reinterpret_cast<float*>(&temp1[x * padded_height + y + i]);
            float* temp2_ptr =
                reinterpret_cast<float*>(&temp2[i * paddedWidth + x]);
            std::copy_n(temp1_ptr, 2, temp2_ptr);
          }
        }
      }

      // Perform 1D C2R FFT over rows
      for (size_t i = 0; i < kUnroll; i++) {
        if ((y + i) < y_end) {
          fftwf_complex* temp2_ptr = &temp2[i * paddedWidth];
          fftwf_execute_dft_c2r(plan_c2r, temp2_ptr,
                                reinterpret_cast<float*>(temp2_ptr));
        }
      }

      // Copy output
      for (size_t i = 0; i < kUnroll; i++) {
        if ((y + i) < y_end) {
          float* temp2_ptr = reinterpret_cast<float*>(&temp2[i * paddedWidth]);
          std::copy_n(temp2_ptr, image_width, &out[(y + i) * image_width]);
        }
      }
    }

    fftwf_free(temp2);
  });

  fftwf_free(temp1);
}
}  // namespace fft
}  // namespace schaapcommon