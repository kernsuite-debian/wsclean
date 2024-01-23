// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "convolution.h"

#include <aocommon/image.h>
#include <iostream>

namespace {
constexpr size_t kWidth = 4;
constexpr size_t kHeight = 4;
constexpr size_t kThreadCount = 2;
}  // namespace

BOOST_AUTO_TEST_SUITE(fft_convolution)

BOOST_AUTO_TEST_CASE(prepare_kernel) {
  // Values are chosen such that a translation to the origin yields the values
  // in ascending order.
  aocommon::Image kernel_in(kWidth, kHeight,
                            {10.0, 11.0, 8.0, 9.0, 14.0, 15.0, 12.0, 13.0, 2.0,
                             3.0, 0.0, 1.0, 6.0, 7.0, 4.0, 5.0});
  aocommon::Image kernel_out(kWidth, kHeight);
  schaapcommon::fft::PrepareConvolutionKernel(
      kernel_out.Data(), kernel_in.Data(), kWidth, kHeight, kThreadCount);
  for (size_t i = 0; i != kernel_out.Size(); ++i) {
    BOOST_CHECK_CLOSE(kernel_out[i], static_cast<float>(i), 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(prepare_small_kernel) {
  BOOST_CHECK_EQUAL(kWidth, kHeight);
  aocommon::Image kernel_in(kWidth / 2, kHeight / 2, {16.0, 13.0, 4.0, 1.0});
  aocommon::Image kernel_out(kWidth, kHeight, 0.0);
  BOOST_CHECK_THROW(schaapcommon::fft::PrepareSmallConvolutionKernel(
                        kernel_out.Data(), kWidth, kHeight, kernel_in.Data(),
                        kWidth * 2, kThreadCount),
                    std::runtime_error);

  schaapcommon::fft::PrepareSmallConvolutionKernel(kernel_out.Data(), kWidth,
                                                   kHeight, kernel_in.Data(),
                                                   kWidth / 2, kThreadCount);

  // kernel_out should have non-zero corner values only
  for (size_t i = 0; i != kernel_out.Size(); ++i) {
    if (i != 0 && i != (kWidth - 1) && i != (kWidth * (kHeight - 1)) &&
        i != kWidth * kHeight - 1) {
      BOOST_CHECK_CLOSE(kernel_out[i], 0.0f, 1e-4);
    } else {
      BOOST_CHECK_CLOSE(kernel_out[i], static_cast<float>(i + 1), 1e-4);
    }
  }
}

BOOST_AUTO_TEST_CASE(convolve) {
  const float dirac_scale = 0.5;
  aocommon::Image image(kWidth, kHeight);
  for (size_t i = 0; i != image.Size(); ++i) {
    image[i] = i;
  }
  aocommon::Image kernel(kWidth, kHeight, 0.0);
  aocommon::Image image_ref = image;

  kernel[kWidth * kHeight / 2 + kWidth / 2] = dirac_scale * 1.0f;

  schaapcommon::fft::MakeFftwfPlannerThreadSafe();
  BOOST_CHECK_THROW(schaapcommon::fft::ResizeAndConvolve(
                        image.Data(), kWidth, kHeight, kernel.Data(),
                        kWidth * 2, kThreadCount),
                    std::runtime_error);

  schaapcommon::fft::ResizeAndConvolve(image.Data(), kWidth, kHeight,
                                       kernel.Data(), kWidth, kThreadCount);

  // Convolution with (scaled) dirac kernel should return same signal,
  // scaled by a constant factor.
  for (size_t i = 0; i != image.Size(); ++i) {
    BOOST_CHECK_CLOSE(image[i], image_ref[i] * dirac_scale, 1e-4);
  }
}

BOOST_AUTO_TEST_SUITE_END()