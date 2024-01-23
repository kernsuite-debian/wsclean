// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "resampler.h"

#include <cmath>

#include <aocommon/image.h>

using schaapcommon::fft::Resampler;

namespace {
constexpr size_t kInputWidth = 16;
constexpr size_t kInputHeight = 16;

constexpr size_t kFactorLow = 2;
constexpr size_t kOutputWidthLow = kInputWidth / kFactorLow;
constexpr size_t kOutputHeightLow = kInputHeight / kFactorLow;

constexpr size_t kFactorHigh = 4;
constexpr size_t kOutputWidthHigh = kInputWidth * kFactorHigh;
constexpr size_t kOutputHeightHigh = kInputWidth * kFactorHigh;

constexpr size_t kThreadCount = 1;

/**
 * @brief Construct a signal on the image. Tune wave number \c k such that
 * signal is periodic.
 */
void ConstructSignal(aocommon::Image& image, float k) {
  for (size_t i = 0; i < image.Size(); ++i) {
    const size_t col = i % image.Width();
    const size_t row = i / image.Width();
    image[i] = std::sin(k * col) + std::cos(k * row);
  }
}

void CheckResampledImage(const aocommon::Image& image_low_resolution,
                         const aocommon::Image& image_high_resolution,
                         bool scale_results) {
  BOOST_CHECK_EQUAL(image_low_resolution.Width(),
                    image_low_resolution.Height());
  BOOST_CHECK_EQUAL(image_high_resolution.Width(),
                    image_high_resolution.Height());
  BOOST_CHECK(image_high_resolution.Width() >= image_low_resolution.Width());

  const size_t scale_factor =
      image_high_resolution.Width() / image_low_resolution.Width();

  size_t dummy_index = 0;
  for (size_t i = 0; i != image_high_resolution.Size(); ++i) {
    const size_t row = i / image_high_resolution.Width();
    const size_t col = i % image_high_resolution.Width();

    if ((row % scale_factor == 0) && (col % scale_factor == 0)) {
      if (scale_results) {
        BOOST_CHECK_SMALL(
            image_high_resolution[i] * scale_factor * scale_factor -
                image_low_resolution[dummy_index],
            static_cast<float>(1e-5));
      } else {
        BOOST_CHECK_SMALL(
            image_high_resolution[i] - image_low_resolution[dummy_index],
            static_cast<float>(1e-5));
      }
      ++dummy_index;
    }
  }
  BOOST_CHECK_EQUAL(dummy_index, image_low_resolution.Size());
}

}  // namespace

BOOST_AUTO_TEST_SUITE(fft_resampler)

BOOST_AUTO_TEST_CASE(downsample_regular_window) {
  BOOST_CHECK(kFactorLow >= 1);

  Resampler resampler(kInputWidth, kInputHeight, kOutputWidthLow,
                      kOutputHeightLow, kThreadCount);

  aocommon::Image input_image(kInputWidth, kInputHeight);
  aocommon::Image output_image(kOutputWidthLow, kOutputHeightLow, 0);
  ConstructSignal(input_image, 0.125 * M_PI);

  resampler.Resample(input_image.Data(), output_image.Data());

  const bool scale_output = true;
  CheckResampledImage(output_image, input_image, scale_output);
}

BOOST_AUTO_TEST_CASE(downsample_tukey_window) {
  BOOST_CHECK(kFactorLow >= 1);

  Resampler resampler(kInputWidth, kInputHeight, kOutputWidthLow,
                      kOutputHeightLow, kThreadCount);

  aocommon::Image input_image(kInputWidth, kInputHeight);
  aocommon::Image output_image(kOutputWidthLow, kOutputHeightLow, 0);
  ConstructSignal(input_image, 0.125 * M_PI);

  resampler.SetTukeyWindow(kInputWidth, kInputWidth);
  resampler.Resample(input_image.Data(), output_image.Data());

  const bool scale_output = false;
  CheckResampledImage(output_image, input_image, scale_output);
}

BOOST_AUTO_TEST_CASE(upsample_tukey_window) {
  Resampler resampler(kInputWidth, kInputHeight, kOutputWidthHigh,
                      kOutputHeightHigh, kThreadCount);

  aocommon::Image input_image(kInputWidth, kInputHeight);
  aocommon::Image output_image(kOutputWidthHigh, kOutputHeightHigh, 0);
  ConstructSignal(input_image, 0.125 * M_PI);

  resampler.SetTukeyWindow(kInputWidth, kInputWidth);
  resampler.Resample(input_image.Data(), output_image.Data());

  const bool scale_output = false;
  CheckResampledImage(input_image, output_image, scale_output);
}

BOOST_AUTO_TEST_CASE(upsample_regular_window) {
  Resampler resampler(kInputWidth, kInputHeight, kOutputWidthHigh,
                      kOutputHeightHigh, kThreadCount);

  aocommon::Image input_image(kInputWidth, kInputHeight);
  aocommon::Image output_image(kOutputWidthHigh, kOutputHeightHigh, 0);
  ConstructSignal(input_image, 0.125 * M_PI);

  resampler.Resample(input_image.Data(), output_image.Data());

  const bool scale_output = false;
  CheckResampledImage(input_image, output_image, scale_output);
}

BOOST_AUTO_TEST_SUITE_END()