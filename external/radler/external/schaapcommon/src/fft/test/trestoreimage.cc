// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "restoreimage.h"

#include <cmath>

#include <aocommon/image.h>
#include <aocommon/imagecoordinates.h>

#include "convolution.h"

#include <iostream>

using schaapcommon::fft::MakeFftwfPlannerThreadSafe;
using schaapcommon::fft::RestoreImage;

namespace {
constexpr size_t kWidth = 16;
constexpr size_t kHeight = 16;
constexpr double kFluxDensity = 5.0;
constexpr long double kPixelScaleL = M_PI / (180.0 * 60.0);  // 1 amin
constexpr long double kPixelScaleM = kPixelScaleL;           // 1 amin
constexpr long double kBeamAxis = 2 * kPixelScaleL;
constexpr size_t kThreadCount = 2;

struct ImageFixture {
  ImageFixture() : model(kWidth, kHeight, 0.0), restored(kWidth, kHeight, 0.0) {
    BOOST_CHECK_EQUAL(kWidth, kHeight);
    // Assign center pixel value of 0.
    model[kHeight / 2 * kWidth + kWidth / 2] = kFluxDensity;
  }

  aocommon::Image model;
  aocommon::Image restored;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(restore_image)

BOOST_FIXTURE_TEST_CASE(restore_no_beam, ImageFixture) {
  RestoreImage(restored.Data(), model.Data(), kWidth, kHeight, 0.0, 0.0, 0.0,
               kPixelScaleL, kPixelScaleM, kThreadCount);

  for (size_t i = 0; i != restored.Size(); ++i) {
    BOOST_CHECK_SMALL(std::abs(restored[i] - model[i]), 1e-7f);
  }
}

BOOST_FIXTURE_TEST_CASE(restore_circular_beam, ImageFixture) {
  MakeFftwfPlannerThreadSafe();

  RestoreImage(restored.Data(), model.Data(), kWidth, kHeight, kBeamAxis,
               kBeamAxis, 0.0, kPixelScaleL, kPixelScaleM, kThreadCount);

  // Compute stddev using the full width half maximum (FWHM) formula for a
  // Gaussian:
  const long double sigma =
      kBeamAxis / (2.0L * std::sqrt(2.0L * std::log(2.0L)));

  for (size_t i = 0; i != restored.Size(); ++i) {
    const size_t x = i % restored.Width();
    const size_t y = i / restored.Width();
    // Convert to l,m coordinates
    double l;
    double m;
    aocommon::ImageCoordinates::XYToLM<double>(x, y, kPixelScaleL, kPixelScaleM,
                                               kWidth, kHeight, l, m);
    // Computed distance scaled by stddev
    const double dist_squared = (l * l + m * m) / (sigma * sigma);
    const double gaussian = kFluxDensity * std::exp(-0.5 * dist_squared);
    BOOST_CHECK_SMALL(std::abs(restored[i] - gaussian), 1e-7);

    if (x == kWidth / 2 && y == kHeight / 2) {
      BOOST_CHECK_CLOSE(restored[i], kFluxDensity, 1e-8);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(restore_elliptical_beam, ImageFixture) {
  MakeFftwfPlannerThreadSafe();

  const long double beam_axis_major = 1.5 * kBeamAxis;
  const long double beam_axis_minor = kBeamAxis;
  RestoreImage(restored.Data(), model.Data(), kWidth, kHeight, beam_axis_major,
               beam_axis_minor, 0.0, kPixelScaleL, kPixelScaleM, kThreadCount);

  // Compute stddev using the full width half maximum (FWHM) formula for a
  // Gaussian:
  const long double sigma_major =
      beam_axis_major / (2.0L * std::sqrt(2.0L * std::log(2.0L)));
  const long double sigma_minor =
      beam_axis_minor / (2.0L * std::sqrt(2.0L * std::log(2.0L)));

  for (size_t i = 0; i != restored.Size(); ++i) {
    const size_t x = i % restored.Width();
    const size_t y = i / restored.Width();
    // Convert to l,m coordinates
    double l;
    double m;
    aocommon::ImageCoordinates::XYToLM<double>(x, y, kPixelScaleL, kPixelScaleM,
                                               kWidth, kHeight, l, m);
    // Computed distance scaled by stddev
    const double dist_squared = (l / sigma_minor) * (l / sigma_minor) +
                                (m / sigma_major) * (m / sigma_major);
    const double gaussian = kFluxDensity * std::exp(-0.5 * dist_squared);
    BOOST_CHECK_SMALL(std::abs(restored[i] - gaussian), 3e-7);

    if (x == kWidth / 2 && y == kHeight / 2) {
      BOOST_CHECK_CLOSE(restored[i], kFluxDensity, 1e-8);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()