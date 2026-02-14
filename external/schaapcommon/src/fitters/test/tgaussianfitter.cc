// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "gaussianfitter.h"

#include <aocommon/image.h>
#include <aocommon/threadpool.h>

#include "../../../include/schaapcommon/math/convolution.h"
#include "../../../include/schaapcommon/math/restoreimage.h"

namespace {
constexpr size_t kThreadCount = 2;
constexpr size_t kWidth = 64;
constexpr size_t kHeight = 64;
constexpr long double kPixelSize = 1 /*amin*/ * (M_PI / 180.0 / 60.0);
}  // namespace

using schaapcommon::math::Ellipse;

BOOST_AUTO_TEST_SUITE(gaussian_fitter)

BOOST_AUTO_TEST_CASE(fit) {
  aocommon::ThreadPool::GetInstance().SetNThreads(kThreadCount);
  for (size_t beam_position_angle_index = 0; beam_position_angle_index != 10;
       ++beam_position_angle_index) {
    const size_t width = 512, height = 512;
    aocommon::Image model(width, height, 0.0);
    aocommon::Image restored(width, height, 0.0);
    model[((height / 2) * width) + (width / 2)] = 1.0;
    const long double kPixelSize = 1.0L /*amin*/ * (M_PI / 180.0 / 60.0);
    const long double beam_major = 20.0L * kPixelSize;
    const long double beam_minor = 5.0L * kPixelSize;
    const long double beam_position_angle =
        beam_position_angle_index * M_PI / 10.0;

    schaapcommon::math::MakeFftwfPlannerThreadSafe();
    schaapcommon::math::RestoreImage(
        restored.Data(), model.Data(), width, height, beam_major, beam_minor,
        beam_position_angle, kPixelSize, kPixelSize);

    const Ellipse ellipse = schaapcommon::fitters::Fit2DGaussianCentred(
        restored.Data(), width, height, 5.0, 10.0, false);

    BOOST_CHECK_CLOSE_FRACTION(ellipse.major, 20.0, 1.0e-3);
    BOOST_CHECK_CLOSE_FRACTION(ellipse.minor, 5.0, 1.0e-3);
    const double fit_position_angle =
        std::fmod((ellipse.position_angle + 2.0 * M_PI), M_PI);
    BOOST_CHECK_CLOSE_FRACTION(fit_position_angle, beam_position_angle, 1.0e-3);
  }
}

BOOST_AUTO_TEST_CASE(fit_full) {
  aocommon::ThreadPool::GetInstance().SetNThreads(kThreadCount);
  const size_t width = 512;
  const size_t height = 512;
  const size_t source_x = 200;
  const size_t source_y = 300;
  for (size_t beam_position_angle_index = 0; beam_position_angle_index != 10;
       ++beam_position_angle_index) {
    aocommon::Image model(width, height, 0.0);
    aocommon::Image restored(width, height, 0.0);
    model[(source_y * width) + source_x] = 2.0;
    const long double kPixelSize = 1.0L /*amin*/ * (M_PI / 180.0 / 60.0);
    long double beam_major = 20.0L * kPixelSize;
    long double beam_minor = 5.0L * kPixelSize;
    long double beam_position_angle = beam_position_angle_index * M_PI / 10.0;

    schaapcommon::math::MakeFftwfPlannerThreadSafe();
    schaapcommon::math::RestoreImage(
        restored.Data(), model.Data(), width, height, beam_major, beam_minor,
        beam_position_angle, kPixelSize, kPixelSize);

    Ellipse fit_ellipse;
    // The following four parameters are used as initial values.
    fit_ellipse.major = 10.0;
    double fit_amplitude = 1.0;
    double fit_x = source_x + 3;
    double fit_y = source_y - 3;
    schaapcommon::fitters::Fit2DGaussianFull(
        restored.Data(), width, height, fit_amplitude, fit_x, fit_y,
        fit_ellipse.major, fit_ellipse.minor, fit_ellipse.position_angle,
        nullptr);

    BOOST_CHECK_CLOSE_FRACTION(fit_amplitude, 2.0, 1.0e-3);
    BOOST_CHECK_CLOSE_FRACTION(fit_x, source_x, 1.0e-3);
    BOOST_CHECK_CLOSE_FRACTION(fit_y, source_y, 1.0e-3);
    BOOST_CHECK_CLOSE_FRACTION(fit_ellipse.major, 20.0, 1.0e-3);
    BOOST_CHECK_CLOSE_FRACTION(fit_ellipse.minor, 5.0, 1.0e-3);
    const double fit_position_angle =
        std::fmod((fit_ellipse.position_angle + 2.0 * M_PI), M_PI);
    BOOST_CHECK_CLOSE_FRACTION(fit_position_angle, beam_position_angle, 1.0e-3);
  }
}

BOOST_AUTO_TEST_CASE(insufficient_data_fit) {
  const size_t width = 1, height = 1;
  aocommon::Image model(width, height, 0.0);
  aocommon::Image restored(width, height, 0.0);
  model[((height / 2) * width) + (width / 2)] = 1.0;
  BOOST_CHECK_NO_THROW(schaapcommon::fitters::Fit2DGaussianCentred(
      restored.Data(), width, height, 1.0, 10.0, false));
}

BOOST_AUTO_TEST_CASE(fit_circular) {
  aocommon::Image model(kWidth, kHeight, 0.0);
  aocommon::Image restored(kWidth, kHeight, 0.0);

  model[((kHeight / 2) * kWidth) + (kWidth / 2)] = 1.0;

  const long double beam_major = 4.0L * kPixelSize;
  const long double beam_minor = 4.0L * kPixelSize;
  const long double beam_position_angle = 0.0;
  const long double estimated_beam_pixel = 1.0;  // this is on purpose way off
  schaapcommon::math::MakeFftwfPlannerThreadSafe();
  schaapcommon::math::RestoreImage(restored.Data(), model.Data(), kWidth,
                                   kHeight, beam_major, beam_minor,
                                   beam_position_angle, kPixelSize, kPixelSize);

  const Ellipse result = schaapcommon::fitters::Fit2DGaussianCentred(
      restored.Data(), restored.Width(), restored.Height(),
      estimated_beam_pixel, 10.0, false);

  BOOST_CHECK_CLOSE_FRACTION(result.major, 4.0, 1.0e-4);
  BOOST_CHECK_CLOSE_FRACTION(result.minor, 4.0, 1.0e-4);
  BOOST_CHECK_SMALL(std::abs(result.position_angle -
                             static_cast<double>(beam_position_angle)),
                    1.0e-4);
}

BOOST_AUTO_TEST_CASE(little_data_circular_fit) {
  const size_t width = 1, height = 1;
  aocommon::Image model(width, height, 0.0);
  aocommon::Image restored(width, height, 0.0);
  model[((height / 2) * width) + (width / 2)] = 1.0;
  double fit = 0.0;
  BOOST_CHECK_NO_THROW(schaapcommon::fitters::Fit2DCircularGaussianCentred(
      restored.Data(), width, height, fit));
}

BOOST_AUTO_TEST_CASE(fit_small_beam) {
  aocommon::ThreadPool::GetInstance().SetNThreads(kThreadCount);
  aocommon::Image model(kWidth, kHeight, 0.0);
  aocommon::Image restored(kWidth, kHeight, 0.0);

  model[((kHeight / 2) * kWidth) + (kWidth / 2)] = 1.0;

  const long double beam_major = 4.0L * kPixelSize;
  const long double beam_minor = 0.5L * kPixelSize;
  const long double beam_position_angle = 0.0;
  const long double estimated_beam_pixel = 1.0;  // this is on purpose way off

  schaapcommon::math::MakeFftwfPlannerThreadSafe();
  schaapcommon::math::RestoreImage(restored.Data(), model.Data(), kWidth,
                                   kHeight, beam_major, beam_minor,
                                   beam_position_angle, kPixelSize, kPixelSize);

  Ellipse ellipse = schaapcommon::fitters::Fit2DGaussianCentred(
      restored.Data(), restored.Width(), restored.Height(),
      estimated_beam_pixel, 10.0, false);

  BOOST_CHECK_CLOSE_FRACTION(ellipse.major, 4.0, 1.0e-4);
  BOOST_CHECK_CLOSE_FRACTION(ellipse.minor, 0.5, 1.0e-4);
  BOOST_CHECK_SMALL(std::abs(ellipse.position_angle -
                             static_cast<double>(beam_position_angle)),
                    1.0e-4);
}

BOOST_AUTO_TEST_CASE(deconvolve) {
  using schaapcommon::fitters::DeconvolveGaussian;
  Ellipse result;

  // Evaluate all four kwadrants.
  for (double phi = 0.0; phi < 2.0 * M_PI; phi += M_PI / 4.0) {
    // Simple case: equal PAs and zero minor axis for second ellipse, basically
    // a 1D deconvolution.
    const double input_phi = 0.5 + phi;
    result = DeconvolveGaussian(Ellipse{2.0, 1.0, input_phi},
                                Ellipse{0.2, 0.0, input_phi});
    BOOST_CHECK_CLOSE_FRACTION(result.major, std::sqrt(2.0 * 2.0 - 0.2 * 0.2),
                               1e-5);
    BOOST_CHECK_CLOSE_FRACTION(result.minor, 1.0, 1e-5);
    double expected = input_phi;
    while (expected > M_PI * 0.5) expected -= M_PI;
    BOOST_CHECK_CLOSE_FRACTION(result.position_angle, expected, 1e-5);

    // Same as above, but rotate second ellipse by 180 deg (which should not
    // change anything).
    result = DeconvolveGaussian(Ellipse{2.0, 1.0, input_phi},
                                Ellipse{0.2, 0.0, input_phi + M_PI});
    BOOST_CHECK_CLOSE_FRACTION(result.major, std::sqrt(2.0 * 2.0 - 0.2 * 0.2),
                               1e-5);
    BOOST_CHECK_CLOSE_FRACTION(result.minor, 1.0, 1e-5);
    BOOST_CHECK_CLOSE_FRACTION(result.position_angle, expected, 1e-5);
  }

  // Deconvolve with zero size ellipse: should return first ellipse.
  result = DeconvolveGaussian(Ellipse{2.0, 1.0, 0.5}, Ellipse{0.0, 0.0, 0.1});
  BOOST_CHECK_CLOSE_FRACTION(result.major, std::sqrt(2.0 * 2.0 - 0.0 * 0.0),
                             1e-5);
  BOOST_CHECK_CLOSE_FRACTION(result.minor, 1.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(result.position_angle, 0.5, 1e-5);

  // Equal sizes: should return zero-sized ellipse.
  result = DeconvolveGaussian(Ellipse{2.0, 1.0, 0.5}, Ellipse{2.0, 1.0, 0.5});
  BOOST_CHECK_LT(std::abs(result.major), 1e-5);
  BOOST_CHECK_LT(std::abs(result.minor), 1e-5);
  // position angle is somewhat undefined, but should still be finite...
  BOOST_CHECK(std::isfinite(result.position_angle));
  // ...and within bounds.
  BOOST_CHECK_LE(result.position_angle, 2.0 * M_PI);
  BOOST_CHECK_GE(result.position_angle, -2.0 * M_PI);

  // 90 degree different position angle.
  result =
      DeconvolveGaussian(Ellipse{5.0, 2.0, 0.0}, Ellipse{1.0, 0.0, M_PI * 0.5});
  BOOST_CHECK_CLOSE_FRACTION(result.major, 5.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(result.minor, std::sqrt(2.0 * 2.0 - 1.0 * 1.0),
                             1e-5);
  BOOST_CHECK_LT(std::abs(result.position_angle), 1e-5);

  // Circular deconvolution.
  result = DeconvolveGaussian(Ellipse{5.0, 5.0, -0.3}, Ellipse{3.0, 3.0, 0.7});
  BOOST_CHECK_CLOSE_FRACTION(result.major, 4.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(result.minor, 4.0, 1e-5);
  BOOST_CHECK(std::isfinite(result.position_angle));
  BOOST_CHECK_LE(result.position_angle, 2.0 * M_PI);
  BOOST_CHECK_GE(result.position_angle, -2.0 * M_PI);

  // A complex case (numbers were calculated using the code, assuming it is
  // correct).
  result = DeconvolveGaussian(Ellipse{10.0, 8.0, 0.3}, Ellipse{7.0, 5.0, 1.4});
  BOOST_CHECK_CLOSE_FRACTION(result.major, 8.477876113222349, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(result.minor, 4.2574190079030156, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(result.position_angle, 0.11532393547063115, 1e-5);

  // Overflow situation.
  result = DeconvolveGaussian(Ellipse{3.0, 3.0, 0.0}, Ellipse{4.0, 0.0, 0.0});
  BOOST_CHECK(!std::isfinite(result.major));
  BOOST_CHECK(!std::isfinite(result.major));
  BOOST_CHECK(!std::isfinite(result.major));
}

BOOST_AUTO_TEST_SUITE_END()
