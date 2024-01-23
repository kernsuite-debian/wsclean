// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "gaussianfitter.h"

#include <aocommon/image.h>

#include "../../../include/schaapcommon/fft/convolution.h"
#include "../../../include/schaapcommon/fft/restoreimage.h"

namespace {
constexpr size_t kThreadCount = 2;
constexpr size_t kWidth = 64;
constexpr size_t kHeight = 64;
constexpr long double kPixelSize = 1 /*amin*/ * (M_PI / 180.0 / 60.0);
}  // namespace

BOOST_AUTO_TEST_SUITE(gaussian_fitter)

BOOST_AUTO_TEST_CASE(fit) {
  for (size_t beam_phase_angle_index = 0; beam_phase_angle_index != 10;
       ++beam_phase_angle_index) {
    const size_t width = 512, height = 512;
    aocommon::Image model(width, height, 0.0);
    aocommon::Image restored(width, height, 0.0);
    model[((height / 2) * width) + (width / 2)] = 1.0;
    const long double kPixelSize = 1.0L /*amin*/ * (M_PI / 180.0 / 60.0);
    const long double beam_major = 20.0L * kPixelSize;
    const long double beam_minor = 5.0L * kPixelSize;
    const long double beam_phase_angle = beam_phase_angle_index * M_PI / 10.0;

    schaapcommon::fft::MakeFftwfPlannerThreadSafe();
    schaapcommon::fft::RestoreImage(
        restored.Data(), model.Data(), width, height, beam_major, beam_minor,
        beam_phase_angle, kPixelSize, kPixelSize, kThreadCount);

    schaapcommon::fitters::GaussianFitter fitter;
    // Check that Fit2DGaussianCentred updates these variables.
    double fit_major = 0.0;  // Should become 20.0.
    double fit_minor = 0.0;  // Should become 5.0.
    double fit_phase_angle =
        -1.0;  // Should become beam_phase_angle, which can be 0.0.
    fitter.Fit2DGaussianCentred(restored.Data(), width, height, 5.0, fit_major,
                                fit_minor, fit_phase_angle, 10.0, false);
    fit_phase_angle = std::fmod((fit_phase_angle + 2.0 * M_PI), M_PI);

    BOOST_CHECK_CLOSE_FRACTION(fit_major, 20.0, 1.0e-3);
    BOOST_CHECK_CLOSE_FRACTION(fit_minor, 5.0, 1.0e-3);
    BOOST_CHECK_CLOSE_FRACTION(fit_phase_angle, beam_phase_angle, 1.0e-3);
  }
}

BOOST_AUTO_TEST_CASE(fit_circular) {
  aocommon::Image model(kWidth, kHeight, 0.0);
  aocommon::Image restored(kWidth, kHeight, 0.0);

  model[((kHeight / 2) * kWidth) + (kWidth / 2)] = 1.0;

  const long double beam_major = 4.0L * kPixelSize;
  const long double beam_minor = 4.0L * kPixelSize;
  const long double beam_phase_angle = 0.0;
  const long double estimated_beam_pixel = 1.0;  // this is on purpose way off
  schaapcommon::fft::MakeFftwfPlannerThreadSafe();
  schaapcommon::fft::RestoreImage(
      restored.Data(), model.Data(), kWidth, kHeight, beam_major, beam_minor,
      beam_phase_angle, kPixelSize, kPixelSize, kThreadCount);

  schaapcommon::fitters::GaussianFitter fitter;
  // Check that Fit2DGaussianCentred updates these variables.
  double fit_major = 0.0;         // Should become 4.0.
  double fit_minor = 0.0;         // Should become 4.0.
  double fit_phase_angle = -1.0;  // Should become beam_phase_angle (0.0).
  fitter.Fit2DGaussianCentred(
      restored.Data(), restored.Width(), restored.Height(),
      estimated_beam_pixel, fit_major, fit_minor, fit_phase_angle, 10.0, false);

  BOOST_CHECK_CLOSE_FRACTION(fit_major, 4.0, 1.0e-4);
  BOOST_CHECK_CLOSE_FRACTION(fit_minor, 4.0, 1.0e-4);
  BOOST_CHECK_SMALL(
      std::abs(fit_phase_angle - static_cast<double>(beam_phase_angle)),
      1.0e-4);
}

BOOST_AUTO_TEST_CASE(fit_small_beam) {
  aocommon::Image model(kWidth, kHeight, 0.0);
  aocommon::Image restored(kWidth, kHeight, 0.0);

  model[((kHeight / 2) * kWidth) + (kWidth / 2)] = 1.0;

  const long double beam_major = 4.0L * kPixelSize;
  const long double beam_minor = 0.5L * kPixelSize;
  const long double beam_phase_angle = 0.0;
  const long double estimated_beam_pixel = 1.0;  // this is on purpose way off

  schaapcommon::fft::MakeFftwfPlannerThreadSafe();
  schaapcommon::fft::RestoreImage(
      restored.Data(), model.Data(), kWidth, kHeight, beam_major, beam_minor,
      beam_phase_angle, kPixelSize, kPixelSize, kThreadCount);

  schaapcommon::fitters::GaussianFitter fitter;
  // Check that Fit2DGaussianCentred updates these variables.
  double fit_major = 0.0;         // Should become 4.0.
  double fit_minor = 0.0;         // Should become 0.5.
  double fit_phase_angle = -1.0;  // Should become beam_phase_angle (0.0).
  fitter.Fit2DGaussianCentred(
      restored.Data(), restored.Width(), restored.Height(),
      estimated_beam_pixel, fit_major, fit_minor, fit_phase_angle, 10.0, false);

  BOOST_CHECK_CLOSE_FRACTION(fit_major, 4.0, 1.0e-4);
  BOOST_CHECK_CLOSE_FRACTION(fit_minor, 0.5, 1.0e-4);
  BOOST_CHECK_SMALL(
      std::abs(fit_phase_angle - static_cast<double>(beam_phase_angle)),
      1.0e-4);
}

BOOST_AUTO_TEST_SUITE_END()
