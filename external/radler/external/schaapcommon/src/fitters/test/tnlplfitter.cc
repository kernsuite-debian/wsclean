// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "nlplfitter.h"

#include <aocommon/uvector.h>

#include <boost/test/unit_test.hpp>

using schaapcommon::fitters::NonLinearPowerLawFitter;

BOOST_AUTO_TEST_SUITE(nlplfitter)

BOOST_AUTO_TEST_CASE(first_order) {
  const double x_factor = 1.0e1;

  NonLinearPowerLawFitter fitter;
  std::vector<float> terms{1.0, -0.7};

  for (size_t x = 1; x != 10; ++x) {
    float y = NonLinearPowerLawFitter::Evaluate(x * x_factor, terms);
    fitter.AddDataPoint(x * x_factor, y);
  }

  float e = 0.0;
  float factor = 0.0;
  fitter.Fit(e, factor);
  BOOST_CHECK_SMALL(std::fabs(1.0 - factor), 1.0e-6);
  BOOST_CHECK_SMALL(std::fabs(-0.7 - e), 1.0e-6);
}

BOOST_AUTO_TEST_CASE(second_order_zero) {
  const float x_factor = 1.0e1;

  NonLinearPowerLawFitter fitter;
  std::vector<float> terms{1.0, -0.7};

  for (size_t x = 1; x != 10; ++x) {
    float y = NonLinearPowerLawFitter::Evaluate(x * x_factor, terms);
    fitter.AddDataPoint(x * x_factor, y);
  }
  std::vector<float> fitted;
  fitter.Fit(fitted, 3);
  BOOST_CHECK_SMALL(std::fabs(1.0f - fitted[0]), 1.0e-6f);
  BOOST_CHECK_SMALL(std::fabs(-0.7f - fitted[1]), 1.0e-6f);
  BOOST_CHECK_SMALL(std::fabs(0.0f - fitted[2]), 1.0e-6f);
}

BOOST_AUTO_TEST_CASE(first_order_stability) {
  const float x_factor = 1.0e1;

  NonLinearPowerLawFitter fitter;
  std::vector<float> terms{1.0, -0.7, -0.01};
  for (size_t x = 1; x != 10; ++x) {
    float y = NonLinearPowerLawFitter::Evaluate(x * x_factor, terms);
    fitter.AddDataPoint(x * x_factor, y);
  }
  float e = 0.0;
  float factor = 0.0;
  fitter.Fit(e, factor);
  BOOST_CHECK_SMALL(std::fabs(terms[0] - factor), 0.1f);
  BOOST_CHECK_SMALL(std::fabs(terms[1] - e), 0.1f);
}

BOOST_AUTO_TEST_CASE(second_order_nonzero) {
  const float x_factor = 1.0e1;

  NonLinearPowerLawFitter fitter;
  std::vector<float> terms{1.0, -0.7, -0.01};
  for (size_t x = 1; x != 10; ++x) {
    float y = NonLinearPowerLawFitter::Evaluate(x * x_factor, terms);
    fitter.AddDataPoint(x * x_factor, y);
  }
  std::vector<float> fitted;
  fitter.Fit(fitted, 3);
  BOOST_CHECK_SMALL(std::fabs(terms[0] - fitted[0]), 1.0e-3f);
  BOOST_CHECK_SMALL(std::fabs(terms[1] - fitted[1]), 1.0e-3f);
  BOOST_CHECK_SMALL(std::fabs(terms[2] - fitted[2]), 1.0e-3f);
}

BOOST_AUTO_TEST_CASE(third_order) {
  const float x_factor = 1.0e1;

  NonLinearPowerLawFitter fitter;
  std::vector<float> terms{1.0, -0.7, -0.01, 0.05};
  for (size_t x = 1; x != 10; ++x) {
    float y = NonLinearPowerLawFitter::Evaluate(x * x_factor, terms);
    fitter.AddDataPoint(x * x_factor, y);
  }
  std::vector<float> fitted;
  fitter.Fit(fitted, 4);
  BOOST_CHECK_SMALL(std::fabs(terms[0] - fitted[0]), 1.0e-2f);
  BOOST_CHECK_SMALL(std::fabs(terms[1] - fitted[1]), 1.0e-2f);
  BOOST_CHECK_SMALL(std::fabs(terms[2] - fitted[2]), 1.0e-2f);
  BOOST_CHECK_SMALL(std::fabs(terms[3] - fitted[3]), 1.0e-2f);
}

BOOST_AUTO_TEST_SUITE_END()
