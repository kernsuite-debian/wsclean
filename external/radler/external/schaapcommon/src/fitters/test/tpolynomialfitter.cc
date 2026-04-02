// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "polynomialfitter.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

using schaapcommon::fitters::PolynomialFitter;

BOOST_AUTO_TEST_SUITE(polynomial_fitter)

BOOST_AUTO_TEST_CASE(fit) {
  PolynomialFitter fitter;
  std::vector<float> terms;
  fitter.AddDataPoint(0.0, 0.0, 1.0);
  fitter.AddDataPoint(1.0, 0.0, 1.0);
  fitter.AddDataPoint(2.0, 1.0, 1.0);
  fitter.AddDataPoint(3.0, 2.0, 1.0);
  fitter.Fit(terms, 2);

  BOOST_CHECK_CLOSE_FRACTION(terms[0], -0.3, 1.0e-3);
  BOOST_CHECK_CLOSE_FRACTION(terms[1], 0.7, 1.0e-3);
}

BOOST_AUTO_TEST_CASE(fit_weighted) {
  PolynomialFitter fitter;
  std::vector<float> terms;
  fitter.AddDataPoint(0.0, 0.0, 1.5);
  fitter.AddDataPoint(1.0, 0.0, 1.5);
  fitter.AddDataPoint(2.0, 1.0, 0.5);
  fitter.AddDataPoint(3.0, 1.0, 0.5);
  fitter.Fit(terms, 1);

  BOOST_CHECK_CLOSE_FRACTION(terms[0], 0.25, 1.0e-3);
}

BOOST_AUTO_TEST_CASE(evaluate) {
  PolynomialFitter fitter;
  std::vector<float> terms;
  terms.push_back(-0.3);
  terms.push_back(0.7);

  BOOST_CHECK_CLOSE_FRACTION(PolynomialFitter::Evaluate(0.0, terms), -0.3,
                             1.0e-3);
  BOOST_CHECK_CLOSE_FRACTION(PolynomialFitter::Evaluate(1.0, terms), 0.4,
                             1.0e-3);
  BOOST_CHECK_CLOSE_FRACTION(PolynomialFitter::Evaluate(2.0, terms), 1.1,
                             1.0e-3);
  BOOST_CHECK_CLOSE_FRACTION(PolynomialFitter::Evaluate(3.0, terms), 1.8,
                             1.0e-3);
}

BOOST_AUTO_TEST_CASE(convert_nlpl_to_polynomial) {
  std::vector<float> pl_terms = {1.0, 2.0, -1.0};
  const float pl_reference_frequency = 15.0;

  const float polynomial_reference_frequency = 30;
  const float minimum_frequency = 25.0;
  const float maximum_frequency = 35.0;

  std::vector<float> terms(2, 0.0);
  schaapcommon::fitters::PowerLawToPolynomialCoefficients(
      terms, pl_terms, pl_reference_frequency, polynomial_reference_frequency,
      minimum_frequency, maximum_frequency);

  const std::vector<float> expected_terms = {3.2, 4.5};
  BOOST_CHECK_CLOSE_FRACTION(terms[0], expected_terms[0], 0.1);
  BOOST_CHECK_CLOSE_FRACTION(terms[1], expected_terms[1], 0.1);

  const float single_term_minimum_frequency = 140.0;
  const float single_term_maximum_frequency = 160.0;

  std::vector<float> single_term(1, 0.0);
  schaapcommon::fitters::PowerLawToPolynomialCoefficients(
      single_term, pl_terms, pl_reference_frequency,
      polynomial_reference_frequency, single_term_minimum_frequency,
      single_term_maximum_frequency);

  const float expected_term = 10;
  BOOST_CHECK_CLOSE_FRACTION(single_term[0], expected_term, 0.1);
}

BOOST_AUTO_TEST_CASE(shift_polynomial_reference_frequency) {
  std::vector<float> terms = {10.0, 0.9, 3.0, -0.1};
  const float input_reference_frequency = 100;
  const float output_reference_frequency = 200.0;

  schaapcommon::fitters::ShiftPolynomialReferenceFrequency(
      terms, input_reference_frequency, output_reference_frequency);

  const std::vector<float> expected_terms = {13.8, 13.2, 10.8, -0.8};
  BOOST_CHECK_CLOSE_FRACTION(terms[0], expected_terms[0], 0.1);
  BOOST_CHECK_CLOSE_FRACTION(terms[1], expected_terms[1], 0.1);
  BOOST_CHECK_CLOSE_FRACTION(terms[2], expected_terms[2], 0.1);
  BOOST_CHECK_CLOSE_FRACTION(terms[3], expected_terms[3], 0.1);
}

BOOST_AUTO_TEST_SUITE_END()
