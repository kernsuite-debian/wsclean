// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "spectralfitter.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

using schaapcommon::fitters::SpectralFitter;
using schaapcommon::fitters::SpectralFittingMode;

namespace {
const double kDefaultReferenceFrequency = 150.0e6;
}

BOOST_AUTO_TEST_SUITE(spectral_fitter)

BOOST_AUTO_TEST_CASE(constructor_default_arguments) {
  const SpectralFittingMode kMode = SpectralFittingMode::kLogPolynomial;
  const size_t kNTerms = 42;

  const SpectralFitter fitter(kMode, kNTerms);
  // BOOST_TEST needs an operator<<(ostream&) for SpectralFittingMode.
  BOOST_TEST((fitter.Mode() == kMode));
  BOOST_TEST(fitter.NTerms() == kNTerms);
  BOOST_TEST(fitter.Frequencies().empty());
  BOOST_TEST(fitter.Weights().empty());
  BOOST_TEST(fitter.ReferenceFrequency() == kDefaultReferenceFrequency);
}

BOOST_AUTO_TEST_CASE(constructor_frequencies_and_weights) {
  const std::vector<double> kFrequencies{10.0, 15.0, 19.0};
  const std::vector<float> kWeights{2.0, 0.0, 1.0};
  const double kReferenceFrequency = 13.0;  // (2*10.0 + 1*19.0) / (2+1)

  const SpectralFitter fitter(SpectralFittingMode::kNoFitting, 0, kFrequencies,
                              kWeights);

  BOOST_CHECK_EQUAL_COLLECTIONS(kFrequencies.begin(), kFrequencies.end(),
                                fitter.Frequencies().begin(),
                                fitter.Frequencies().end());
  BOOST_CHECK_EQUAL_COLLECTIONS(kWeights.begin(), kWeights.end(),
                                fitter.Weights().begin(),
                                fitter.Weights().end());
  BOOST_TEST(fitter.ReferenceFrequency() == kReferenceFrequency);
}

BOOST_AUTO_TEST_CASE(constructor_negative_weight) {
  const std::vector<double> kFrequencies{42.0};
  const std::vector<float> kWeights{-1.5};

  const SpectralFitter fitter(SpectralFittingMode::kNoFitting, 0, kFrequencies,
                              kWeights);

  BOOST_CHECK_EQUAL_COLLECTIONS(kFrequencies.begin(), kFrequencies.end(),
                                fitter.Frequencies().begin(),
                                fitter.Frequencies().end());
  BOOST_CHECK_EQUAL_COLLECTIONS(kWeights.begin(), kWeights.end(),
                                fitter.Weights().begin(),
                                fitter.Weights().end());
  BOOST_TEST(fitter.ReferenceFrequency() == kDefaultReferenceFrequency);
}

BOOST_AUTO_TEST_CASE(constructor_frequencies_weights_mismatch) {
  const std::vector<double> kFrequencies;
  const std::vector<float> kWeights{1.0};

  BOOST_CHECK_THROW(
      const SpectralFitter fitter(SpectralFittingMode::kPolynomial, 0,
                                  kFrequencies, kWeights),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(set_forced_terms) {
  using ImageVector = std::vector<aocommon::Image>;
  const size_t kNTerms = 5;
  const ImageVector kTerms(kNTerms - 1);
  const ImageVector kNotEnoughTerms(kNTerms - 2);

  SpectralFitter fitter(SpectralFittingMode::kForcedTerms, kNTerms);
  BOOST_CHECK_NO_THROW(fitter.SetForcedTerms(ImageVector(kTerms)));
  BOOST_CHECK_THROW(fitter.SetForcedTerms(ImageVector(kNotEnoughTerms)),
                    std::invalid_argument);

  SpectralFitter incorrect_mode(SpectralFittingMode::kLogPolynomial, kNTerms);
  BOOST_CHECK_THROW(incorrect_mode.SetForcedTerms(ImageVector(kTerms)),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(forced_fit) {
  const size_t kNTerms = 2;
  SpectralFitter fitter(SpectralFittingMode::kForcedTerms, kNTerms,
                        {100e6, 110e6}, {1.0, 1.0});
  const size_t kWidth = 10;
  const size_t kHeight = 20;
  std::vector<aocommon::Image> spectral_terms(kNTerms - 1);
  spectral_terms[0] = aocommon::Image(kWidth, kHeight, 1.0);
  fitter.SetForcedTerms(std::move(spectral_terms));
  const std::array<float, 2> values = {7.0, 7.0};
  std::vector<float> fitted_terms(2);
  // Fit the last pixel of the image
  fitter.Fit(fitted_terms, values.data(), kWidth - 1, kHeight - 1);
  BOOST_CHECK_CLOSE_FRACTION(fitted_terms[0], 7.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(fitted_terms[1], 1.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(forced_fit_and_evaluate) {
  const size_t kNTerms = 2;
  SpectralFitter fitter(SpectralFittingMode::kForcedTerms, kNTerms,
                        {100e6, 110e6}, {1.0, 1.0});
  const size_t kWidth = 5;
  const size_t kHeight = 5;
  std::vector<aocommon::Image> spectral_terms(kNTerms - 1);
  spectral_terms[0] = aocommon::Image(kWidth, kHeight, 1.0);
  // Set the term at position (0, 0) to 0
  spectral_terms[0][0] = 0.0;
  fitter.SetForcedTerms(std::move(spectral_terms));
  std::array<float, 2> values = {40.0, 44.0};
  std::vector<float> fitted_terms(2);
  // Fit and evaluate a pixel at position (0, 0)
  fitter.FitAndEvaluate(values.data(), 0, 0, fitted_terms);
  BOOST_CHECK_CLOSE_FRACTION(fitted_terms[0], 42.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(fitted_terms[1], 0.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(values[0], 42.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(values[1], 42.0, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
