// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "polynomialfitter.h"
#include "nlplfitter.h"

#include <stdexcept>
#include <cassert>

#include <gsl/gsl_multifit.h>

namespace schaapcommon {
namespace fitters {

void PolynomialFitter::Fit(std::vector<NumT>& terms, size_t nTerms) {
  size_t n = data_points_.size();
  terms.assign(nTerms, 0.0);

  if (nTerms > n) {
    nTerms = n;
  }

  gsl_multifit_linear_workspace* workspace =
      gsl_multifit_linear_alloc(n, nTerms);

  gsl_matrix* xData = gsl_matrix_alloc(n, nTerms);
  gsl_matrix* cov = gsl_matrix_alloc(nTerms, nTerms);
  gsl_vector* wData = gsl_vector_alloc(n);
  gsl_vector* yData = gsl_vector_alloc(n);
  gsl_vector* resultTerms = gsl_vector_alloc(nTerms);
  double chisq;

  for (size_t i = 0; i != n; ++i) {
    NumT x = data_points_[i][0];
    NumT y = data_points_[i][1];
    NumT w = data_points_[i][2];

    NumT f = 1.0;
    gsl_matrix_set(xData, i, 0, f);
    for (size_t j = 1; j != nTerms; ++j) {
      f *= x;  // f = x^j
      gsl_matrix_set(xData, i, j, f);
    }
    gsl_vector_set(yData, i, y);
    gsl_vector_set(wData, i, w);
  }

  int result = gsl_multifit_wlinear(xData, wData, yData, resultTerms, cov,
                                    &chisq, workspace);

  for (size_t j = 0; j != nTerms; ++j) {
    terms[j] = gsl_vector_get(resultTerms, j);
  }

  gsl_vector_free(resultTerms);
  gsl_vector_free(yData);
  gsl_vector_free(wData);

  gsl_matrix_free(cov);
  gsl_matrix_free(xData);

  gsl_multifit_linear_free(workspace);

  if (result) {
    throw std::runtime_error("Polynomial fit failed");
  }
}

void PowerLawToPolynomialCoefficients(std::vector<float>& terms,
                                      const std::vector<float>& pl_terms,
                                      float pl_reference_frequency_hz,
                                      float polynomial_reference_frequency_hz,
                                      float min_frequency_hz,
                                      float max_frequency_hz) {
  assert(max_frequency_hz >= min_frequency_hz);
  assert(!terms.empty());

  const std::size_t n_terms = terms.size();

  std::vector<float> frequencies;
  const float frequency_step =
      terms.size() == 1
          ? 0.5 * (min_frequency_hz + max_frequency_hz)
          : (max_frequency_hz - min_frequency_hz) / (terms.size() - 1);
  for (std::size_t term_index = 0; term_index != n_terms; ++term_index) {
    frequencies.emplace_back(min_frequency_hz + frequency_step * term_index);
  }

  std::vector<float> values;
  for (float frequency : frequencies) {
    const float value = NonLinearPowerLawFitter::Evaluate(
        frequency, pl_terms, pl_reference_frequency_hz);
    values.push_back(value);
  }

  PolynomialFitter fitter;
  for (std::size_t i = 0; i != frequencies.size(); ++i) {
    fitter.AddDataPoint(
        frequencies[i] / polynomial_reference_frequency_hz - 1.0, values[i], 1);
  }

  fitter.Fit(terms, n_terms);
}

void ShiftPolynomialReferenceFrequency(std::vector<float>& terms,
                                       float input_reference_frequency_hz,
                                       float output_reference_frequency_hz) {
  assert(!terms.empty());

  const float min_frequency = 0.5 * std::min(input_reference_frequency_hz,
                                             output_reference_frequency_hz);
  const float max_frequency = 2.0 * std::max(input_reference_frequency_hz,
                                             output_reference_frequency_hz);

  const std::size_t n_terms = terms.size();

  std::vector<float> frequencies;
  const float frequency_step =
      n_terms == 1 ? 0.5 * (min_frequency + max_frequency)
                   : (max_frequency - min_frequency) / (n_terms - 1);
  for (std::size_t term_index = 0; term_index != n_terms; ++term_index) {
    frequencies.emplace_back(min_frequency + frequency_step * term_index);
  }

  std::vector<float> values;
  for (float frequency : frequencies) {
    const float value = PolynomialFitter::Evaluate(
        frequency / input_reference_frequency_hz - 1.0f, terms);
    values.push_back(value);
  }

  PolynomialFitter fitter;
  for (std::size_t i = 0; i != frequencies.size(); ++i) {
    fitter.AddDataPoint(frequencies[i] / output_reference_frequency_hz - 1.0f,
                        values[i], 1);
  }

  fitter.Fit(terms, n_terms);
}

}  // namespace fitters
}  // namespace schaapcommon
