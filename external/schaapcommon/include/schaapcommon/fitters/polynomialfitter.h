// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_FITTERS_POLYNOMIAL_FITTER_H_
#define SCHAAPCOMMON_FITTERS_POLYNOMIAL_FITTER_H_

#include <array>
#include <vector>

namespace schaapcommon {
namespace fitters {

class PolynomialFitter {
 public:
  using NumT = float;

  void Clear() { data_points_.clear(); }

  void AddDataPoint(NumT x, NumT y, NumT w) {
    data_points_.emplace_back(std::array<NumT, 3>{{x, y, w}});
  }

  /**
   * @param [out] terms The resulting terms.
   * Using a pre-allocated vector instead of a return value avoids
   * memory allocations in this performance-critical function.
   */
  void Fit(std::vector<NumT>& terms, std::size_t nTerms);

  static NumT Evaluate(NumT x, const std::vector<NumT>& terms) {
    NumT val = terms[0];
    NumT f = 1.0;
    for (std::size_t i = 1; i != terms.size(); ++i) {
      f *= x;
      val += f * terms[i];
    }
    return val;
  }

  std::size_t Size() const { return data_points_.size(); }

 private:
  std::vector<std::array<NumT, 3>> data_points_;
};

/**
 * Convert logarithmic polynomial coefficients to those of a (normal) polynomial
 * of a given order.
 * @param [in,out] terms The size of this vector determines on input the number
 * of coefficients fitted. On output, it contains the calculated fit polynomial
 * terms.
 */
void PowerLawToPolynomialCoefficients(std::vector<float>& terms,
                                      const std::vector<float>& pl_terms,
                                      float pl_reference_frequency_hz,
                                      float polynomial_reference_frequency_hz,
                                      float min_frequency_hz,
                                      float max_frequency_hz);

/**
 * Shift (normal) polynomial to a different reference frequency.
 * @param [in,out] terms The coefficients of the input polynomial. On output, it
 * contains the calculated fit polynomial terms.
 */
void ShiftPolynomialReferenceFrequency(std::vector<float>& terms,
                                       float input_reference_frequency_hz,
                                       float output_reference_frequency_hz);

}  // namespace fitters
}  // namespace schaapcommon

#endif
