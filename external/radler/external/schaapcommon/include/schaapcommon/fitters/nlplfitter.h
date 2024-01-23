// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_FITTERS_NLPL_FITTER_H_
#define SCHAAPCOMMON_FITTERS_NLPL_FITTER_H_

#include <cmath>
#include <memory>
#include <vector>

namespace schaapcommon {
namespace fitters {

/**
 * This class fits a power law to a set of points. Note that there is a
 * linear solution for this problem, but the linear solution requires
 * all values to be positive, which is not the case for e.g. spectral
 * energy distributions, because these have noise.
 * This fitter does not have this requirement.
 */
class NonLinearPowerLawFitter {
 public:
  using NumT = float;

  NonLinearPowerLawFitter();

  ~NonLinearPowerLawFitter();

  void AddDataPoint(NumT x, NumT y);

  void Fit(NumT& exponent, NumT& factor);

  void Fit(NumT& a, NumT& b, NumT& c);

  /**
   * @param [out] terms The resulting terms.
   * Using a pre-allocated vector instead of a return value avoids
   * memory allocations in this performance-critical function.
   */
  void Fit(std::vector<NumT>& terms, size_t n_terms);
  void FitStable(std::vector<NumT>& terms, size_t n_terms);

  void FastFit(NumT& exponent, NumT& factor);

  static NumT Evaluate(NumT x, const std::vector<NumT>& terms,
                       NumT reference_frequency_hz = 1.0);

  static long double Evaluate(long double factor, long double exponent,
                              long double frequency_hz) {
    return factor * std::pow(frequency_hz, exponent);
  }

 private:
  void FitImplementation(std::vector<NumT>& terms, size_t n_terms);

  std::unique_ptr<class NLPLFitterData> data_;
};

}  // namespace fitters
}  // namespace schaapcommon

#endif
