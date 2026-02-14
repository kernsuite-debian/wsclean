// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_FITTERS_SPECTRAL_FITTER_H_
#define SCHAAPCOMMON_FITTERS_SPECTRAL_FITTER_H_

#include <vector>

#include <aocommon/image.h>

namespace schaapcommon {
namespace fitters {

enum class SpectralFittingMode {
  kNoFitting,     /*!< No fitting, each channel gets a separate solution. */
  kPolynomial,    /*!< Use polynomial for spectral fitting. */
  kLogPolynomial, /*!< Use double log polynomial for spectral fitting. */
  kForcedTerms    /*!< Use forced terms for spectral fitting. */
};

class SpectralFitter {
 public:
  using NumT = float;

  SpectralFitter(SpectralFittingMode mode, size_t n_terms,
                 std::vector<double> frequencies = {},
                 std::vector<NumT> weights = {});

  /**
   * Fit an array of values to a curve.
   *
   * The type of the curve is set in the constructor or with @ref SetMode().
   * The coordinates are used in case the forced term fitting mode is used, in
   * which case it is used to look up the spectral index (or other terms) from
   * a specified image.
   *
   * @param [out] terms will hold the fitted terms. The meaning of these terms
   * depends on the fitted curve type, and are relative to the reference
   * frequency. Using a pre-allocated vector instead of a return value avoids
   * memory allocations in this performance-critical function.
   * @param values array of size @ref NFrequencies() with the values to be
   * fitted. values[i] should correspond Frequency(i) and Weight(i).
   * @param x a pixel index giving the horizontal position
   * @param y a pixel index giving the vertical position
   */
  void Fit(std::vector<NumT>& terms, const NumT* values, size_t x,
           size_t y) const;

  /**
   * Evaluate the curve at the initialized frequencies.
   *
   * @param values array of size @ref NFrequencies() that will be filled with
   * curve values.
   * @param terms array of size @ref NTerms() with previously fitted terms.
   */
  void Evaluate(NumT* values, const std::vector<NumT>& terms) const;

  /**
   * Evaluate the curve at a specified frequency.
   *
   * @param terms array of size @ref NTerms() with previously fitted terms.
   * @param frequency Frequency in Hz.
   */
  NumT Evaluate(const std::vector<NumT>& terms, double frequency) const;

  /**
   * Fit an array of values to a curve, and replace those values
   * with the curve values. This function combines @ref Fit()
   * and @ref Evaluate().
   *
   * @param terms is a vector of any size, that is used to store the terms.
   * Having this parameter explicitly is useful to avoid repeated allocation,
   * to temporarily store the terms: This function is used in reasonably
   * critical loops inside deconvolution. It will be resized to @ref NTerms().
   */
  void FitAndEvaluate(NumT* values, size_t x, size_t y,
                      std::vector<NumT>& terms) const {
    Fit(terms, values, x, y);
    Evaluate(values, terms);
  }

  SpectralFittingMode Mode() const { return mode_; }

  size_t NTerms() const { return n_terms_; }

  const std::vector<double>& Frequencies() const { return frequencies_; };

  const std::vector<NumT>& Weights() const { return weights_; }

  double ReferenceFrequency() const { return reference_frequency_; }

  /**
   * Update the forced terms, when using forced terms for spectral fitting.
   *
   * @param terms New terms. Force using move semantics for this argument, since
   *        images are typically large.
   * @throw std::runtime_error If the mode is not kForcedTerms.
   * @throw std::invalid_argument If terms does not have enough elements.
   */
  void SetForcedTerms(std::vector<aocommon::Image>&& terms);

 private:
  void ForcedFit(std::vector<NumT>& terms, const NumT* values, size_t x,
                 size_t y) const;

  SpectralFittingMode mode_;
  size_t n_terms_;
  std::vector<double> frequencies_;
  std::vector<NumT> weights_;
  double reference_frequency_;
  std::vector<aocommon::Image> forced_terms_;
};

}  // namespace fitters
}  // namespace schaapcommon

#endif
