// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_FITTERS_SPECTRAL_FITTER_H_
#define SCHAAPCOMMON_FITTERS_SPECTRAL_FITTER_H_

#include <optional>
#include <vector>

#include <aocommon/image.h>

#include "rmfitter.h"

namespace schaapcommon {
namespace fitters {

enum class SpectralFittingMode {
  kNoFitting,      /*!< No fitting, each channel gets a separate solution. */
  kPolynomial,     /*!< Use polynomial for spectral fitting. */
  kLogPolynomial,  /*!< Use double log polynomial for spectral fitting. */
  kForcedTerms,    /*!< Use forced terms for spectral fitting. */
  kRotationMeasure /*!< Fit a rotation measure for RM synthesis. */
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
   * fitted. values[i] should correspond to Frequency(i) and Weight(i).
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

  /**
   * This function behaves like @ref FitAndEvaluate(), but expects an array
   * where the values for different polarizations are interlaced. The first
   * n_pol values are for the first channel, the next n_pol for the second
   * channel, etc. For fitting methods in which the polarizations are
   * independently fit (e.g. polynomial fit), the relevant values are selected
   * and passed to FitAndEvaluate(). For RM fitting, n_polarizations is expected
   * to be 2 (for Stokes Q and U), and the RM fitter is called.
   */
  void MultiPolarizationFitAndEvaluate(NumT* values, size_t x, size_t y,
                                       std::vector<NumT>& fitting_scratch,
                                       size_t n_polarizations) {
    const size_t n = frequencies_.size();
    if (mode_ == schaapcommon::fitters::SpectralFittingMode::kRotationMeasure) {
      assert(n_polarizations == 2);
      rm_fitter_->Fit(std::span(values, n * 2));
    } else {
      for (size_t p = 0; p != n_polarizations; ++p) {
        // Values are ordered by pol, so reshuffle so all frequencies are
        // together. It's somewhat like an (inplace) transpose, but then for
        // only one column.
        for (size_t ch = 0; ch != n; ++ch) {
          std::swap(values[ch * n_polarizations + p], values[ch]);
        }
        FitAndEvaluate(values, x, y, fitting_scratch);
        // placing channel values back should be in reversed order to
        // undo multiple moves of a single value that might have happened
        for (size_t i = 0; i != n; ++i) {
          const size_t ch = n - i - 1;
          std::swap(values[ch * n_polarizations + p], values[ch]);
        }
      }
    }
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
  std::optional<RmFitter> rm_fitter_;
};

}  // namespace fitters
}  // namespace schaapcommon

#endif
