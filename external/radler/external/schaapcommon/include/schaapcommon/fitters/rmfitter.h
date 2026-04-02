#ifndef RADLER_MATH_RM_FITTER_H_
#define RADLER_MATH_RM_FITTER_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <numeric>
#include <span>
#include <vector>

namespace schaapcommon::fitters {
namespace rm {

/**
 * Convert a frequency in Hz into a wavelength squared number
 * in meter^2.
 */
inline double FrequencyToWavelengthSquared(double frequency) {
  const double wavelength = 299792458.0 / frequency;
  return wavelength * wavelength;
}

/**
 * Get the phase factor corresponding to a given Faraday depth (phi) and
 * wavelength squared value (in m^2).
 */
inline std::complex<float> GetPhaseFactor(double phi,
                                          double wavelength_squared) {
  double phase = -2.0 * phi * wavelength_squared;
  return {static_cast<float>(std::cos(phase)),
          static_cast<float>(std::sin(phase))};
}

/**
 * Convert an array of frequencies in Hz to a vector of wavelength-squared
 * values in meter^2.
 */
inline std::vector<double> MakeWavelengthSquaredVector(
    std::span<const double> frequencies) {
  std::vector<double> wavelength_squared_;
  wavelength_squared_.reserve(frequencies.size());
  for (double f : frequencies) {
    wavelength_squared_.emplace_back(FrequencyToWavelengthSquared(f));
  }
  return wavelength_squared_;
}

/**
 * Determine appropriate values for the Faraday depth grid. The returned value
 * consists of the number of grid points n and the maximum Faraday depth p_max.
 * The Faraday depth of a grid point at index i can then be calculated using:
 *   p_max * (i / n) - p_max / 2.
 * @returns a pair where first is the number of grid points and second is the
 * maximum Faraday depth.
 */
inline std::pair<size_t, double> GetPhiGridParameters(
    std::span<const double> wavelength_squared) {
  double max_phi;
  size_t n_phi;
  if (wavelength_squared.size() > 1) {
    double max_distance =
        wavelength_squared.front() - wavelength_squared.back();
    assert(max_distance > 0);
    double min_distance = std::numeric_limits<double>::max();
    for (size_t i = 0; i + 1 != wavelength_squared.size(); ++i) {
      const double delta = wavelength_squared[i] - wavelength_squared[i + 1];
      assert(delta > 0);
      min_distance = std::min(delta, min_distance);
    }
    // The max frequency is 1 / min_distance. A frequency 'turn' is 2 pi, but
    // the phi-to-phase factor is 2 (without pi). Hence the factor of pi.
    max_phi = M_PI / min_distance;
    // (1/min_distance) / (1/max_distance)
    // = max_phi_ * max_distance
    // The factor 4 is to oversample the phi grid, so that the maximum phi
    // value can be more accurately measured.
    n_phi = std::ceil(4.0 * max_distance * max_phi);

    // Prefer the grid to have an even size, so that rm 0 is
    // exactly on a grid point.
    if (n_phi % 2 != 0) ++n_phi;
  } else {
    max_phi = 0.0;
    n_phi = 1;
  }
  return {n_phi, max_phi};
}

}  // namespace rm

/**
 * Class that can perform rotation measure fitting. It's interface is tailored
 * for usage inside the spectral fitting of clean.
 */
class RmFitter {
 public:
  /**
   * Construct a fitter from an array of frequency values in Hz. The frequencies
   * are used to determine an appropriate Faraday depth grid.
   */
  RmFitter(std::span<const double> frequencies, std::span<const float> weights)
      : weights_(weights.begin(), weights.end()) {
    sum_of_weights_ = std::accumulate(weights.begin(), weights.end(), 0.0f);
    assert(weights_.size() == frequencies.size());
    wavelength_squared_ = rm::MakeWavelengthSquaredVector(frequencies);
    std::tie(n_phi_, max_phi_) = rm::GetPhiGridParameters(wavelength_squared_);
    phasors_ = MakePhasorsMatrix();
  }

  /**
   * Finds the strongest signal in the data, and replaces the data with only
   * that signal. The data array is an array of real and imaginary values.
   * It should have a total of '2 x n_frequencies' values. The data has type
   * float instead of complex to match with how the data is stored in the
   * deconvolution procedures.
   */
  void Fit(std::span<float> data) const {
    const std::pair<size_t, std::complex<float>> value = FindPeakRmIndex(data);
    EvaluateRm(data, value.first, value.second);
  }

  /**
   * Function that returns the Faraday depth and magnitude of the strongest
   * signal. See @ref Fit() for the structure of @p data.
   */
  std::pair<float, float> FindPeakRm(std::span<const float> data) const {
    const std::pair<size_t, std::complex<float>> value = FindPeakRmIndex(data);
    return {GetPhi(value.first), std::abs(value.second)};
  }

 private:
  /**
   * Returns the grid index of the RM value with the largest norm, together with
   * the value of that RM value.
   */
  std::pair<size_t, std::complex<float>> FindPeakRmIndex(
      std::span<const float> data) const {
    assert(data.size() == wavelength_squared_.size() * 2);
    std::vector<std::complex<float>> rm_values(n_phi_, 0.0);
    std::vector<std::complex<float>>::const_iterator phasor = phasors_.begin();
    for (size_t i = 0; i != n_phi_; ++i) {
      for (size_t j = 0; j != wavelength_squared_.size(); ++j) {
        const std::complex<float> value(data[j * 2], data[j * 2 + 1]);
        if (weights_[j] != 0.0) rm_values[i] += value * *phasor * weights_[j];
        ++phasor;
      }
    }
    const std::vector<std::complex<float>>::const_iterator phi_iterator =
        std::max_element(rm_values.begin(), rm_values.end(),
                         [](std::complex<float> lhs, std::complex<float> rhs) {
                           return std::norm(lhs) < std::norm(rhs);
                         });
    return {phi_iterator - rm_values.begin(), *phi_iterator / sum_of_weights_};
  }

  void EvaluateRm(std::span<float> data, size_t phi_index,
                  std::complex<float> rm_value) const {
    std::vector<std::complex<float>>::const_iterator phasor =
        phasors_.begin() + phi_index * wavelength_squared_.size();
    for (size_t j = 0; j != wavelength_squared_.size(); ++j) {
      const std::complex<float> value = std::conj(*phasor) * rm_value;
      data[j * 2] = value.real();
      data[j * 2 + 1] = value.imag();
      ++phasor;
    }
  }

  float GetPhi(size_t grid_index) const {
    return max_phi_ * grid_index / n_phi_ - max_phi_ * 0.5;
  }

  std::vector<std::complex<float>> MakePhasorsMatrix() const {
    std::vector<std::complex<float>> phasors;
    phasors.reserve(n_phi_ * wavelength_squared_.size());
    for (size_t phi_index = 0; phi_index != n_phi_; ++phi_index) {
      for (double w : wavelength_squared_) {
        phasors.emplace_back(rm::GetPhaseFactor(GetPhi(phi_index), w));
      }
    }
    return phasors;
  }

  size_t n_phi_;
  float max_phi_;
  std::vector<double> wavelength_squared_;
  std::vector<float> weights_;
  float sum_of_weights_;
  // matrix of n_channels x n_phi
  std::vector<std::complex<float>> phasors_;
};

}  // namespace schaapcommon::fitters

#endif
