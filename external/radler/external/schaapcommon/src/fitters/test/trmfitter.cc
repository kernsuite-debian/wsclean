
// SPDX-License-Identifier: LGPL-3.0-only

#include <boost/test/unit_test.hpp>

#include "rmfitter.h"

#include <aocommon/uvector.h>

#include <random>

namespace schaapcommon::fitters {
namespace {
inline std::vector<double> LofarFrequencies() {
  std::vector<double> frequencies;
  for (size_t i = 0; i != 200; ++i) {
    // 120 - 160 MHz
    frequencies.emplace_back(double(i) * 0.2e6 + 120e6);
  }
  return frequencies;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(rm_fitter)

BOOST_AUTO_TEST_CASE(conversion) {
  BOOST_CHECK_CLOSE_FRACTION(1.0, rm::FrequencyToWavelengthSquared(299.8e6),
                             1e-3);
  BOOST_CHECK_CLOSE_FRACTION(16.4126, rm::FrequencyToWavelengthSquared(74e6),
                             1e-3);
}

BOOST_AUTO_TEST_CASE(grid_parameters) {
  const std::vector<double> frequencies = LofarFrequencies();
  const std::vector<double> wavelengths =
      rm::MakeWavelengthSquaredVector(frequencies);
  const std::pair<size_t, float> result = rm::GetPhiGridParameters(wavelengths);
  // Current implementation return 3876, but any reasonable value is okay:
  BOOST_CHECK_GT(result.first, 2000);
  BOOST_CHECK_LT(result.first, 5000);
  // Similarly, implementation return 355.928589, but any reasonable value is
  // okay:
  BOOST_CHECK_GT(result.second, 250);
  BOOST_CHECK_LT(result.second, 600);
}

BOOST_AUTO_TEST_CASE(single_channel) {
  const std::vector<double> frequencies = LofarFrequencies();
  const double frequency = 120e6;
  const float weight = 1.0;
  const RmFitter fitter(std::span(&frequency, 1), std::span(&weight, 1));
}

BOOST_AUTO_TEST_CASE(find_peak_rm_and_fit) {
  const std::vector<double> frequencies = LofarFrequencies();
  const std::vector<double> wavelengths =
      rm::MakeWavelengthSquaredVector(frequencies);
  const RmFitter fitter(frequencies,
                        std::vector<float>(frequencies.size(), 7.0f));
  // Construct an array with two signals with different Faraday depths and
  // brightnesses.
  std::vector<float> data;
  const float magnitude_a = 25.0;
  const float rm_a = 8.0;
  const float magnitude_b = 10.0;
  const float rm_b = 3.0;
  for (float w : wavelengths) {
    const float theta_a = 2.0 * w * rm_a;
    const float theta_b = 2.0 * w * rm_b;
    std::complex<float> value =
        std::polar(magnitude_a, theta_a) + std::polar(magnitude_b, theta_b);
    data.emplace_back(value.real());
    data.emplace_back(value.imag());
  }
  float result_rm;
  float result_magnitude;
  // FindPeakRm should find the strongest signal
  std::tie(result_rm, result_magnitude) = fitter.FindPeakRm(data);
  BOOST_CHECK_GT(result_rm, rm_a - 0.1);
  BOOST_CHECK_LT(result_rm, rm_a + 0.1);
  BOOST_CHECK_GT(result_magnitude, magnitude_a - 0.5);
  BOOST_CHECK_LT(result_magnitude, magnitude_a + 0.5);

  // Fitting the data should replace the data with data from the
  // strongest signal.
  std::vector<float> fitted_data = data;
  fitter.Fit(fitted_data);

  // A first simple test if fitting worked: the new data should still
  // have the strongest signal
  std::tie(result_rm, result_magnitude) = fitter.FindPeakRm(fitted_data);
  BOOST_CHECK_GT(result_rm, rm_a - 0.1);
  BOOST_CHECK_LT(result_rm, rm_a + 0.1);
  BOOST_CHECK_GT(result_magnitude, magnitude_a - 0.5);
  BOOST_CHECK_LT(result_magnitude, magnitude_a + 0.5);

  // Now subtract the fitted data with the strongest signal from the
  // original data with two signals. The result should be that the
  // weaker signal is left.
  for (size_t i = 0; i != data.size(); ++i) {
    data[i] -= fitted_data[i];
  }
  std::tie(result_rm, result_magnitude) = fitter.FindPeakRm(data);
  BOOST_CHECK_GT(result_rm, rm_b - 0.1);
  BOOST_CHECK_LT(result_rm, rm_b + 0.1);
  BOOST_CHECK_GT(result_magnitude, magnitude_b - 0.5);
  BOOST_CHECK_LT(result_magnitude, magnitude_b + 0.5);

  // Test an RM 0 signal
  const float rm_0_magnitude = 19.82;
  data.clear();
  for (size_t i = 0; i != wavelengths.size(); ++i) {
    data.emplace_back(rm_0_magnitude);
    data.emplace_back(0.0);
  }
  std::tie(result_rm, result_magnitude) = fitter.FindPeakRm(data);
  BOOST_CHECK_LT(std::abs(result_rm), 0.1);
  BOOST_CHECK_GT(result_magnitude, rm_0_magnitude - 0.2);
  BOOST_CHECK_LT(result_magnitude, rm_0_magnitude + 0.2);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace schaapcommon::fitters
