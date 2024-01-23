// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_FITTERS_POLYNOMIAL_CHANNEL_FITTER_H_
#define SCHAAPCOMMON_FITTERS_POLYNOMIAL_CHANNEL_FITTER_H_

#include <vector>

namespace schaapcommon {
namespace fitters {

class PolynomialChannelFitter {
 public:
  void Clear() {
    channels_.clear();
    data_points_.clear();
  }

  void AddChannel(double start_frequency, double end_frequency) {
    channels_.emplace_back(start_frequency, end_frequency);
  }

  void AddDataPoint(std::size_t channel, double y) {
    data_points_.emplace_back(channel, y);
  }

  void Fit(std::vector<double>& terms, std::size_t n_terms);

  static double Evaluate(double x, const std::vector<double>& terms);

 private:
  /**
   * Start and end frequencies of the channels
   */
  std::vector<std::pair<double, double>> channels_;
  std::vector<std::pair<std::size_t, double>> data_points_;
};

}  // namespace fitters
}  // namespace schaapcommon

#endif
