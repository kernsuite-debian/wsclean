// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_H5PARM_H5CACHE_H_
#define SCHAAPCOMMON_H5PARM_H5CACHE_H_

#include <limits>

#include <xtensor/xtensor.hpp>

#include "soltab.h"

namespace schaapcommon {
namespace h5parm {

/**
 * This class is meant to be used as a drop-in replacement
 * of the parent class SolTab for quick access to H5 parameters.
 *
 * However, it is not possible to replace any arbitrary SolTab
 * by H5Cache. A soltab to be cacheable must have at least these
 * axes: ant, time, freq, dir, and pol, and their sizes multiplied
 * must be the size of the val and weight datasets.
 *
 * In the constructor, this class reads all the H5 parameters
 * from the file using parent class methods, flags them and
 * puts them all in the in-memory array 'flagged_values_'.
 * This makes the access to the H5 parameters much faster, at the
 * expense of more memory usage.
 *
 * In a mutithreading system, the constructor must be invoked
 * only once. 'GetValues(..)' could be called concurrently
 * since it only reads data.
 */
class H5Cache : public SolTab {
 public:
  explicit H5Cache(H5::Group group);

  /**
   * This class is meant to store data in an in-memory
   * structure for quick access — potentially large data.
   * Hence, copying this object is disabled.
   */
  H5Cache(const H5Cache&) = delete;             // non construction-copyable
  H5Cache& operator=(const H5Cache&) = delete;  // non copyable

  /**
   * Returns a vector of values already flagged by weights,
   * selected from all the H5 parameters according to direction,
   * times and frequencies. These values are read from the in-memory
   * array flagged_values_, and not from the H5 file, which makes
   * a difference with respect to the parent class.
   */
  std::vector<double> GetValues(const std::string& antenna_name,
                                const std::vector<double>& times,
                                const std::vector<double>& frequencies,
                                size_t polarization, size_t direction,
                                bool nearest) const override;

 private:
  /**
   * Reads all the H5 parameters from the file using parent class
   * methods and put them all in the in-memory array.
   */
  void FillFlaggedValues();

  /**
   * This method does not have the argument 'val_or_weight',
   * as the method with the same name in the parent class,
   * because this class stores the values already flagged, i.e.
   * the weights were read and applied.
   */
  std::vector<double> GetSubArray(const std::string& antenna_name,
                                  size_t start_time, size_t n_times,
                                  size_t time_step, size_t start_freq,
                                  size_t n_freqs, size_t freq_step,
                                  size_t polarization, size_t direction) const;

  void MapAxes();
  void GetAxisSizes();
  bool HasCanonicalOrder();
  void CopyAndReorder(const std::vector<double>& values);
  /**
   * In-memory array to keep all flagged values of the H5
   * solution file. The dimensions are:
   * [time, frequency, antenna, direction, polarization].
   */
  xt::xtensor<double, 5> flagged_values_;

  /**
   * Constant to initialize variables.
   */
  const size_t kInvalidAxis = std::numeric_limits<size_t>::max();

  /**
   * These variables keep the order of the axes
   * of the file data. For instance, if file data
   * is ordered this way: [ant, dir, time, pol, freq]
   * then
   *   time_axis_         = 2;
   *   frequency_axis_    = 4;
   *   antenna_axis_      = 0;
   *   direction_axis_    = 1;
   *   polarization_axis_ = 3;
   * @{
   */
  size_t time_axis_ = kInvalidAxis;
  size_t frequency_axis_ = kInvalidAxis;
  size_t antenna_axis_ = kInvalidAxis;
  size_t direction_axis_ = kInvalidAxis;
  size_t polarization_axis_ = kInvalidAxis;
  /** @} */

  /**
   * This variable is used with the *_axis_ variables
   * to "point" to an specific axis index of the file data.
   * For instance: axis_index_[direction_axis_]
   */
  std::array<size_t, 5> axis_index_;

  /**
   * This variable is used with the *_axis_ variables
   * to "point" to the the size of an specific axis.
   * For instante: axis_size_[antenna_axis_].
   */
  std::array<size_t, 5> axis_size_;
};

}  // namespace h5parm
}  // namespace schaapcommon

#endif
