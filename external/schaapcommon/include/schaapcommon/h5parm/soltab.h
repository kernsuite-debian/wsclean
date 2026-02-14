// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_H5PARM_SOLTAB_H_
#define SCHAAPCOMMON_H5PARM_SOLTAB_H_

#include <string>
#include <vector>
#include <complex>
#include <map>

#include <H5Cpp.h>

namespace schaapcommon {
namespace h5parm {

/// A name and the length of an exis, e.g. ('freq', 800) for 800 frequencies
struct AxisInfo {
  std::string name;
  size_t size;
};

struct TimesAndFrequencies {
  std::vector<double> time_axis;
  size_t start_time_index;
  size_t num_times;
  std::vector<double> freq_axis;
  size_t start_freq_index;
  size_t num_freqs;
};

/**
 * Determine the axis values (e.g. times) that fall within a requested
 * range. The unit tests give a number of examples which may help
 * understanding.
 * @param axis a sorted non-empty vector with axis values. This vector is
 * reused as return value to avoid an allocation.
 * @param requested_start First axis value that is needed.
 * @param requested_end Last axis value that is needed.
 * @param nearest If @c true, the requested values are rounded to the nearest
 * axis values. If @c false, it is assumed that interpolation is required.
 * This makes sure that there are samples available at the beginning and the
 * end of the result for interpolation. Particularly, the first selected
 * value is equal or lower than @p requested_start, and the last
 * selected value is equal or higher than @p requested_end, if such values
 * exist. Otherwise, the most extreme values are returned.
 * @param [out] start_index index of the first axis value returned.
 * @returns sorted vector with axis values that correspond with the requested
 * selection.
 */
std::vector<double> SelectAxisValues(std::vector<double>&& axis,
                                     double requested_start,
                                     double requested_end, bool nearest,
                                     size_t& start_index);

/// @brief SolTab is a solution table as defined in the H5Parm standard. It
/// contains one solution, e.g. all TEC values, with different axes
/// for that solution (e.g. time, freq, pol).
class SolTab : private H5::Group {
 public:
  SolTab() = default;

  /// Create a new soltab, add it to its parent
  SolTab(H5::Group group, const std::string& type,
         const std::vector<AxisInfo>& axes  /// Axes, fastest varying last
  );

  /// Create a soltab from a H5::Group (for reading existing files)
  explicit SolTab(H5::Group& group);

  /// The destructor could check for valid subtables
  ~SolTab() override;

  /// Add a version stamp in the attributes of the group
  static void AddVersionStamp(H5::Group& node);

  std::vector<AxisInfo>& GetAxes() { return axes_; }

  const std::vector<AxisInfo>& GetAxes() const { return axes_; }

  AxisInfo GetAxis(size_t i) const;

  /// Get an axis, throw an exception if it does not exist
  AxisInfo GetAxis(const std::string& axis_name) const;

  size_t NumAxes() const { return axes_.size(); }

  bool HasAxis(const std::string& axis_name) const;

  /// Get the index of an axis
  size_t GetAxisIndex(const std::string& axis_name) const;

  void SetAntennas(const std::vector<std::string>& sol_antennas);

  void SetSources(const std::vector<std::string>& sol_sources);

  void SetPolarizations(const std::vector<std::string>& polarizations);

  void SetFreqs(const std::vector<double>& freqs);

  /// Get the values of a real-valued axis (e.g. "time" or "freq")
  std::vector<double> GetRealAxis(const std::string& axis_name) const;

  /// Get the values of a string-valued axis (e.g. "dir" or "pol")
  /// @param axis_name Axis name. Only "ant" and "dir" are supported.
  /// @return The requested values, in their original order.
  /// @throw std::runtime_error If the axis name is not supported.
  const std::vector<std::string>& GetStringAxis(
      const std::string& axis_name) const;

  /// Get the index of freq, using nearest neighbor.
  /// This function assumes that the frequencies are in increasing order.
  /// @throw std::runtime_error If the frequency is less than a full frequency
  /// width below the lowest or above the highest frequency
  hsize_t GetFreqIndex(double freq) const;

  /// Get the index of a time. Matches with 0.5*timeInterval
  /// This assumes that all times are regularly spaced
  hsize_t GetTimeIndex(double time) const;

  /// Get the index for an antenna name.
  hsize_t GetAntIndex(const std::string& ant_name) const;

  /// Get the index for a direction name.
  hsize_t GetDirIndex(const std::string& direction_name) const;

  /// Gets the interval (in s.) between a time slot (default first) and
  /// the next. Throws error if there is only one time slot.
  double GetTimeInterval(size_t start = 0) const {
    return GetInterval("time", start);
  }

  /// Gets the interval (in s.) between a channel (default first) and
  /// the next. Throws error if there is only one frequency.
  double GetFreqInterval(size_t start = 0) const {
    return GetInterval("freq", start);
  }

  void SetTimes(const std::vector<double>& times);

  /// Set metadata about an axis (like freq or time))
  void SetAxisMeta(const std::string& meta_name,
                   const std::vector<double>& meta_vals);

  /// Set metadata about an axis (like polarization, direction)
  void SetAxisMeta(const std::string& meta_name, size_t str_len,
                   const std::vector<std::string>& meta_vals);

  /// Adds a real solution.
  /// If weights are emtpy, write ones everywhere
  void SetValues(const std::vector<double>& vals,
                 const std::vector<double>& weights,
                 const std::string& history = "");

  /// Add a complex solution, taking either amplitude or phase
  void SetComplexValues(const std::vector<std::complex<double>>& vals,
                        const std::vector<double>& weights, bool to_amplitudes,
                        const std::string& history = "");

  /// Get the name of this SolTab
  std::string GetName() const;

  std::string GetType() const { return type_; }

  /// Get the values of this SolTab for a given antenna.
  std::vector<double> GetValues(const std::string& ant_name,
                                size_t starttimeslot, size_t ntime,
                                size_t timestep, size_t startfreq, size_t nfreq,
                                size_t freqstep, size_t pol, size_t dir) const {
    return GetSubArray("val", ant_name, starttimeslot, ntime, timestep,
                       startfreq, nfreq, freqstep, pol, dir);
  }

  /// Get the weights of this SolTab for a given antenna.
  std::vector<double> GetWeights(const std::string& ant_name,
                                 size_t starttimeslot, size_t ntime,
                                 size_t timestep, size_t startfreq,
                                 size_t nfreq, size_t freqstep, size_t pol,
                                 size_t dir) const {
    return GetSubArray("weight", ant_name, starttimeslot, ntime, timestep,
                       startfreq, nfreq, freqstep, pol, dir);
  }

  /// Get the values of this SolTab for a given antenna.
  /// Also reads the corresponding weights, which should be one or zero.
  /// If a weight is zero, the corresponding value becomes NaN.
  virtual std::vector<double> GetValues(const std::string& antenna_name,
                                        const std::vector<double>& times,
                                        const std::vector<double>& frequencies,
                                        size_t polarization, size_t direction,
                                        bool nearest) const;

  // Find the nearest time and frequency indices from H5 file for the
  // time and frequencies in the arguments.
  TimesAndFrequencies GetTimesAndFrequencies(const std::vector<double>& times,
                                             const std::vector<double>& freqs,
                                             size_t pol, size_t dir,
                                             bool nearest) const;

  /// Get the values or weights of this SolTab for a given antenna given an
  /// antenna name, a direction index, and a (range of) times and frequencies.
  /// In the returned vector, the freq will be the fastest changing index,
  /// irrespective of the axis ordering in the underlying h5 data structure.
  std::vector<double> GetSubArray(const std::string& val_or_weight,
                                  const std::string& antenna_name,
                                  size_t start_time_index, size_t n_times,
                                  size_t time_step, size_t start_freq_index,
                                  size_t n_freqs, size_t freq_step, size_t pol,
                                  size_t dir) const;

  std::vector<double> GetCompleteArray(
      const std::string& value_or_weight, const size_t n_axes,
      const std::array<size_t, 5>& size_axis) const;

  std::vector<double> GetValuesOrWeights(const std::string& val_or_weight,
                                         const std::string& antenna_name,
                                         const std::vector<double>& times,
                                         const std::vector<double>& frequencies,
                                         size_t polarization, size_t direction,
                                         bool nearest) const;

  /**
   * Depending on 'weights' some 'values' are discarded,
   * i.e. replaced by not-a-number.
   */
  static void ApplyFlags(std::vector<double>& values,
                         const std::vector<double>& weights);

 private:
  static double TakeAbs(std::complex<double> c) { return std::abs(c); }
  static double TakeArg(std::complex<double> c) { return std::arg(c); }

  void ReadAxes();

  void FillCache(std::vector<std::string>& list,
                 std::map<std::string, hsize_t>& map,
                 const std::string& table_name) const;

  /// Get the interval of the axis axis_name
  double GetInterval(const std::string& axis_name, size_t start = 0) const;
  hsize_t GetNamedIndex(std::vector<std::string>& list,
                        std::map<std::string, hsize_t>& map,
                        const std::string& table_name,
                        const std::string& element_name) const;

  std::string type_;
  std::vector<AxisInfo> axes_;

  // The entries below are mutable since they implement caching.
  mutable std::vector<std::string> ant_list_;
  mutable std::vector<std::string> dir_list_;
  mutable std::vector<std::string> pol_list_;
  mutable std::map<std::string, hsize_t> ant_map_;
  mutable std::map<std::string, hsize_t> dir_map_;
  mutable std::map<std::string, hsize_t> pol_map_;
};
}  // namespace h5parm
}  // namespace schaapcommon

#endif
