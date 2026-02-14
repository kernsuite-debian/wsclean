// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "h5cache.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <hdf5.h>

#include <xtensor/xview.hpp>

#include "gridinterpolate.h"

namespace schaapcommon {
namespace h5parm {

H5Cache::H5Cache(H5::Group group) : SolTab(group) { FillFlaggedValues(); }

std::vector<double> H5Cache::GetValues(const std::string& antenna_name,
                                       const std::vector<double>& times,
                                       const std::vector<double>& frequencies,
                                       size_t polarization, size_t direction,
                                       bool nearest) const {
  TimesAndFrequencies times_and_frequencies = GetTimesAndFrequencies(
      times, frequencies, polarization, direction, nearest);

  const std::vector<double> values =
      GetSubArray(antenna_name, times_and_frequencies.start_time_index,
                  times_and_frequencies.num_times, 1,
                  times_and_frequencies.start_freq_index,
                  times_and_frequencies.num_freqs, 1, polarization, direction);

  MemoryLayout mem_layout = MemoryLayout::kRowMajor;
  // If the frequency index is lower than the time index, time will be the
  // fastest changing index. The ordering needs to be swapped, to ensure that
  // the frequency will be the fastest changing index
  if (HasAxis("freq") && HasAxis("time") &&
      GetAxisIndex("freq") < GetAxisIndex("time")) {
    mem_layout = MemoryLayout::kColumnMajor;
  }

  return GridNearestNeighbor(times_and_frequencies.time_axis,
                             times_and_frequencies.freq_axis, times,
                             frequencies, values, mem_layout, nearest);
}

std::vector<double> H5Cache::GetSubArray(const std::string& antenna_name,
                                         size_t start_time, size_t n_times,
                                         size_t time_step, size_t start_freq,
                                         size_t n_freqs, size_t freq_step,
                                         size_t polarization,
                                         size_t direction) const {
  assert(time_step >= 1);
  assert(freq_step >= 1);

  std::vector<double> sub_array;
  sub_array.resize(n_times * n_freqs);

  size_t sub_index = 0;
  const int antenna = GetAntIndex(antenna_name);
  for (size_t time = start_time; time < (start_time + n_times * time_step);
       time += time_step) {
    for (size_t frequency = start_freq;
         frequency < (start_freq + n_freqs * freq_step);
         frequency += freq_step) {
      sub_array[sub_index] =
          flagged_values_(time, frequency, antenna, direction, polarization);
      sub_index++;
    }
  }

  return sub_array;
}

/**
 * The in-memory array stores the data in the canonical axis order:
 *   [time, freq, ant, dir, pol]
 *
 * If the file data has the canonical order of axes then copying the
 * data from file to the in-memory array would be equivalent to:
 *
 *   iterate over all indexes:
 *   in-memory-array[time, frequency, antenna, direction, polarization] =
 *     file_data[time, frequency, antenna, direction, polarization]
 *
 * The actual code is simpler, since the underlying storage
 * is copied directly from file to memory. That is possible because
 * they have the same order.
 *
 * If file data have a different order we would need to do something
 * equivalent to:
 *
 *   iterate over all indexes:
 *   in-memory-array[time, frequency, antenna, direction, polarization] =
 *     file_data[antenna, time, direction, polarization, frequency]
 *
 * Notice that indexes do not have the same place in each array.
 * In this is a particular case 'ant' is the first axis, 'time' the second,
 * etc. To apply this mechanism for any possible order of the file data we
 * use a mapping of variables that have the name of an axis and it "points"
 * to the actual order of the axis. For instance, in this particular case:
 *
 *   time_axis_         = 1;
 *   frequency_axis_    = 4;
 *   antenna_axis_      = 0;
 *   direction_axis_    = 2;
 *   polarization_axis_ = 3;
 *
 * Applying this idea to copy the data from file to in-memory array would
 * be equivalent to:
 *
 *   size_t& time = axis_index_[time_axis_];
 *   size_t& frequency = axis_index_[frequency_axis_];
 *   size_t& antenna = axis_index_[antenna_axis_];
 *   size_t& direction = axis_index_[direction_axis_];
 *   size_t& polarization = axis_index_[polarization_axis_];
 *
 *   iterate over all indexes:
 *   in-memory-array[time, frequency, antenna, direction, polarization] =
 *     file_data[axis_index_[0], axis_index_[1], axis_index_[2],
 *              axis_index_[3], axis_index_[4]]
 *
 * Assigning the order of axis in file data to *_axis_ variables
 * is done in the method MapAxes().
 *
 * The same idea is applied to keep the size of each axis. The actual value
 * of each dimension is set in the method GetAxisSizes().
 */
void H5Cache::FillFlaggedValues() {
  // Wiring of pointers to file indexes and initialize variable values.
  MapAxes();
  GetAxisSizes();
  const bool canonical_order = HasCanonicalOrder();

  // Read complete arrays and apply flags.
  const size_t n_axes = NumAxes();
  std::vector<double> values =
      SolTab::GetCompleteArray("val", n_axes, axis_size_);
  const std::vector<double> weights =
      SolTab::GetCompleteArray("weight", n_axes, axis_size_);
  ApplyFlags(values, weights);

  // Resize the in-memory array.
  const std::array<size_t, 5> shape = {
      axis_size_[time_axis_], axis_size_[frequency_axis_],
      axis_size_[antenna_axis_], axis_size_[direction_axis_],
      axis_size_[polarization_axis_]};
  flagged_values_.resize(shape);

  // Copy the data from the local container to the in-memory array.
  if (canonical_order) {
    std::copy(values.begin(), values.end(), flagged_values_.begin());
  } else {
    CopyAndReorder(values);
  }
}

void H5Cache::MapAxes() {
  if (NumAxes() < 4) {
    throw std::runtime_error(
        "The solution should have at least 4 axis: time, freq, ant and dir.");
  }

  bool time_found = false;
  bool frequency_found = false;
  bool antenna_found = false;
  bool direction_found = false;
  bool polarization_found = false;

  for (size_t i = 0; i < NumAxes(); i++) {
    const std::string axis_name = GetAxis(i).name;
    if (axis_name == "time") {
      time_axis_ = i;
      time_found = true;
    } else if (axis_name == "freq") {
      frequency_axis_ = i;
      frequency_found = true;
    } else if (axis_name == "ant") {
      antenna_axis_ = i;
      antenna_found = true;
    } else if (axis_name == "dir") {
      direction_axis_ = i;
      direction_found = true;
    } else if (axis_name == "pol") {
      polarization_axis_ = i;
      polarization_found = true;
    } else {
      const std::string message =
          "H5Cache: unrecognized axis name:'" + axis_name + "'";
      throw std::runtime_error(message);
    }
  }

  // If the file data does not have 'pol' axis, then
  // it will be assigned the fifth place.
  if ((NumAxes() == 4) && (!polarization_found)) {
    polarization_axis_ = 4;
    polarization_found = true;
  }

  if (!time_found || !frequency_found || !antenna_found || !direction_found ||
      !polarization_found) {
    throw std::runtime_error("H5Cache: missing one or more axes.");
  }
}

void H5Cache::GetAxisSizes() {
  const std::vector<std::string>& antenna_list = GetStringAxis("ant");

  axis_size_[time_axis_] = GetRealAxis("time").size();
  axis_size_[frequency_axis_] = GetRealAxis("freq").size();
  axis_size_[antenna_axis_] = antenna_list.size();
  axis_size_[direction_axis_] = GetStringAxis("dir").size();

  axis_size_[polarization_axis_] = 2;  // default value
  if (!HasAxis("pol")) {
    axis_size_[polarization_axis_] = 1;
  } else {
    try {
      axis_size_[polarization_axis_] = GetStringAxis("pol").size();
    } catch (const std::exception& e) {
      // There is no dataset corresponding to 'pol' axis, even though
      // the axis exists. Ignore and use the default value.
    }
  }
}

bool H5Cache::HasCanonicalOrder() {
  const std::vector<std::string> canonical_names = {"time", "freq", "ant",
                                                    "dir", "pol"};
  for (size_t i = 0; i < NumAxes(); i++) {
    const std::string axis_name = GetAxis(i).name;
    if (axis_name != canonical_names[i]) {
      return false;
    }
  }

  return true;
}

void H5Cache::CopyAndReorder(const std::vector<double>& values) {
  size_t& time = axis_index_[time_axis_];
  size_t& frequency = axis_index_[frequency_axis_];
  size_t& antenna = axis_index_[antenna_axis_];
  size_t& direction = axis_index_[direction_axis_];
  size_t& polarization = axis_index_[polarization_axis_];

  for (time = 0; time < axis_size_[time_axis_]; ++time) {
    for (frequency = 0; frequency < axis_size_[frequency_axis_]; ++frequency) {
      for (antenna = 0; antenna < axis_size_[antenna_axis_]; ++antenna) {
        for (direction = 0; direction < axis_size_[direction_axis_];
             ++direction) {
          for (polarization = 0; polarization < axis_size_[polarization_axis_];
               ++polarization) {
            const size_t index =
                axis_index_[0] * axis_size_[1] * axis_size_[2] * axis_size_[3] *
                    axis_size_[4] +
                axis_index_[1] * axis_size_[2] * axis_size_[3] * axis_size_[4] +
                axis_index_[2] * axis_size_[3] * axis_size_[4] +
                axis_index_[3] * axis_size_[4] + axis_index_[4];

            flagged_values_(time, frequency, antenna, direction, polarization) =
                values[index];
          }
        }
      }
    }
  }
}

}  // namespace h5parm
}  // namespace schaapcommon
