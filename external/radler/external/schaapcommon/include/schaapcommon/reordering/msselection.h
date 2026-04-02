// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_MSSELECTION_H_
#define SCHAAPCOMMON_MSSELECTION_H_

#include <aocommon/io/serialstreamfwd.h>
#include <aocommon/multibanddata.h>

#include <cstring>
#include <vector>

namespace schaapcommon::reordering {

class MSSelection {
 public:
  enum EvenOddSelection { kAllTimesteps, kEvenTimesteps, kOddTimesteps };

  const static size_t kAllFields;

  bool HasChannelRange() const { return end_channel_ != 0; }
  bool HasInterval() const { return end_timestep_ != 0; }
  bool HasMinUVWInM() const { return min_uvw_in_m_ != 0.0; }
  bool HasMaxUVWInM() const { return max_uvw_in_m_ != 0.0; }

  size_t BandId() const { return band_id_; }

  size_t ChannelRangeStart() const { return start_channel_; }
  size_t ChannelRangeEnd() const { return end_channel_; }

  size_t IntervalStart() const { return start_timestep_; }
  size_t IntervalEnd() const { return end_timestep_; }

  const std::vector<size_t>& FieldIds() const { return field_ids_; }

  double MinUVWInM() const { return min_uvw_in_m_; }
  double MaxUVWInM() const { return max_uvw_in_m_; }

  bool IsSelected(size_t field_id, size_t timestep, size_t antenna1,
                  size_t antenna2, const double* uvw) const {
    if (HasMinUVWInM() || HasMaxUVWInM()) {
      double u = uvw[0], v = uvw[1], w = uvw[2];
      return IsSelected(field_id, timestep, antenna1, antenna2,
                        std::sqrt(u * u + v * v + w * w));
    } else {
      return IsSelected(field_id, timestep, antenna1, antenna2, 0.0);
    }
  }

  bool IsSelected(size_t field_id, size_t timestep, size_t antenna1,
                  size_t antenna2, double uvw_in_meters) const {
    if (!IsFieldSelected(field_id)) {
      return false;
    } else if (HasInterval() &&
               (timestep < start_timestep_ || timestep >= end_timestep_)) {
      return false;
    } else if (!auto_correlations_ && (antenna1 == antenna2)) {
      return false;
    } else if (HasMinUVWInM() && uvw_in_meters < min_uvw_in_m_) {
      return false;
    } else if (HasMaxUVWInM() && uvw_in_meters > max_uvw_in_m_) {
      return false;
    } else if (even_odd_selection_ != kAllTimesteps) {
      if (even_odd_selection_ == kEvenTimesteps && timestep % 2 != 0) {
        return false;
      } else if (even_odd_selection_ == kOddTimesteps && timestep % 2 != 1) {
        return false;
      }
    }
    return true;
  }

  bool IsFieldSelected(size_t field_id) const {
    return std::find(field_ids_.begin(), field_ids_.end(), field_id) !=
               field_ids_.end() ||
           field_ids_[0] == kAllFields;
  }

  bool IsTimeSelected(size_t timestep) {
    if (HasInterval() &&
        (timestep < start_timestep_ || timestep >= end_timestep_)) {
      return false;
    } else if (even_odd_selection_ != kAllTimesteps) {
      if (even_odd_selection_ == kEvenTimesteps && timestep % 2 != 0) {
        return false;
      } else if (even_odd_selection_ == kOddTimesteps && timestep % 2 != 1) {
        return false;
      }
    }
    return true;
  }

  void SetFieldIds(const std::vector<size_t>& field_ids) {
    field_ids_ = field_ids;
  }
  void SetBandId(size_t band_id) { band_id_ = band_id; }
  void SetChannelRange(size_t start_channel, size_t end_channel) {
    start_channel_ = start_channel;
    end_channel_ = end_channel;
  }
  void SetNoChannelRange() {
    start_channel_ = 0;
    end_channel_ = 0;
  }
  void SetInterval(size_t start_timestep, size_t end_timestep) {
    start_timestep_ = start_timestep;
    end_timestep_ = end_timestep;
  }
  void SetMinUVWInM(double min_uvw) { min_uvw_in_m_ = min_uvw; }
  void SetMaxUVWInM(double max_uvw) { max_uvw_in_m_ = max_uvw; }
  void SetEvenOrOddTimesteps(EvenOddSelection even_or_odd) {
    even_odd_selection_ = even_or_odd;
  }
  EvenOddSelection EvenOrOddTimesteps() const { return even_odd_selection_; }

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);

 private:
  std::vector<size_t> field_ids_ = {0};
  size_t band_id_ = 0;
  size_t start_channel_ = 0;
  size_t end_channel_ = 0;
  size_t start_timestep_ = 0;
  size_t end_timestep_ = 0;
  double min_uvw_in_m_ = 0.0;
  double max_uvw_in_m_ = 0;
  bool auto_correlations_ = false;
  enum EvenOddSelection even_odd_selection_ = kAllTimesteps;
};

}  // namespace schaapcommon::reordering

#endif
