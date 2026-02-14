// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "msselection.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include <aocommon/logger.h>

#include <limits>

namespace schaapcommon::reordering {

const size_t MSSelection::kAllFields = std::numeric_limits<size_t>::max();

void MSSelection::Serialize(aocommon::SerialOStream& stream) const {
  stream.VectorUInt64(field_ids_)
      .UInt64(band_id_)
      .UInt64(start_channel_)
      .UInt64(end_channel_)
      .UInt64(start_timestep_)
      .UInt64(end_timestep_)
      .Double(min_uvw_in_m_)
      .Double(max_uvw_in_m_)
      .Bool(auto_correlations_)
      .UInt32(even_odd_selection_);
}

void MSSelection::Unserialize(aocommon::SerialIStream& stream) {
  stream.VectorUInt64(field_ids_)
      .UInt64(band_id_)
      .UInt64(start_channel_)
      .UInt64(end_channel_)
      .UInt64(start_timestep_)
      .UInt64(end_timestep_)
      .Double(min_uvw_in_m_)
      .Double(max_uvw_in_m_)
      .Bool(auto_correlations_)
      .UInt32(even_odd_selection_);
}

}  // namespace schaapcommon::reordering
