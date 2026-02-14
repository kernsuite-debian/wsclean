// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "reorderedhandle.h"
#include "reordering.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include <aocommon/logger.h>

using aocommon::Logger;

namespace schaapcommon::reordering {

ReorderedHandleData::~ReorderedHandleData() {
  if (!(is_copy_ || keep_temporary_files_)) {
    // We can't throw inside destructor, so catch potential exceptions that
    // occur during writing the measurement sets.
    try {
      // Skip writing back data if we are in the middle of handling an exception
      // (stack unwinding)
      if (std::uncaught_exceptions()) {
        Logger::Info << "An exception occurred, writing back will be skipped. "
                        "Cleaning up...\n";
      } else {
        if (model_update_required_) cleanup_callback_(*this);
        Logger::Info << "Cleaning up temporary files...\n";
      }

      std::set<size_t> removed_meta_files;
      for (size_t part = 0; part != channels_.size(); ++part) {
        for (aocommon::PolarizationEnum p : polarizations_) {
          std::string prefix =
              GetPartPrefix(ms_path_, part, p, channels_[part].data_desc_id,
                            temporary_directory_);
          std::remove((prefix + ".tmp").c_str());
          std::remove((prefix + "-w.tmp").c_str());
          std::remove((prefix + "-m.tmp").c_str());
        }
        const size_t data_desc_id = channels_[part].data_desc_id;
        if (removed_meta_files.count(data_desc_id) == 0) {
          removed_meta_files.insert(data_desc_id);
          std::string meta_file =
              GetMetaFilename(ms_path_, temporary_directory_, data_desc_id);
          std::remove(meta_file.c_str());
        }
      }
    } catch (std::exception& exception) {
      Logger::Error << "Error occurred while finishing IO task: "
                    << exception.what()
                    << "\nMeasurement set might not have been updated.\n";
    }
  }
}

void ReorderedHandleData::Serialize(aocommon::SerialOStream& stream) const {
  stream.String(ms_path_)
      .String(data_column_name_)
      .String(model_column_name_)
      .String(temporary_directory_)
      .UInt64(channels_.size());
  for (const ChannelRange& range : channels_) {
    stream.UInt64(range.data_desc_id).UInt64(range.start).UInt64(range.end);
  }
  stream.Bool(initial_model_required_)
      .Bool(model_update_required_)
      .UInt64(polarizations_.size());
  for (aocommon::PolarizationEnum p : polarizations_) {
    stream.UInt32(p);
  }
  selection_.Serialize(stream);
  stream.UInt64(n_antennas_);
}

void ReorderedHandleData::Unserialize(aocommon::SerialIStream& stream) {
  is_copy_ = true;
  stream.String(ms_path_)
      .String(data_column_name_)
      .String(model_column_name_)
      .String(temporary_directory_);
  channels_.resize(stream.UInt64());
  for (ChannelRange& range : channels_) {
    stream.UInt64(range.data_desc_id).UInt64(range.start).UInt64(range.end);
  }
  stream.Bool(initial_model_required_).Bool(model_update_required_);
  size_t n_pol = stream.UInt64();
  polarizations_.clear();
  for (size_t i = 0; i != n_pol; ++i) {
    polarizations_.emplace((aocommon::PolarizationEnum)stream.UInt32());
  }
  selection_.Unserialize(stream);
  stream.UInt64(n_antennas_);
}

size_t GetMaxChannels(const std::vector<ChannelRange>& channel_ranges) {
  size_t max_channels = 0;
  for (const ChannelRange& range : channel_ranges) {
    max_channels = std::max(max_channels, range.end - range.start);
  }
  return max_channels;
}

std::map<size_t, size_t> GetDataDescIdMap(
    const std::vector<ChannelRange>& channels) {
  std::map<size_t, size_t> data_desc_ids;
  size_t spw_index = 0;
  for (const ChannelRange& range : channels) {
    if (data_desc_ids.count(range.data_desc_id) == 0) {
      data_desc_ids.emplace(range.data_desc_id, spw_index);
      ++spw_index;
    }
  }
  return data_desc_ids;
}

}  // namespace schaapcommon::reordering
