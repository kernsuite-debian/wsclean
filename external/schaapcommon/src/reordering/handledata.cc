// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "handledata.h"
#include "reordering.h"

#include <exception>

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include <aocommon/logger.h>

using aocommon::Logger;

namespace schaapcommon::reordering {

HandleData::~HandleData() {
  // Skip writing back data if we are in the middle of handling an exception
  // (stack unwinding)
  if (std::uncaught_exceptions()) {
    Logger::Info << "An exception occurred, writing back will be skipped.\n";
  } else {
    // We can't throw inside destructor, so catch potential exceptions that
    // occur during writing the measurement sets.
    try {
      if (!is_copy_ && model_update_required_) cleanup_callback_(*this);
    } catch (std::exception& exception) {
      Logger::Error << "Error occurred while finishing IO task: "
                    << exception.what()
                    << "\nMeasurement set might not have been updated.\n";
    }
  }

  if (!is_copy_ && !keep_temporary_files_) {
    Logger::Info << "Cleaning up temporary files...\n";
    std::set<size_t> meta_file_ids;
    for (size_t part = 0; part != channels_.size(); ++part) {
      for (aocommon::PolarizationEnum p : polarizations_) {
        std::string prefix =
            GetPartPrefix(ms_path_, part, p, temporary_directory_);
        std::remove((prefix + ".tmp").c_str());
        std::remove((prefix + "-w.tmp").c_str());
        std::remove((prefix + "-m.tmp").c_str());
      }
      meta_file_ids.emplace(metadata_indices_[part]);
    }
    for (size_t meta_data_id : meta_file_ids) {
      std::string meta_file =
          GetMetaFilename(ms_path_, temporary_directory_, meta_data_id);
      std::remove(meta_file.c_str());
    }
  }
}

void HandleData::Serialize(aocommon::SerialOStream& stream) const {
  stream.String(ms_path_)
      .String(data_column_name_)
      .String(model_column_name_)
      .UInt32(static_cast<uint32_t>(model_storage_manager_))
      .String(temporary_directory_)
      .UInt64(channels_.size());
  for (const aocommon::VectorMap<ChannelRange>& file_ranges : channels_) {
    stream.UInt64(file_ranges.Size());
    for (const ChannelRange& range : file_ranges) {
      stream.UInt64(range.data_desc_id).UInt64(range.start).UInt64(range.end);
    }
  }
  stream.Bool(initial_model_required_)
      .Bool(model_update_required_)
      .UInt64(polarizations_.size());
  for (aocommon::PolarizationEnum p : polarizations_) {
    stream.UInt32(p);
  }
  selection_.Serialize(stream);
  stream.ObjectVector(bands_per_part_);
  stream.UInt64(metadata_indices_.size());
  for (const std::pair<const size_t, size_t>& item : metadata_indices_)
    stream.UInt64(item.first).UInt64(item.second);
  stream.UInt64(n_antennas_);
  // skip is_copy_; gets set by Unserialize()
  stream.Bool(keep_temporary_files_);
}

void HandleData::Unserialize(aocommon::SerialIStream& stream) {
  stream.String(ms_path_)
      .String(data_column_name_)
      .String(model_column_name_)
      .UInt32(model_storage_manager_)
      .String(temporary_directory_);
  channels_.resize(stream.UInt64());
  for (aocommon::VectorMap<ChannelRange>& file_ranges : channels_) {
    file_ranges.Resize(stream.UInt64());
    for (ChannelRange& range : file_ranges) {
      stream.UInt64(range.data_desc_id).UInt64(range.start).UInt64(range.end);
    }
  }
  stream.Bool(initial_model_required_).Bool(model_update_required_);
  size_t n_pol = stream.UInt64();
  polarizations_.clear();
  for (size_t i = 0; i != n_pol; ++i) {
    polarizations_.emplace((aocommon::PolarizationEnum)stream.UInt32());
  }
  selection_.Unserialize(stream);
  stream.ObjectVector(bands_per_part_);
  size_t n_indices = stream.UInt64();
  metadata_indices_.clear();
  for (size_t i = 0; i != n_indices; ++i) {
    const size_t part = stream.UInt64();
    const size_t metafile_index = stream.UInt64();
    metadata_indices_.emplace(part, metafile_index);
  }
  stream.UInt64(n_antennas_);
  is_copy_ = true;
  stream.Bool(keep_temporary_files_);
  // cleanup_callback_ is not serialized as it is a signal only relevant for the
  // original data.
}

}  // namespace schaapcommon::reordering
