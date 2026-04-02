#include "channelrange.h"

#include <algorithm>
#include <set>

namespace schaapcommon::reordering {
namespace {
bool IsNotEmpty(const ChannelRange& range) { return !range.Empty(); };

std::set<size_t> GetDataDescIds(
    const std::vector<aocommon::VectorMap<ChannelRange>>& channels) {
  std::set<size_t> data_desc_ids;
  for (const aocommon::VectorMap<ChannelRange>& ranges : channels) {
    for (size_t data_desc_id = 0; data_desc_id != ranges.Size();
         ++data_desc_id) {
      if (!ranges[data_desc_id].Empty()) data_desc_ids.emplace(data_desc_id);
    }
  }
  return data_desc_ids;
}

bool IsFullyRegular(
    const std::vector<aocommon::VectorMap<ChannelRange>>& channels) {
  for (const aocommon::VectorMap<ChannelRange>& ranges : channels) {
    const size_t n_data_desc_ids = ranges.CountKeysIf(IsNotEmpty);
    if (n_data_desc_ids > 1) {
      return false;
    }
  }
  return true;
}

}  // namespace

size_t GetMaxChannels(const aocommon::VectorMap<ChannelRange>& ranges) {
  size_t max_channels = 0;
  for (const ChannelRange& range : ranges) {
    max_channels = std::max(max_channels, range.end - range.start);
  }
  return max_channels;
}

size_t GetMaxChannels(
    const std::vector<aocommon::VectorMap<ChannelRange>>& channel_ranges) {
  size_t max_channels = 0;
  for (const aocommon::VectorMap<ChannelRange>& file_ranges : channel_ranges) {
    max_channels = std::max(max_channels, GetMaxChannels(file_ranges));
  }
  return max_channels;
}

bool GetNRanges(const aocommon::VectorMap<ChannelRange>& ranges_map,
                size_t data_desc_id) {
  return std::count_if(
      ranges_map.begin(), ranges_map.end(),
      [](const ChannelRange& range) { return range.start != range.end; });
}

std::pair<std::map<size_t, size_t>,
          std::vector<aocommon::OptionalNumber<size_t>>>
MakeMetaFilesMap(
    const std::vector<aocommon::VectorMap<ChannelRange>>& channels) {
  const std::set<size_t> data_desc_ids = GetDataDescIds(channels);
  const bool is_fully_regular = IsFullyRegular(channels);

  // The map links the part index to meta_file index.
  // The vector links file index to data_desc_id (if the metafile
  // is associated with only one data_desc_id).
  std::pair<std::map<size_t, size_t>,
            std::vector<aocommon::OptionalNumber<size_t>>>
      result;
  std::map<size_t, size_t>& meta_file_indices = result.first;
  std::vector<aocommon::OptionalNumber<size_t>>& data_desc_id_per_file =
      result.second;
  if (is_fully_regular) {
    // Link all parts with the same data_desc_id to the same meta file,
    // i.e. they will share the same meta data. This is possible because
    // data with the same data desc id comes from the same rows.
    for (const size_t data_desc_id : data_desc_ids) {
      bool has_data = false;
      for (size_t part = 0; part != channels.size(); ++part) {
        const aocommon::VectorMap<ChannelRange>& ranges = channels[part];
        if (ContainsDataDescId(ranges, data_desc_id)) {
          has_data = true;
          meta_file_indices[part] = data_desc_id_per_file.size();
        }
      }
      if (has_data) {
        data_desc_id_per_file.emplace_back(data_desc_id);
      }
    }
  } else {
    // Parts with the same data_desc_id, where all use only that data_desc_id
    // can share their meta file.
    std::vector<size_t> parts;
    for (const size_t data_desc_id : data_desc_ids) {
      parts.clear();
      bool can_share = true;
      for (size_t part = 0; part != channels.size(); ++part) {
        const aocommon::VectorMap<ChannelRange>& ranges = channels[part];
        if (ContainsDataDescId(ranges, data_desc_id)) {
          if (ranges.CountKeysIf(IsNotEmpty) > 1) {
            can_share = false;
            break;
          } else {
            parts.emplace_back(part);
          }
        }
      }
      if (can_share) {
        for (size_t part : parts)
          meta_file_indices[part] = data_desc_id_per_file.size();
        data_desc_id_per_file.emplace_back(data_desc_id);
      }
    }
    for (size_t part = 0; part != channels.size(); ++part) {
      if (!meta_file_indices.contains(part)) {
        meta_file_indices[part] = data_desc_id_per_file.size();
        data_desc_id_per_file.emplace_back();
      }
    }
  }
  return result;
}

}  // namespace schaapcommon::reordering
