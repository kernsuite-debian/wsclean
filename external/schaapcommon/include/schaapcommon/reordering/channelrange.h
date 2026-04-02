#ifndef SCHAAPCOMMON_REORDERED_CHANNEL_RANGE_
#define SCHAAPCOMMON_REORDERED_CHANNEL_RANGE_

#include <cstring>
#include <fstream>
#include <map>

#include <aocommon/optionalnumber.h>
#include <aocommon/uvector.h>
#include <aocommon/vectormap.h>

namespace schaapcommon::reordering {

struct ChannelRange {
  size_t data_desc_id;
  size_t start, end;
  bool operator<(const ChannelRange& rhs) const {
    if (data_desc_id < rhs.data_desc_id) return true;
    if (data_desc_id > rhs.data_desc_id) return false;
    if (start < rhs.start) return true;
    if (start > rhs.start) return false;
    return end < rhs.end;
  }
  bool Empty() const { return start == end; }
};

size_t GetMaxChannels(const aocommon::VectorMap<ChannelRange>& ranges_map);
size_t GetMaxChannels(
    const std::vector<aocommon::VectorMap<ChannelRange>>& channel_ranges);
size_t GetNRanges(const aocommon::VectorMap<ChannelRange>& ranges_map);

inline bool ContainsDataDescId(
    const aocommon::VectorMap<ChannelRange>& ranges_map, size_t data_desc_id) {
  return ranges_map.Size() > data_desc_id && !ranges_map[data_desc_id].Empty();
}

/**
 * Makes two maps: one that links part index to meta file index, and one
 * that links meta file index to data desc id, if the metafile
 * is associated with only one data_desc_id. If a metafile contains multiple
 * data desc ids, the optional value is left unset.
 *
 * While combining these two in one function makes the code a bit dense, one
 * is the by product of the other so it is efficient to make both at once.
 */
std::pair<std::map<size_t, size_t>,
          std::vector<aocommon::OptionalNumber<size_t>>>
MakeMetaFilesMap(
    const std::vector<aocommon::VectorMap<ChannelRange>>& channels);

inline std::vector<aocommon::VectorMap<ChannelRange>> MakeRegularChannelMap(
    const std::vector<ChannelRange>& ranges) {
  std::vector<aocommon::VectorMap<ChannelRange>> result;
  for (const ChannelRange& range : ranges) {
    result.emplace_back().AlwaysEmplace(range.data_desc_id, range);
  }
  return result;
}

}  // namespace schaapcommon::reordering

#endif
