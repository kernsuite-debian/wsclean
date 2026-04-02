// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "channelrange.h"

#include <boost/test/unit_test.hpp>

namespace schaapcommon::reordering {

BOOST_AUTO_TEST_SUITE(reordering_channel_range)

BOOST_AUTO_TEST_CASE(max_channel_range) {
  const std::vector<ChannelRange> kChannelRanges{
      {0, 50, 100},
      {0, 0, 50},
      {0, 100, 200},
      {0, 200, 500},
  };
  const size_t actual = GetMaxChannels(MakeRegularChannelMap(kChannelRanges));
  BOOST_CHECK_EQUAL(actual, 300);
}

BOOST_AUTO_TEST_CASE(max_channel_range_empty_range) {
  const std::vector<ChannelRange> kChannelRanges;
  const size_t actual = GetMaxChannels(MakeRegularChannelMap(kChannelRanges));
  BOOST_CHECK_EQUAL(actual, 0);
}

BOOST_AUTO_TEST_CASE(make_meta_files_map) {
  const std::vector<ChannelRange> kChannelRanges{
      {7, 0, 100}, {8, 0, 50}, {9, 50, 500}};
  std::map<size_t, size_t> part_to_file;
  std::vector<aocommon::OptionalNumber<size_t>> file_to_data_desc_id;
  std::tie(part_to_file, file_to_data_desc_id) =
      MakeMetaFilesMap(MakeRegularChannelMap(kChannelRanges));

  for (size_t i = 0; i != kChannelRanges.size(); ++i) {
    BOOST_CHECK_EQUAL(*file_to_data_desc_id[i], kChannelRanges[i].data_desc_id);
  }

  // links metafile index to file index
  const std::map<size_t, size_t> expected_part_to_file = {
      {0, 0}, {1, 1}, {2, 2}};
  BOOST_CHECK_EQUAL(part_to_file.size(), expected_part_to_file.size());
  for (const std::pair<const size_t, size_t>& entry : expected_part_to_file) {
    const size_t result_file_index = part_to_file.at(entry.first);
    const size_t expected_file_index = entry.second;
    BOOST_CHECK_EQUAL(result_file_index, expected_file_index);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace schaapcommon::reordering
