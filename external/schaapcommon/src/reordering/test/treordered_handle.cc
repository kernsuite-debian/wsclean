// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "reorderedhandle.h"

#include <boost/test/unit_test.hpp>

using schaapcommon::reordering::ChannelRange;
using schaapcommon::reordering::GetDataDescIdMap;
using schaapcommon::reordering::GetMaxChannels;

BOOST_AUTO_TEST_SUITE(reordered_handle)

BOOST_AUTO_TEST_CASE(max_channel_range) {
  const std::vector<ChannelRange> channel_ranges{
      {0, 0, 100},
      {0, 0, 50},
      {0, 100, 200},
      {0, 50, 500},
  };
  const size_t actual = GetMaxChannels(channel_ranges);
  BOOST_CHECK_EQUAL(actual, 450);
}

BOOST_AUTO_TEST_CASE(max_channel_range_empty_range) {
  const std::vector<ChannelRange> channel_ranges;
  const size_t actual = GetMaxChannels(channel_ranges);
  BOOST_CHECK_EQUAL(actual, 0);
}

BOOST_AUTO_TEST_CASE(data_desc_id_map) {
  const std::vector<ChannelRange> channel_ranges{
      {2, 50, 500},
      {0, 0, 100},
      {1, 0, 50},
  };
  const std::map<size_t, size_t> expected = {{0, 1}, {1, 2}, {2, 0}};
  const std::map<size_t, size_t> actual_ddi_map =
      GetDataDescIdMap(channel_ranges);

  for (const auto& ddi_entry : actual_ddi_map) {
    BOOST_CHECK_EQUAL(ddi_entry.second, expected.find(ddi_entry.first)->second);
  }
}

BOOST_AUTO_TEST_SUITE_END()
