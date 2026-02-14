// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "reordering.h"

#include <boost/test/unit_test.hpp>

#include <aocommon/polarization.h>

#include <vector>
#include <cstddef>

using aocommon::Polarization;
using schaapcommon::reordering::GetFilenamePrefix;
using schaapcommon::reordering::GetMetaFilename;
using schaapcommon::reordering::GetPartPrefix;

BOOST_AUTO_TEST_SUITE(reordering_utils)

BOOST_AUTO_TEST_CASE(filename_prefix) {
  const std::string filename_prefix = GetFilenamePrefix("test.ms", "");
  BOOST_CHECK_EQUAL(filename_prefix, "test.ms");
}

BOOST_AUTO_TEST_CASE(filename_prefix_tmp_dir) {
  const std::string filename_prefix = GetFilenamePrefix("tmp/test.ms", "tmp");
  BOOST_CHECK_EQUAL(filename_prefix, "tmp/test.ms");
}

BOOST_AUTO_TEST_CASE(filenameprefix_remove_trailing_separator) {
  const std::string filename_prefix = GetFilenamePrefix("test.ms/", "");
  BOOST_CHECK_EQUAL(filename_prefix, "test.ms");
}

BOOST_AUTO_TEST_CASE(metafilename) {
  const std::string meta_filename = GetMetaFilename("test.ms", "", 0);
  BOOST_CHECK_EQUAL(meta_filename, "test.ms-spw0-parted-meta.tmp");
}

BOOST_AUTO_TEST_CASE(metafilename_with_tmp) {
  const std::string meta_filename = GetMetaFilename("test.ms", "tmp", 0);
  BOOST_CHECK_EQUAL(meta_filename, "tmp/test.ms-spw0-parted-meta.tmp");
}

BOOST_AUTO_TEST_CASE(metafilename_ddi) {
  const std::string meta_filename = GetMetaFilename("test.ms", "", 1);
  BOOST_CHECK_EQUAL(meta_filename, "test.ms-spw1-parted-meta.tmp");
}

BOOST_AUTO_TEST_CASE(partprefix) {
  const std::string partprefix =
      GetPartPrefix("test.ms", 0, Polarization::StokesI, 0, "");
  BOOST_CHECK_EQUAL(partprefix, "test.ms-part0000-I-b0");
}

BOOST_AUTO_TEST_SUITE_END()
