# Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

include(${SCHAAPCOMMON_SOURCE_DIR}/cmake/unittest.cmake)

add_unittest(
  reordering
  runtests.cc
  tchannelrange.cc
  tfilewriter.cc
  treordering.cc
  treordering_extraction.cc
  treordering_store.cc)
