# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

include(${SCHAAPCOMMON_SOURCE_DIR}/cmake/unittest.cmake)

add_unittest(
  fitters
  runtests.cc
  tgaussianfitter.cc
  tnlplfitter.cc
  tpolynomialchannelfitter.cc
  tpolynomialfitter.cc
  tspectralfitter.cc)
