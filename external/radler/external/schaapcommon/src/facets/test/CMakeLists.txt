# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

include(${SCHAAPCOMMON_SOURCE_DIR}/cmake/unittest.cmake)

add_unittest(facets runtests.cc tfacet.cc tfacetimage.cc tds9facetfile.cc
             tboundingbox.cc)

# Copy test resources.
foreach(file empty.reg foursources.reg)
  configure_file(resources/${file} ${CMAKE_BINARY_DIR}/src/facets/test/${file}
                 COPYONLY)
endforeach()
