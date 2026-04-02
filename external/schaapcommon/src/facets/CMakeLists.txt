# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Find and include Boost (geometry) headers
find_package(Boost REQUIRED)
target_link_libraries(${SCHAAPCOMMON_PROJECT_NAME} Boost::boost)

# Boost_VERSION has 0's instead of periods for older versions of boost
string(REPLACE 0 . Boost_VERSION_TMP ${Boost_VERSION})
if(Boost_VERSION_TMP VERSION_LESS "1.66")
  target_compile_definitions(${SCHAAPCOMMON_PROJECT_NAME}
                             PRIVATE -DHAVE_BOOST_LT_166)
endif()

get_filename_component(MODULE ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(PUBLIC_HEADER_DIR ${SCHAAPCOMMON_SOURCE_DIR}/include/schaapcommon/${MODULE})

set(PUBLIC_HEADERS
    "${PUBLIC_HEADER_DIR}/facet.h" "${PUBLIC_HEADER_DIR}/facetimage.h"
    "${PUBLIC_HEADER_DIR}/boundingbox.h" "${PUBLIC_HEADER_DIR}/ds9facetfile.h")

target_sources(
  ${SCHAAPCOMMON_PROJECT_NAME}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/boundingbox.cc
          ${CMAKE_CURRENT_SOURCE_DIR}/facet.cc
          ${CMAKE_CURRENT_SOURCE_DIR}/facetimage.cc)

# Simplifies including the public headers.
target_include_directories(${SCHAAPCOMMON_PROJECT_NAME}
                           PRIVATE "$<BUILD_INTERFACE:${PUBLIC_HEADER_DIR}>")

# Install headers and add test directory when built as stand-alone.
if(SCHAAPCOMMON_MASTER_PROJECT)
  install(FILES ${PUBLIC_HEADERS}
          DESTINATION "include/${CMAKE_PROJECT_NAME}/${MODULE}")

  if(BUILD_TESTING)
    add_subdirectory(test)
  endif()
endif()
