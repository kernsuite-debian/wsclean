# Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

find_package(
  Boost
  COMPONENTS filesystem date_time
  REQUIRED)

target_include_directories(${SCHAAPCOMMON_PROJECT_NAME} SYSTEM
                           PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(${SCHAAPCOMMON_PROJECT_NAME} ${Boost_FILESYSTEM_LIBRARY}
                      ${Boost_DATE_TIME_LIBRARY} Boost::boost)

get_filename_component(MODULE ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(PUBLIC_HEADER_DIR ${SCHAAPCOMMON_SOURCE_DIR}/include/schaapcommon/${MODULE})

set(PUBLIC_HEADERS
    ${PUBLIC_HEADER_DIR}/msselection.h
    ${PUBLIC_HEADER_DIR}/reorderedfilewriter.h
    ${PUBLIC_HEADER_DIR}/reorderedhandle.h ${PUBLIC_HEADER_DIR}/reordering.h)

target_sources(
  ${SCHAAPCOMMON_PROJECT_NAME}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/msselection.cc
          ${CMAKE_CURRENT_SOURCE_DIR}/reorderedfilewriter.cc
          ${CMAKE_CURRENT_SOURCE_DIR}/reorderedhandle.cc
          ${CMAKE_CURRENT_SOURCE_DIR}/reordering.cc)

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
