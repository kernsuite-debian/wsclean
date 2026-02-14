# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Top level CMakeLists.txt file for SchaapCommon
cmake_minimum_required(VERSION 3.13) # CMP0079 needs CMake 3.13

set(SCHAAPCOMMON_ALL_MODULES facets fitters h5parm math ducc0 reordering)

# Define a couple of dedicated variables, to prevent that standard CMAKE
# variables get screwed up when schaapcommon included in another package
set(SCHAAPCOMMON_PROJECT_NAME schaapcommon)
set(SCHAAPCOMMON_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

project(${SCHAAPCOMMON_PROJECT_NAME})

option(BUILD_TESTING "Build tests" OFF)

set(CMAKE_MODULE_PATH ${SCHAAPCOMMON_SOURCE_DIR}/cmake
                      ${SCHAAPCOMMON_SOURCE_DIR}/cmake/external)

# CMP0063 is needed to hide symbols in static libraries
cmake_policy(SET CMP0063 NEW)
# Suppress warnings related to ROOT variables
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

# Enable loading dependencies in subdirectories.
cmake_policy(SET CMP0079 NEW)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_POSITION_INDEPENDENT_CODE True)
add_compile_options(
  -Wall
  -Wnon-virtual-dtor
  -Wzero-as-null-pointer-constant
  -Wduplicated-branches
  -Wundef
  -Wvla
  -Wpointer-arith
  -Wextra
  -Wno-unused-parameter)

# Let schaapcommon be a STATIC library, to reduce the risk on linking problems
add_library(schaapcommon STATIC "")

find_package(Threads REQUIRED)
target_link_libraries(${SCHAAPCOMMON_PROJECT_NAME} Threads::Threads)

# Convenience conditional to distinguish between stand-alone vs dependency builds
if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  set(SCHAAPCOMMON_MASTER_PROJECT TRUE)
else()
  set(SCHAAPCOMMON_MASTER_PROJECT FALSE)
endif()

if(SCHAAPCOMMON_MASTER_PROJECT)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  set(CLANG_TIDY_EXE
      "clang-tidy"
      CACHE STRING "clang-tidy executable.")
  set(CLANG_TIDY_FIXES_FILE
      ""
      CACHE STRING "Write clang-tidy fixes to yml file with specified name.")
  option(CLANG_TIDY_WERROR OFF)
  if(CLANG_TIDY_WERROR)
    set(CLANG_TIDY_WARNINGS_AS_ERRORS "*")
  endif()

  file(GLOB_RECURSE ALL_CXX_SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/*.cc)
  add_custom_target(
    check-clang-tidy
    COMMAND
      ${CLANG_TIDY_EXE} --export-fixes=${CLANG_TIDY_FIXES_FILE}
      --warnings-as-errors=${CLANG_TIDY_WARNINGS_AS_ERRORS} -p .
      ${ALL_CXX_SOURCE_FILES})

  # Require AOCOMMON when compiled in standalone.
  set(AOCOMMON_INCLUDE_DIR
      ""
      CACHE FILEPATH "Path to aocommon include")
  if("${AOCOMMON_INCLUDE_DIR}" STREQUAL "")
    message(
      FATAL_ERROR
        "AOCOMMON_INCLUDE_DIR is not defined!\n"
        "Please use 'cmake -DAOCOMMON_INCLUDE_DIR=<aocommon include path>'")
  endif()
  if(NOT EXISTS ${AOCOMMON_INCLUDE_DIR}/aocommon/uvector.h)
    message(
      FATAL_ERROR
        "Can't find uvector.h in AOCOMMON_INCLUDE_DIR (\"${AOCOMMON_INCLUDE_DIR}\")"
    )
  endif()
  include_directories(${AOCOMMON_INCLUDE_DIR})
  include(${AOCOMMON_INCLUDE_DIR}/../CMake/FetchXTensor.cmake)

  # Build everything
  set(SCHAAPCOMMON_MODULES ${SCHAAPCOMMON_ALL_MODULES})

  if(BUILD_TESTING)
    include(CTest)
    find_package(
      Boost
      COMPONENTS unit_test_framework
      REQUIRED)
  endif()
else()
  # Projects using schaapcommon and XTensor should include FetchXTensor.cmake
  # *before* schaapcommon, since:
  # - FetchXTensor.cmake should apply compiler directives to the entire project.
  # - FetchXTensor.cmake does not support duplicate inclusion.

  # Fetch XTensor if it has not been fetched before.
  if(NOT TARGET xtensor)
    include(${CMAKE_SOURCE_DIR}/external/aocommon/CMake/FetchXTensor.cmake)
  endif()

  # Compile in release mode if schaapcommon is included as external dependency.
  set(CMAKE_BUILD_TYPE Release)
endif()

# For backwards compatibility: Build everything if SCHAAPCOMMON_MODULES is unset.
if(NOT SCHAAPCOMMON_MODULES)
  set(SCHAAPCOMMON_MODULES ${SCHAAPCOMMON_ALL_MODULES})
endif()

message(STATUS "Enabled schaapcommon modules: ${SCHAAPCOMMON_MODULES}")

# Source directories
foreach(MODULE ${SCHAAPCOMMON_MODULES})
  if(MODULE STREQUAL "ducc0")
    add_subdirectory(external/${MODULE})
  else()
    add_subdirectory(src/${MODULE})
  endif()
endforeach()
