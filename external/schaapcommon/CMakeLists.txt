# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Top level CMakeLists.txt file for SchaapCommon
cmake_minimum_required(VERSION 3.8)

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

# Convenience conditional to distinguish between stand-alone vs dependency builds
set(SCHAAPCOMMON_MASTER_PROJECT OFF)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(SCHAAPCOMMON_MASTER_PROJECT ON)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  set(CLANG_TIDY_EXE
      "clang-tidy"
      CACHE STRING "clang-tidy executable.")
  set(CLANG_TIDY_FIXES_FILE
      ""
      CACHE STRING "Write clang-tidy fixes to yml file with specified name.")
  option(CLANG_TIDY_WERROR OFF)
endif()

if(SCHAAPCOMMON_MASTER_PROJECT)
  # Require the AOCOMMON headers when compiled in standalone.
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

  if(BUILD_TESTING)
    include(CTest)
    find_package(
      Boost
      COMPONENTS unit_test_framework
      REQUIRED)
  endif()
endif()

# Compile in release mode if schaapcommon is included from external
# dependency
if(NOT SCHAAPCOMMON_MASTER_PROJECT)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Spawn all compiler warnings
add_compile_options(-Wall -Wzero-as-null-pointer-constant)

# Load dependencies, could be done a bit neater with CMP0079, i.e.
# cmake_policy(SET CMP0079 NEW)
# but requires CMake >= 3.13

# Find and include HDF5
find_package(
  HDF5
  COMPONENTS C CXX
  REQUIRED)
add_definitions(${HDF5_DEFINITIONS})
include_directories(SYSTEM ${HDF5_INCLUDE_DIR})

# Find and include FFTW
find_library(FFTW3F_LIB fftw3f REQUIRED HINTS ENV FFTW3_LIB)
find_library(FFTW3F_THREADS_LIB fftw3f_threads REQUIRED HINTS ENV FFTW3_LIB)
find_path(
  FFTW3_INCLUDE_DIR
  NAMES fftw3.h
  HINTS ENV FFTW3_INCLUDE)

include_directories(SYSTEM ${FFTW3_INCLUDE_DIR})

# Temporary include casacore until this is replaced with aocommon::MC2x2
find_package(
  Casacore
  COMPONENTS casa
  REQUIRED)
include_directories(SYSTEM ${CASACORE_INCLUDE_DIR})

# Find and include Boost (geometry) headers
find_package(Boost REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# Boost_VERSION has 0's instead of periods for older versions of boost
string(REPLACE 0 . Boost_VERSION_TMP ${Boost_VERSION})
if(Boost_VERSION_TMP VERSION_LESS "1.66")
  add_definitions(-DHAVE_BOOST_LT_166)
endif()

find_package(Threads REQUIRED)

find_package(GSL REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})

# Let schaapcommon be a STATIC library, to reduce the risk on linking problems
add_library(schaapcommon STATIC "")
target_link_libraries(
  ${SCHAAPCOMMON_PROJECT_NAME}
  ${HDF5_LIBRARIES}
  ${HDF5_CXX_LIBRARIES}
  ${CASACORE_LIBRARIES}
  ${FFTW3F_LIB}
  ${FFTW3F_THREADS_LIB}
  ${CMAKE_THREAD_LIBS_INIT}
  ${GSL_LIBRARIES})

if(SCHAAPCOMMON_MASTER_PROJECT)
  if(CMAKE_EXPORT_COMPILE_COMMANDS)
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
  endif()

  # Install targets
  install(
    TARGETS schaapcommon
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)
endif()

# Source directories
add_subdirectory(src/facets)
add_subdirectory(src/fft)
add_subdirectory(src/fitters)
add_subdirectory(src/h5parm)
