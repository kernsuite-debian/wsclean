# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Find HDF5 before CasaCore. Since CasaCore also finds HDF5, HDF5 may not have
# CXX support if CasaCore finds it first.
find_package(
  HDF5
  COMPONENTS C CXX
  REQUIRED)

# In CMake 3.20, the hdf5::hdf5_cpp target can replace the lines below.
target_compile_definitions(${SCHAAPCOMMON_PROJECT_NAME}
                           PUBLIC ${HDF5_DEFINITIONS})
target_include_directories(${SCHAAPCOMMON_PROJECT_NAME} SYSTEM
                           PUBLIC ${HDF5_INCLUDE_DIR})
target_link_libraries(${SCHAAPCOMMON_PROJECT_NAME} ${HDF5_LIBRARIES}
                      ${HDF5_CXX_LIBRARIES})

# Temporarily include casacore until this is replaced with aocommon::MC2x2.
find_package(
  Casacore
  COMPONENTS casa
  REQUIRED)

target_include_directories(${SCHAAPCOMMON_PROJECT_NAME} SYSTEM
                           PUBLIC ${CASACORE_INCLUDE_DIR})
target_link_libraries(${SCHAAPCOMMON_PROJECT_NAME} ${CASACORE_LIBRARIES}
                      xtensor)
get_directory_property(H5PARM_COMPILE_DEFINITIONS COMPILE_DEFINITIONS)
if("XTENSOR_USE_XSIMD" IN_LIST H5PARM_COMPILE_DEFINITIONS)
  target_link_libraries(${SCHAAPCOMMON_PROJECT_NAME} xsimd)
endif()

get_filename_component(MODULE ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(PUBLIC_HEADER_DIR ${SCHAAPCOMMON_SOURCE_DIR}/include/schaapcommon/${MODULE})

set(PUBLIC_HEADERS
    ${PUBLIC_HEADER_DIR}/h5cache.h ${PUBLIC_HEADER_DIR}/h5parm.h
    ${PUBLIC_HEADER_DIR}/jonesparameters.h ${PUBLIC_HEADER_DIR}/soltab.h)

target_sources(
  ${SCHAAPCOMMON_PROJECT_NAME}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/h5cache.cc
          ${CMAKE_CURRENT_SOURCE_DIR}/h5parm.cc
          ${CMAKE_CURRENT_SOURCE_DIR}/jonesparameters.cc
          ${CMAKE_CURRENT_SOURCE_DIR}/soltab.cc
          ${CMAKE_CURRENT_SOURCE_DIR}/gridinterpolate.cc)

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
