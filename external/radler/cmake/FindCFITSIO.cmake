# * Try to find CFITSIO. Variables used by this module: CFITSIO_ROOT_DIR     -
#   CFITSIO root directory Variables defined by this module: CFITSIO_FOUND -
#   system has CFITSIO CFITSIO_INCLUDE_DIR  - the CFITSIO include directory
#   (cached) CFITSIO_INCLUDE_DIRS - the CFITSIO include directories (identical
#   to CFITSIO_INCLUDE_DIR) CFITSIO_LIBRARY      - the CFITSIO library (cached)
#   CFITSIO_LIBRARIES    - the CFITSIO libraries (identical to CFITSIO_LIBRARY)
#   CFITSIO_VERSION_STRING the found version of CFITSIO, padded to 3 digits

# SPDX-License-Identifier: LGPL-3.0-only

if(NOT CFITSIO_FOUND)

  find_path(
    CFITSIO_INCLUDE_DIR fitsio.h
    HINTS ${CFITSIO_ROOT_DIR}
    PATH_SUFFIXES include include/cfitsio include/libcfitsio0)

  if(CFITSIO_INCLUDE_DIR)
    file(READ "${CFITSIO_INCLUDE_DIR}/fitsio.h" CFITSIO_H)
    set(CFITSIO_VERSION_REGEX
        ".*#define CFITSIO_VERSION[^0-9]*([0-9]+)\\.([0-9]+).*")
    if("${CFITSIO_H}" MATCHES ${CFITSIO_VERSION_REGEX})
      # Pad CFITSIO minor version to three digit because 3.181 is older than
      # 3.35
      string(REGEX REPLACE ${CFITSIO_VERSION_REGEX} "\\1.\\200"
                           CFITSIO_VERSION_STRING "${CFITSIO_H}")
      string(SUBSTRING ${CFITSIO_VERSION_STRING} 0 5 CFITSIO_VERSION_STRING)
      string(REGEX REPLACE "^([0-9]+)[.]([0-9]+)" "\\1" CFITSIO_VERSION_MAJOR
                           ${CFITSIO_VERSION_STRING})
      # CFITSIO_VERSION_MINOR will contain 80 for 3.08, 181 for 3.181 and 200
      # for 3.2
      string(REGEX REPLACE "^([0-9]+)[.]0*([0-9]+)" "\\2" CFITSIO_VERSION_MINOR
                           ${CFITSIO_VERSION_STRING})
    else()
      set(CFITSIO_VERSION_STRING "Unknown")
    endif()
  endif(CFITSIO_INCLUDE_DIR)

  find_library(
    CFITSIO_LIBRARY cfitsio
    HINTS ${CFITSIO_ROOT_DIR}
    PATH_SUFFIXES lib)
  find_library(M_LIBRARY m)
  mark_as_advanced(CFITSIO_INCLUDE_DIR CFITSIO_LIBRARY M_LIBRARY)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(
    CFITSIO
    REQUIRED_VARS CFITSIO_LIBRARY M_LIBRARY CFITSIO_INCLUDE_DIR
    VERSION_VAR CFITSIO_VERSION_STRING)

  set(CFITSIO_INCLUDE_DIRS ${CFITSIO_INCLUDE_DIR})
  set(CFITSIO_LIBRARIES ${CFITSIO_LIBRARY} ${M_LIBRARY})

endif(NOT CFITSIO_FOUND)
