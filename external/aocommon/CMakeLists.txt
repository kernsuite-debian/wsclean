# 'FetchContent' in FetchXTensor.cmake requires CMake 3.11.
# 'add_compile_definitions' in FetchXTensor.cmake requires CMake 3.12.
cmake_minimum_required(VERSION 3.12)

project(aocommon LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(
  -O3
  -Wall
  -Wnon-virtual-dtor
  -Wzero-as-null-pointer-constant
  -Wduplicated-branches
  -Wundef
  -Wvla
  -Wpointer-arith
  -Wextra
  -Wno-unused-parameter
  -ggdb)

# Use may optionally set `TARGET_CPU` if `PORTABLE=OFF`
option(PORTABLE "Build portable binaries (with slightly decreased performance)"
       OFF)
include(CMake/SetTargetCPU.cmake)

# add target to generate API documentation with Doxygen
find_package(Threads REQUIRED)
find_package(Doxygen)

if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                 ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(
    doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
endif(DOXYGEN_FOUND)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Casacore has a separate CMake file in this directory
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits)
include_directories(SYSTEM ${CASACORE_INCLUDE_DIRS})

find_package(CFITSIO REQUIRED)
include_directories(SYSTEM ${CFITSIO_INCLUDE_DIR})

include(CMake/FetchXTensor.cmake)

add_executable(
  runtests
  tests/runtests.cpp
  tests/tangle.cpp
  tests/tbanddata.cpp
  tests/tbarrier.cpp
  tests/tvectormap.cpp
  tests/tcounting_semaphore.cpp
  tests/tdynamicfor.cpp
  tests/tfits.cpp
  tests/tfluxdensity.cpp
  tests/thmatrix4x4.cpp
  tests/timage.cpp
  tests/timagecoordinates.cpp
  tests/tlane.cpp
  tests/tlogger.cpp
  tests/tmatrix2x2.cpp
  tests/tmatrix2x2diag.cpp
  tests/tmatrix4x4.cpp
  tests/tmultibanddata.cpp
  tests/toptionalnumber.cpp
  tests/toverlappingtaskprocessor.cpp
  tests/tpolarization.cpp
  tests/tpyuniqueptr.cpp
  tests/tqueue.cpp
  tests/tradeccoord.cpp
  tests/trecursivefor.cpp
  tests/tstaticfor.cpp
  tests/tserialstream.cpp
  tests/tspan.cpp
  tests/ttaskqueue.cpp
  tests/tthreadpool.cpp
  tests/tthrowruntimeerror.cpp
  tests/ttransform_if.cpp
  tests/tutensor.cpp
  tests/tuvector.cpp
  tests/avx/tDiagonalMatrixComplexDouble2x2.cpp
  tests/avx/tDiagonalMatrixComplexFloat2x2.cpp
  tests/avx/tMatrixComplexDouble2x2.cpp
  tests/avx/tMatrixComplexFloat2x2.cpp)
target_link_libraries(runtests ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}
                      ${CASACORE_LIBRARIES} ${CFITSIO_LIBRARY} xtensor)

add_custom_target(execute_runtests_target ALL DEPENDS execute_runtests)
add_custom_command(
  OUTPUT execute_runtests
  COMMAND runtests
  DEPENDS runtests)

add_custom_target(
  coverage
  COMMAND gcovr -r ${CMAKE_SOURCE_DIR} -e '.*/tests/.*' -e '.*/CompilerIdCXX/.*'
          -e '_deps/.*' --txt --xml coverage.xml ${CMAKE_BINARY_DIR}
  DEPENDS execute_runtests)

# CMAKE_CXX_FLAGS does not contain all flags.
# For example, CMake derives C++ version flags from CMAKE_CXX_STANDARD etc.
message(STATUS "Flags passed to C++ compiler (incomplete): " ${CMAKE_CXX_FLAGS})
