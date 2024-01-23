#Allow overriding XTensor versions, e.g., for testing a new version in DP3.
#For avoiding ODR violations, repositories that use aocommon should not override
#these versions in their master branch. That way, the XTensor versions will
#be equal in all repositories.
if (NOT xtl_GIT_TAG)
  set(xtl_GIT_TAG b3d0091a77af52f1b479b5b768260be4873aa8a7)
endif()
if (NOT xsimd_GIT_TAG)
  set(xsimd_GIT_TAG 2f5eddf8912c7e2527f0c50895c7560b964d29af)
endif()
if (NOT xtensor_GIT_TAG)
  set(xtensor_GIT_TAG 0.24.2)
endif()
if (NOT xtensor-blas_GIT_TAG)
  set(xtensor-blas_GIT_TAG 0.20.0)
endif()
if (NOT xtensor-fftw_GIT_TAG)
  set(xtensor-fftw_GIT_TAG e6be85a376624da10629b6525c81759e02020308)
endif()

# By default, only load the basic 'xtensor' and 'xtl' modules.
if (NOT XTENSOR_LIBRARIES)
  set(XTENSOR_LIBRARIES xtl xtensor) # Load xtl first, since xtensor uses it.
endif()

include(FetchContent)

foreach(LIB ${XTENSOR_LIBRARIES})
  set(XT_GIT_TAG "${${LIB}_GIT_TAG}")
  if (NOT XT_GIT_TAG)
    message(FATAL_ERROR "Unknown git tag for XTensor library '${LIB}'")
  endif()

  # Checking out a specific git commit hash does not (always) work when
  # GIT_SHALLOW is TRUE. See the documentation for GIT_TAG in
  # https://cmake.org/cmake/help/latest/module/ExternalProject.html
  # -> If the GIT_TAG is a commit hash, use a non-shallow clone.
  string(LENGTH "${XT_GIT_TAG}" XT_TAG_LENGTH)
  set(XT_SHALLOW TRUE)
  if(XT_TAG_LENGTH EQUAL 40 AND XT_GIT_TAG MATCHES "^[0-9a-f]+$")
    set(XT_SHALLOW FALSE)
  endif()

  FetchContent_Declare(
    ${LIB}
    GIT_REPOSITORY https://github.com/xtensor-stack/${LIB}.git
    GIT_SHALLOW ${XT_SHALLOW}
    GIT_TAG ${XT_GIT_TAG})

  if ("${LIB}" STREQUAL "xtensor-fftw")
    # Unlike the other libraries, xtensor-fftw does not define a CMake target.
    # Its CMakeLists.txt also loads FFTW using custom options.
    # -> Do not build this library, and define an INTERFACE target manually.
    FetchContent_GetProperties(${LIB})
    if(NOT ${${LIB}_POPULATED})
      FetchContent_Populate(${LIB})
    endif()
    add_library(${LIB} INTERFACE)
    target_include_directories(${LIB} SYSTEM INTERFACE "${${LIB}_SOURCE_DIR}/include")
  else()
    FetchContent_MakeAvailable(${LIB})
    # Ensure XTensor headers are included as system headers.
    get_target_property(IID ${LIB} INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(${LIB} PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${IID}")
  endif()

endforeach()
