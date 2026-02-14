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

# By default, only load the basic 'xtensor' library.
if (NOT XTENSOR_LIBRARIES)
  set(XTENSOR_LIBRARIES xtensor)
endif()

# The 'xtensor' library requires the 'xtl' library.
if (NOT xtl IN_LIST XTENSOR_LIBRARIES)
  list(APPEND XTENSOR_LIBRARIES xtl)
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

  # FetchContent_MakeAvailable makes ${LIB} part of the project.
  # Headers from ${LIB} are then installed along with the project.
  # However, most projects only use ${LIB} internally, at compile time,
  # and should not install ${LIB}, including its headers:
  # - For libraries, XTensor shouldn't be part the public API.
  # - For applications, installing headers isn't needed at all.
  #
  # Instead of FetchContent_MakeAvailable, we therefore use
  # FetchContent_Populate and define an INTERFACE target manually.
  # This approach also supports xtensor-fftw, which does not define a CMake
  # target and also loads FFTW using custom options

  FetchContent_GetProperties(${LIB})
  if(NOT ${${LIB}_POPULATED})
    FetchContent_Populate(${LIB})
  endif()
  add_library(${LIB} INTERFACE)
  target_include_directories(${LIB} SYSTEM INTERFACE "${${LIB}_SOURCE_DIR}/include")
endforeach()

# Since xtensor uses xtl and possibly xsimd headers, create dependencies.
target_link_libraries(xtensor INTERFACE xtl)
if (xsimd IN_LIST XTENSOR_LIBRARIES)
  target_link_libraries(xtensor INTERFACE xsimd)
endif()
