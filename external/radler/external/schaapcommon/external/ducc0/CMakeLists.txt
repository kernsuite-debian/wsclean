set(DUCC0_PATH ${CMAKE_CURRENT_SOURCE_DIR})

set(DUCC0_HEADERS ${DUCC0_PATH}/fft/*.h ${DUCC0_PATH}/infra/*.h
                  ${DUCC0_PATH}/math/*.h ${DUCC0_PATH}/wgridder/*.h)
target_sources(
  ${SCHAAPCOMMON_PROJECT_NAME}
  PRIVATE math/gridding_kernel.cc
          math/gl_integrator.cc
          infra/threading.cc
          infra/mav.cc
          infra/string_utils.cc
          wgridder/wgridder.cc)

target_include_directories(${SCHAAPCOMMON_PROJECT_NAME}
                           PUBLIC ${DUCC0_PATH}/../)

if(SCHAAPCOMMON_MASTER_PROJECT)
  install(
    DIRECTORY ${DUCC0_PATH}
    DESTINATION "include"
    FILES_MATCHING
    PATTERN *.h)
endif()

