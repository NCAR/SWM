if (SWM_CUDA)
  enable_language(CUDA)
  if(CMAKE_VERSION VERSION_LESS 3.20)
     include(AMReX_SetupCUDA)
  endif()
endif()

set(SWM_AMREX_SOURCES
    ${PROJECT_SOURCE_DIR}/swm_amrex/main.cpp
    ${PROJECT_SOURCE_DIR}/swm_amrex/swm_mini_app_utils.cpp
    swm_mini_app_utils_amrex.cpp
)

add_executable(swm_amrex ${SWM_AMREX_SOURCES})

target_sources(swm_amrex PRIVATE ${SWM_AMREX_SOURCES})

if (SWM_CUDA)
  include("AMReXTargetHelpers")
  setup_target_for_cuda_compilation(swm_amrex)
endif()

target_include_directories(swm_amrex PUBLIC ${PROJECT_SOURCE_DIR}/swm_amrex)
target_include_directories(swm_amrex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(swm_amrex PUBLIC AMReX::amrex_2d)
target_link_libraries(swm_amrex PUBLIC m)
