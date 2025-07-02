
set(SWM_AMREX_TARGET swm_amrex)

if(SWM_OMP)
  set(SWM_AMREX_TARGET "${SWM_AMREX_TARGET}_OMP")
endif()

if(SWM_MPI)
  set(SWM_AMREX_TARGET "${SWM_AMREX_TARGET}_MPI")
endif()

if(SWM_CUDA)
  set(SWM_AMREX_TARGET "${SWM_AMREX_TARGET}_CUDA")
endif()

set(SWM_AMREX_SOURCES
    ${PROJECT_SOURCE_DIR}/swm_amrex/main.cpp
    ${PROJECT_SOURCE_DIR}/swm_amrex/swm_mini_app_utils.cpp
    swm_mini_app_utils_amrex.cpp
)

add_executable(${SWM_AMREX_TARGET} ${SWM_AMREX_SOURCES})

target_sources(${SWM_AMREX_TARGET} PRIVATE ${SWM_AMREX_SOURCES})

if (SWM_CUDA)
  include("AMReXTargetHelpers")
  setup_target_for_cuda_compilation(${SWM_AMREX_TARGET})
endif()

target_include_directories(${SWM_AMREX_TARGET} PUBLIC ${PROJECT_SOURCE_DIR}/swm_amrex)
target_include_directories(${SWM_AMREX_TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(${SWM_AMREX_TARGET} PUBLIC AMReX::amrex_2d)
target_link_libraries(${SWM_AMREX_TARGET} PUBLIC m)