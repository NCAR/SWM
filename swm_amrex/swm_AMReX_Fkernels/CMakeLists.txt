set(SWM_AMREX_FKERNELS_SOURCES
    ${PROJECT_SOURCE_DIR}/swm_amrex/common/main.cpp
    ${PROJECT_SOURCE_DIR}/swm_amrex/common/swm_mini_app_utils.cpp
    swm_mini_app_utils_fkernels.cpp
    swm_fortran_amrex_kernels.F90
)

add_executable(swm_amrex_fkernels ${SWM_AMREX_FKERNELS_SOURCES})

target_include_directories(swm_amrex_fkernels PUBLIC ${PROJECT_SOURCE_DIR}/swm_amrex/common)
target_include_directories(swm_amrex_fkernels PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# This is a temporay fix to include the AMReX kernels from the swm_AMReX directory. This version should not be using that. Instead it should use the fortran kernels in this directory.
target_include_directories(swm_amrex_fkernels PUBLIC ${PROJECT_SOURCE_DIR}/swm_amrex/swm_AMReX)

target_link_libraries(swm_amrex_fkernels PUBLIC AMReX::amrex_2d)
target_link_libraries(swm_amrex_fkernels PUBLIC m)
