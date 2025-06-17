set(SWM_AMREX_FKERNELS_SOURCES
    ../main.cpp
    swm_mini_app_utils.cpp
    swm_fortran_amrex_kernels.F90
)

add_executable(swm_amrex_fkernels ${SWM_AMREX_FKERNELS_SOURCES})

target_include_directories(swm_amrex_fkernels PUBLIC ${PROJECT_SOURCE_DIR}/swm_amrex ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(swm_amrex_fkernels PUBLIC AMReX::amrex_2d)
target_link_libraries(swm_amrex_fkernels PUBLIC m)
