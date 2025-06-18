#find_package(AMReX REQUIRED FORTRAN) 

#TODO ADD OpenMP flags. 

set(SWM_AMREX_FSUBROUTINE_OMP_SOURCES
    ${PROJECT_SOURCE_DIR}/swm_amrex/main.cpp
    ${PROJECT_SOURCE_DIR}/swm_amrex/swm_mini_app_utils.cpp
    ../swm_mini_app_utils_fsubroutine.cpp
    swm_fortran_kernelsOMP.F90
)

add_executable(swm_amrex_fsubroutine_omp ${SWM_AMREX_FSUBROUTINE_OMP_SOURCES})

target_include_directories(swm_amrex_fsubroutine_omp PUBLIC ${PROJECT_SOURCE_DIR}/swm_amrex)
target_include_directories(swm_amrex_fsubroutine_omp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(swm_amrex_fsubroutine_omp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(swm_amrex_fsubroutine_omp PUBLIC AMReX::amrex_2d)
target_link_libraries(swm_amrex_fsubroutine_omp PUBLIC m)