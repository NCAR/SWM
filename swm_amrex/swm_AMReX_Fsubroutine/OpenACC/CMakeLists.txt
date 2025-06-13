#find_package(AMReX REQUIRED FORTRAN) 

#TODO ADD OpenACC flags. 

set(SWM_AMREX_FSUBROUTINE_OMP_SOURCES
    ../main.cpp
    ../swm_mini_app_utils.cpp
    ../swm_mini_app_utils.h
    swm_fortran_kernelsACC.F90
)
#    ../funcF.H does not get used in this version?

add_executable(swm_amrex_fsubroutine_omp ${SWM_AMREX_FSUBROUTINE_OMP_SOURCES})

target_link_libraries(swm_amrex_fsubroutine_omp PUBLIC AMReX::amrex_2d)
target_link_libraries(swm_amrex_fsubroutine_omp PUBLIC m)