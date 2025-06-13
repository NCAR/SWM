
#TODO ADD OpenMP flags. 

set(SWM_AMREX_FSUBROUTINE_ACC_SOURCES
    ../main.cpp
    ../swm_mini_app_utils.cpp
    ../swm_mini_app_utils.h
    swm_fortran_kernelsOMP.F90
)

add_executable(swm_amrex_fsubroutine_acc ${SWM_AMREX_FSUBROUTINE_ACC_SOURCES})

target_link_libraries(swm_amrex_fsubroutine_acc PUBLIC AMReX::amrex_2d)
target_link_libraries(swm_amrex_fsubroutine_acc PUBLIC m)