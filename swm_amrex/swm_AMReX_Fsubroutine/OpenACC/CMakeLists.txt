#find_package(AMReX REQUIRED FORTRAN) 

#TODO ADD OpenMP flags. 

set(SWM_AMREX_FSUBROUTINE_ACC_SOURCES
    ../../main.cpp
    ../swm_mini_app_utils.cpp
    ../swm_mini_app_utils.h
    swm_fortran_kernelsACC.F90
)
#    ../funcF.H does not get used in this version?

add_executable(swm_amrex_fsubroutine_acc ${SWM_AMREX_FSUBROUTINE_ACC_SOURCES})

target_include_directories(swm_amrex_fsubroutine_acc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_include_directories(swm_amrex_fsubroutine_acc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(swm_amrex_fsubroutine_acc PUBLIC AMReX::amrex_2d)
target_link_libraries(swm_amrex_fsubroutine_acc PUBLIC m)
