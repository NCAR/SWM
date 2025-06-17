#find_package(AMReX REQUIRED FORTRAN) 

#TODO ADD OpenACC flags. 

set(SWM_AMREX_FSUBROUTINE_ACC_SOURCES
    ${PROJECT_SOURCE_DIR}/swm_amrex/main.cpp
    ../swm_mini_app_utils.cpp
    swm_fortran_kernelsACC.F90
)

add_executable(swm_amrex_fsubroutine_acc ${SWM_AMREX_FSUBROUTINE_ACC_SOURCES})

target_include_directories(swm_amrex_fsubroutine_acc PUBLIC ${PROJECT_SOURCE_DIR}/swm_amrex ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(swm_amrex_fsubroutine_acc PUBLIC AMReX::amrex_2d)
target_link_libraries(swm_amrex_fsubroutine_acc PUBLIC m)