###############################################################################
# Hack to run old Makefile
###############################################################################
add_custom_target(run_gnumakefile ALL
    COMMAND ${CMAKE_MAKE_PROGRAM} -f GNUmakefile
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Building with legacy GNUmakefile"
)

###############################################################################
# Using CMake to buidd the SWM AMReX Fsubroutine with OpenACC support
###############################################################################
##find_package(AMReX REQUIRED FORTRAN) 
#
##TODO ADD OpenMP flags. 
#
#set(SWM_AMREX_FSUBROUTINE_OMP_SOURCES
#    ${PROJECT_SOURCE_DIR}/swm_amrex/common/main.cpp
#    ${PROJECT_SOURCE_DIR}/swm_amrex/common/swm_mini_app_utils.cpp
#    ../swm_mini_app_utils_fsubroutine.cpp
#    swm_fortran_kernelsOMP.F90
#)
#
#add_executable(swm_amrex_fsubroutine_omp ${SWM_AMREX_FSUBROUTINE_OMP_SOURCES})
#
#target_include_directories(swm_amrex_fsubroutine_omp PUBLIC ${PROJECT_SOURCE_DIR}/swm_amrex/common)
#target_include_directories(swm_amrex_fsubroutine_omp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#target_include_directories(swm_amrex_fsubroutine_omp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
#
#target_link_libraries(swm_amrex_fsubroutine_omp PUBLIC AMReX::amrex_2d)
#target_link_libraries(swm_amrex_fsubroutine_omp PUBLIC m)