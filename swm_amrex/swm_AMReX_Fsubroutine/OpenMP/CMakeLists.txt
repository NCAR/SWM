################################################################################
## Hack to run old Makefile
################################################################################
#add_custom_target(swm_amrex_fsubroutine_omp ALL
#    COMMAND ${CMAKE_MAKE_PROGRAM} -f GNUmakefile
#    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#    COMMENT "Building with legacy GNUmakefile"
#)

###############################################################################
# Using CMake to buidd the SWM AMReX Fsubroutine with OpenACC support
###############################################################################
find_package(AMReX REQUIRED FORTRAN FINTERFACES) 

find_package(OpenMP REQUIRED)

# TODO: This is basicaly a copy of the block from the CMakeLists.txt for the C OpenACC build. Maybe combine in a function?
# Set OpenACC flags based on compiler and device
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVIDIA")
    if(SWM_DEVICE STREQUAL "gpu")
        set(OPENMP_COMPILE_FLAGS -mp=gpu)
        set(OPENMP_LINK_FLAGS -mp=gpu)
    endif()
else()
    # Maybe add a warning here that we are not using OpenACC for whatever compiler you are using.
    message(FATAL_ERROR "We have not added support for your compiler (${CMAKE_Fortran_COMPILER_ID}) for OpenACC build. Exiting configuration.")
endif()

set(SWM_AMREX_FSUBROUTINE_OMP_SOURCES
    ${PROJECT_SOURCE_DIR}/swm_amrex/common/main.cpp
    ${PROJECT_SOURCE_DIR}/swm_amrex/common/swm_mini_app_utils.cpp
    ../swm_mini_app_utils_fsubroutine.cpp
    swm_fortran_kernelsOMP.F90
)

add_executable(swm_amrex_fsubroutine_omp ${SWM_AMREX_FSUBROUTINE_OMP_SOURCES})


set_source_files_properties(swm_fortran_kernelsOMP.F90 PROPERTIES COMPILE_OPTIONS "${OPENMP_COMPILE_FLAGS};-Mnomain")

target_link_options(swm_amrex_fsubroutine_omp PRIVATE "${OPENMP_LINK_FLAGS};-Mnomain")

target_link_libraries(swm_amrex_fsubroutine_omp PRIVATE OpenMP::OpenMP_Fortran)

if (SWM_CUDA)
  include("AMReXTargetHelpers")
  setup_target_for_cuda_compilation(swm_amrex_fsubroutine_omp
    SOURCES ${SWM_AMREX_FSUBROUTINE_OMP_SOURCES}
  )
endif()

target_include_directories(swm_amrex_fsubroutine_omp PUBLIC ${PROJECT_SOURCE_DIR}/swm_amrex/common)
target_include_directories(swm_amrex_fsubroutine_omp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(swm_amrex_fsubroutine_omp PUBLIC AMReX::amrex_2d)
target_link_libraries(swm_amrex_fsubroutine_omp PUBLIC m)