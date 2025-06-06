#enable_language(Fortran)

# Set default flags
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
set(CMAKE_Fortran_FLAGS_DEBUG   "-Wall -g")

# Set default real kind flag based on compiler
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(DEFAULT_REAL_FLAG "-fdefault-real-8")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(DEFAULT_REAL_FLAG "-r8")
else()
    set(DEFAULT_REAL_FLAG "")
endif()

add_subdirectory(fortran)
add_subdirectory(fortran_OpenACC)
add_subdirectory(fortran_AMReX_proxy)


## Fortran mock up of AMReX-style driver
#add_executable(swm_fortran_amrex_driver swm_fortran_amrex_driver.F90 swm_fortran_amrex_kernels.F90)
#target_compile_options(swm_fortran_amrex_driver PRIVATE ${DEFAULT_REAL_FLAG})