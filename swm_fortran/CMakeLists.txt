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

# Set OpenACC flags based on compiler
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(OPENACC_FLAGS -fopenacc -foffload=default)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(OPENACC_FLAGS -acc=gpu)
else()
    set(OPENACC_FLAGS "")
endif()

# Plain Fortran CPU
add_executable(swm_fortran swm_fortran.F90 params.F90)
target_compile_options(swm_fortran PRIVATE ${DEFAULT_REAL_FLAG})

# Fortran with OpenACC on CPU 
add_executable(swm_fortran_acc swm_fortran_driver.F90 swm_fortran_kernels.F90 params.F90)
target_compile_options(swm_fortran_acc PRIVATE ${DEFAULT_REAL_FLAG})

# Fortran with OpenACC on GPU target
add_executable(swm_fortran_acc_gpu swm_fortran_driver.F90 swm_fortran_kernels.F90 params.F90)
target_compile_options(swm_fortran_acc_gpu PRIVATE ${DEFAULT_REAL_FLAG} ${OPENACC_FLAGS})
target_link_options(swm_fortran_acc_gpu PRIVATE ${OPENACC_FLAGS})

# Fortran mock up of AMReX-style driver
add_executable(swm_fortran_amrex_driver swm_fortran_amrex_driver.F90 swm_fortran_amrex_kernels.F90)
target_compile_options(swm_fortran_amrex_driver PRIVATE ${DEFAULT_REAL_FLAG})