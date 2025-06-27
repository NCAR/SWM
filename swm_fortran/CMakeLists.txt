
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

if (SWM_DEVICE STREQUAL "cpu")
    add_subdirectory(fortran)
    add_subdirectory(fortran_AMReX_proxy)
endif()

if (SWM_OPENACC)
    add_subdirectory(fortran_OpenACC)
endif()
