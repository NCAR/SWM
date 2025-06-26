cmake_minimum_required(VERSION 3.13)

project(SWM C CXX Fortran)

# Assume we want to build in Release mode by default unless specified otherwise
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()


set(SWM_DEVICE "cpu" CACHE STRING "Select device(s): cpu, gpu, or cpu;gpu for both (semicolon-separated)")
set_property(CACHE SWM_DEVICE PROPERTY STRINGS cpu gpu "cpu;gpu")

option(SWM_C "Enable the C versions of the mini-app" ON)
option(SWM_FORTRAN "Enable the Fortran versions of the mini-app" ON)
option(SWM_AMREX "Enable the AMReX version of the mini-app" OFF)

option(SWM_OPENACC "Enable the OpenACC versions of the mini-app" ON)
option(SWM_OPENMP "Enable the OpenMP versions of the mini-app" OFF)
option(SWM_CUDA "Enable the CUDA versions of the mini-app" OFF)


if("cpu" IN_LIST SWM_DEVICE)
    # Enable CPU target
endif()

if("gpu" IN_LIST SWM_DEVICE)
    # Enable GPU target
endif()

if (SWM_C)
  add_subdirectory(swm_c)
endif()

if (SWM_FORTRAN)
  add_subdirectory(swm_fortran)
endif() 

if(SWM_AMREX)
    add_subdirectory(swm_amrex)
endif()
