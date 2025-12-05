cmake_minimum_required(VERSION 3.13)

project(SWM C CXX Fortran)

# Assume we want to build in Release mode by default unless specified otherwise
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(SWM_DEVICE "cpu" CACHE STRING "Select device: cpu or gpu")
set_property(CACHE SWM_DEVICE PROPERTY STRINGS cpu gpu)

option(SWM_C "Enable the C versions of the mini-app" ON)
option(SWM_FORTRAN "Enable the Fortran versions of the mini-app" ON)
option(SWM_AMREX "Enable the AMReX version of the mini-app" OFF)
option(SWM_KOKKOS "Enable the Kokkos version of the mini-app" OFF)

option(SWM_OPENACC "Enable the OpenACC versions of the mini-app" OFF)
option(SWM_OPENMP "Enable the OpenMP versions of the mini-app" OFF)
option(SWM_MPI "Enable the MPI versions of the mini-app" OFF)
option(SWM_CUDA "Enable the CUDA versions of the mini-app" OFF)
option(SWM_HIP "Enable the HIP versions of the mini-app" OFF)

if(SWM_CUDA AND SWM_HIP)
    message(FATAL_ERROR "SWM_CUDA and SWM_HIP cannot both be ON. Please enable only one.")
endif()

if(SWM_DEVICE STREQUAL "cpu" AND SWM_CUDA)
    message(WARNING "SWM_CUDA is ON but SWM_DEVICE is set to 'cpu'. Disabling SWM_CUDA.")
    set(SWM_CUDA OFF CACHE BOOL "Enable the CUDA versions of the mini-app" FORCE)
endif()

if(SWM_DEVICE STREQUAL "cpu" AND SWM_HIP)
    message(WARNING "SWM_HIP is ON but SWM_DEVICE is set to 'cpu'. Disabling SWM_HIP.")
    set(SWM_HIP OFF CACHE BOOL "Enable the HIP versions of the mini-app" FORCE)
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

if(SWM_KOKKOS)
    add_subdirectory(swm_kokkos)
endif()