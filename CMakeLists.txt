cmake_minimum_required(VERSION 3.13)

project(SWM C CXX Fortran)

# Assume we want to build in Release mode by default unless specified otherwise
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_subdirectory(swm_c)
add_subdirectory(swm_fortran)
add_subdirectory(swm_amrex)
