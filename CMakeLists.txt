# Kokkos minimally requires 3.10 right now,
# but your project can set it higher
cmake_minimum_required(VERSION 3.10)

# Declare your project name and list (white space separated)
# the file types of the source code
project(SWM CXX)

# You need this for using Kokkos_ROOT variable
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.12.0")
  message(STATUS "Setting policy CMP0074 to use <Package>_ROOT variables")
  cmake_policy(SET CMP0074 NEW)
endif()

# Look for an installed Kokkos (using Kokkos_DIR environment variable)
find_package(Kokkos REQUIRED)

# Declare executable name and list files utilized in compilation
add_executable(swm_kokkos shallow_unroll.kokkos.cpp wtime.cpp)

# This is the only thing required to set up compiler/linker flags
target_link_libraries(swm_kokkos Kokkos::kokkos)