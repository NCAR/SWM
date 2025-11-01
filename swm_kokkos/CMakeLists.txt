set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)

find_package(Kokkos REQUIRED)

# Use -O2 optimization for release builds to be consistent with the C implementation
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

# Flags that get used when building debug flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g")

set(SWM_KOKKOS_SOURCES shallow_teampolicy.cpp)

add_executable(swm_kokkos ${SWM_KOKKOS_SOURCES})

target_link_libraries(swm_kokkos Kokkos::kokkos)