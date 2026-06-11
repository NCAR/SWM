# Kokkos various policy implementations

set(KOKKOS_TARGETS
    shallow_rangepolicy
    shallow_rangepolicy_flat
    shallow_mdrangepolicy
    shallow_teampolicy
)

foreach(target ${KOKKOS_TARGETS})
    add_executable(${target} ${target}.cpp)
    target_link_libraries(${target} Kokkos::kokkos)
endforeach()
