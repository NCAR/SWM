set(SWM_C_SOURCES
    shallow_swap.c
    ../common/wtime.c
)

add_executable(swm_c ${SWM_C_SOURCES})

target_compile_definitions(swm_c PRIVATE M=${M})
target_compile_definitions(swm_c PRIVATE N=${N})

# Link math library
target_link_libraries(swm_c m)