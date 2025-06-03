set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Add optimization flag
add_compile_options(-O2)

###############################################################################
# Plain C version
###############################################################################
set(SWM_C_SOURCES
    shallow_swap.c
    wtime.c
)

add_executable(swm_c ${SWM_C_SOURCES})

# Link math library
target_link_libraries(swm_c m)

###############################################################################
# C with OpenACC version
###############################################################################
set(SWM_C_ACC_SWAP_SOURCES
    shallow_swap.acc.c
    wtime.c
)

# TODO: Add OpenACC support if available... via compiler flags or CMake options

add_executable(swm_c_acc ${SWM_C_ACC_SWAP_SOURCES})

target_link_libraries(swm_c_acc m)

###############################################################################
# C with OpenACC version
###############################################################################
set(SWM_C_ACC_TILE_SOURCES
    shallow_swap.acc.Tile.c
    wtime.c
)

# TODO: Add OpenACC support if available... via compiler flags or CMake options

add_executable(swm_c_acc_tile ${SWM_C_ACC_TILE_SOURCES})

target_link_libraries(swm_c_acc_tile m)
