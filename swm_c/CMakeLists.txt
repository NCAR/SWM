set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Add optimization flag
#add_compile_options(-O2)
# Maybe move this to be set if CMAKE_BUILD_TYPE is Release
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2")

# Flags that get used when building  debug flags
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g")

# Additional flags if using the NVIDIA compiler for OpenACC
set(NVC_OPENACC_FLAGS gpu=cc70 -Minfo=accel -Mnofma)

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

add_executable(swm_c_acc ${SWM_C_ACC_SWAP_SOURCES})

target_link_libraries(swm_c_acc m)

if(CMAKE_C_COMPILER_ID STREQUAL "NVIDIA")
    target_compile_options(swm_c_acc PRIVATE ${NVC_OPENACC_FLAGS})
endif()

###############################################################################
# C with OpenACC version
###############################################################################
set(SWM_C_ACC_TILE_SOURCES
    shallow_swap.acc.Tile.c
    wtime.c
)

add_executable(swm_c_acc_tile ${SWM_C_ACC_TILE_SOURCES})
target_link_libraries(swm_c_acc_tile m)

if(CMAKE_C_COMPILER_ID STREQUAL "NVIDIA")
    target_compile_options(swm_c_acc_tile PRIVATE ${NVC_OPENACC_FLAGS})
endif()
