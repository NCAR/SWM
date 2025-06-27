# Set OpenACC flags based on compiler and device
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    if(SWM_DEVICE STREQUAL "cpu")
        set(OPENACC_COMPILE_FLAGS -fopenacc -foffload=disable)
    elseif(SWM_DEVICE STREQUAL "gpu")
        set(OPENACC_COMPILE_FLAGS -fopenacc)
    endif()
    set(OPENACC_LINK_FLAGS -fopenacc)
elseif(CMAKE_C_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_C_COMPILER_ID STREQUAL "NVIDIA")
    if(SWM_DEVICE STREQUAL "cpu")
        set(OPENACC_FLAGS -acc=multicore -Minfo=accel -Mnofma)
        set(OPENACC_LINK_FLAGS -acc=multicore)
    elseif(SWM_DEVICE STREQUAL "gpu")
        set(OPENACC_FLAGS -acc=gpu -Minfo=accel -Mnofma )
        set(OPENACC_LINK_FLAGS -acc=gpu)
    endif()
else()
    message(WARNING "We have not added your support for your compiler (${CMAKE_C_COMPILER_ID}) for OpenACC build. So we are not passing any OpenACC flags to the compiler. It is likely the pragmas will be ignored.")
    set(OPENACC_COMPILE_FLAGS "")
    set(OPENACC_LINK_FLAGS "")
endif()


###############################################################################
# C with OpenACC version
###############################################################################
set(SWM_C_ACC_SWAP_SOURCES
    shallow_swap.acc.c
    ../c/wtime.c
)

#set(SWM_C_ACC_SWAP_TARGET swm_c_acc_${SWM_DEVICE})
set(SWM_C_ACC_SWAP_TARGET swm_c_acc)

add_executable(${SWM_C_ACC_SWAP_TARGET} ${SWM_C_ACC_SWAP_SOURCES})
target_link_libraries(${SWM_C_ACC_SWAP_TARGET} m)
target_compile_options(${SWM_C_ACC_SWAP_TARGET} PRIVATE ${OPENACC_COMPILE_FLAGS})
target_link_options(${SWM_C_ACC_SWAP_TARGET} PRIVATE ${OPENACC_LINK_FLAGS})

###############################################################################
# C with OpenACC and tile version
###############################################################################
set(SWM_C_ACC_TILE_SOURCES
    shallow_swap.acc.Tile.c
    ../c/wtime.c
)

#set(SWM_C_ACC_TILE_TARGET swm_c_acc_tile_${SWM_DEVICE})
set(SWM_C_ACC_TILE_TARGET swm_c_acc_tile)

add_executable(${SWM_C_ACC_TILE_TARGET} ${SWM_C_ACC_TILE_SOURCES})
target_link_libraries(${SWM_C_ACC_TILE_TARGET} m)
target_compile_options(${SWM_C_ACC_TILE_TARGET} PRIVATE ${OPENACC_COMPILE_FLAGS})
target_link_options(${SWM_C_ACC_TILE_TARGET} PRIVATE ${OPENACC_LINK_FLAGS})