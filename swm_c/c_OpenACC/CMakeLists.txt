# Option to select device: cpu or gpu (default: cpu)
set(ACC_DEVICE "cpu" CACHE STRING "Target device for OpenACC: cpu or gpu")
set_property(CACHE ACC_DEVICE PROPERTY STRINGS cpu gpu)

# Set OpenACC flags based on compiler and device
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    if(ACC_DEVICE STREQUAL "cpu")
        set(OPENACC_COMPILE_FLAGS -fopenacc -foffload=disable)
    elseif(ACC_DEVICE STREQUAL "gpu")
        set(OPENACC_COMPILE_FLAGS -fopenacc)
    endif()
    set(OPENACC_LINK_FLAGS -fopenacc)
elseif(CMAKE_C_COMPILER_ID STREQUAL "NVIDIA")
    if(ACC_DEVICE STREQUAL "cpu")
        set(OPENACC_FLAGS -acc=multicore -Minfo=accel -Mnofma)
        set(OPENACC_LINK_FLAGS -acc=multicore)
    elseif(ACC_DEVICE STREQUAL "gpu")
        set(OPENACC_FLAGS -acc=gpu -Minfo=accel -Mnofma )
        set(OPENACC_LINK_FLAGS -acc=gpu)
    endif()
else()
    # Maybe add a warning here that we are not using OpenACC for whatever compiler you are using.
    set(OPENACC_COMPILE_FLAGS "")
    set(OPENACC_LINK_FLAGS "")
endif()

# Maybe add back in the machine specific flags for the gpus we commonly run on.
# For example if we are running using the nvidia compiler compiler and derecho then add gpu=ccXX where XX is the compute capability of the gpu. 
# I am avoiding machine specific flags for now to keep this more portable for now.

###############################################################################
# C with OpenACC version
###############################################################################
set(SWM_C_ACC_SWAP_SOURCES
    shallow_swap.acc.c
    ../c/wtime.c
)

add_executable(swm_c_acc ${SWM_C_ACC_SWAP_SOURCES})
target_link_libraries(swm_c_acc m)
target_compile_options(swm_c_acc PRIVATE ${OPENACC_COMPILE_FLAGS})
target_link_options(swm_c_acc  PRIVATE ${OPENACC_LINK_FLAGS})

###############################################################################
# C with OpenACC version
###############################################################################
set(SWM_C_ACC_TILE_SOURCES
    shallow_swap.acc.Tile.c
    ../c/wtime.c
)

add_executable(swm_c_acc_tile ${SWM_C_ACC_TILE_SOURCES})
target_link_libraries(swm_c_acc_tile m)
target_compile_options(swm_c_acc PRIVATE ${OPENACC_COMPILE_FLAGS})
target_link_options(swm_c_acc  PRIVATE ${OPENACC_LINK_FLAGS})