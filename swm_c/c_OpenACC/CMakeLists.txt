# Option to select device: cpu or gpu (default: cpu)
#set(ACC_DEVICE "cpu" CACHE STRING "Target device for OpenACC: cpu or gpu")
#set_property(CACHE ACC_DEVICE PROPERTY STRINGS cpu gpu)

foreach(device IN LISTS SWM_DEVICE)
    message(STATUS "Building for device: ${device}")
    # Add logic for each device here

    # Set OpenACC flags based on compiler and device
    if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
        if(${device} STREQUAL "cpu")
            set(OPENACC_COMPILE_FLAGS -fopenacc -foffload=disable)
        elseif(${device} STREQUAL "gpu")
            set(OPENACC_COMPILE_FLAGS -fopenacc)
        endif()
        set(OPENACC_LINK_FLAGS -fopenacc)
    elseif(CMAKE_C_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_C_COMPILER_ID STREQUAL "NVIDIA")
        if(${device} STREQUAL "cpu")
            set(OPENACC_FLAGS -acc=multicore -Minfo=accel -Mnofma)
            set(OPENACC_LINK_FLAGS -acc=multicore)
        elseif(${device} STREQUAL "gpu")
            set(OPENACC_FLAGS -acc=gpu -Minfo=accel -Mnofma )
            set(OPENACC_LINK_FLAGS -acc=gpu)
        endif()
    else()
        # Maybe add a warning here that we are not using OpenACC for whatever compiler you are using.
        message(WARNING "We have not added your support for your compiler (${CMAKE_C_COMPILER_ID}) for OpenACC build.")
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
    
    add_executable(swm_c_acc_${device} ${SWM_C_ACC_SWAP_SOURCES})
    target_link_libraries(swm_c_acc_${device} m)
    target_compile_options(swm_c_acc_${device} PRIVATE ${OPENACC_COMPILE_FLAGS})
    target_link_options(swm_c_acc_${device}  PRIVATE ${OPENACC_LINK_FLAGS})

    ###############################################################################
    # C with OpenACC and tile version
    ###############################################################################
    set(SWM_C_ACC_TILE_SOURCES
        shallow_swap.acc.Tile.c
        ../c/wtime.c
    )

    add_executable(swm_c_acc_tile_${device} ${SWM_C_ACC_TILE_SOURCES})
    target_link_libraries(swm_c_acc_tile_${device} m)
    target_compile_options(swm_c_acc_tile_${device} PRIVATE ${OPENACC_COMPILE_FLAGS})
    target_link_options(swm_c_acc_tile_${device}  PRIVATE ${OPENACC_LINK_FLAGS})

endforeach()