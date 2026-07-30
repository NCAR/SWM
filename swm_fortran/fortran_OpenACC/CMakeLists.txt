# TODO: This is basicaly a copy of the block from the CMakeLists.txt for the C OpenACC build. Maybe combine in a function?
# Set OpenACC flags based on compiler and device
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if(SWM_DEVICE STREQUAL "cpu")
        set(OPENACC_COMPILE_FLAGS -fopenacc -foffload=disable)
    elseif(SWM_DEVICE STREQUAL "gpu")
        set(OPENACC_COMPILE_FLAGS -fopenacc)
    endif()
    set(OPENACC_LINK_FLAGS -fopenacc)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVIDIA")
    if(SWM_DEVICE STREQUAL "cpu")
        set(OPENACC_COMPILE_FLAGS -acc=multicore -Minfo=accel -Mnofma)
        set(OPENACC_LINK_FLAGS -acc=multicore)
    elseif(SWM_DEVICE STREQUAL "gpu")
        set(OPENACC_COMPILE_FLAGS -acc=gpu -Minfo=accel -Mnofma)
        set(OPENACC_LINK_FLAGS -acc=gpu)
    endif()
else()
    # Maybe add a warning here that we are not using OpenACC for whatever compiler you are using.
    message(WARNING "We have not added your support for your compiler (${CMAKE_Fortran_COMPILER_ID}) for OpenACC build. So we are not passing any OpenACC flags to the compiler. It is likely the pragmas will be ignored.")
    set(OPENACC_COMPILE_FLAGS "")
    set(OPENACC_LINK_FLAGS "")
endif()

set(SWM_FORTRAN_ACC_TARGET swm_fortran_acc)

add_executable(${SWM_FORTRAN_ACC_TARGET} swm_fortran_driver.F90 swm_fortran_kernels.F90 ../common/params.F90)
target_compile_options(${SWM_FORTRAN_ACC_TARGET} PRIVATE ${DEFAULT_REAL_FLAG} ${OPENACC_COMPILE_FLAGS})
target_link_options(${SWM_FORTRAN_ACC_TARGET} PRIVATE ${OPENACC_LINK_FLAGS})