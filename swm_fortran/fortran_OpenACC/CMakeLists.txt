# Option to select device: cpu or gpu (default: cpu)
set(ACC_DEVICE "cpu" CACHE STRING "Target device for OpenACC: cpu or gpu")
set_property(CACHE ACC_DEVICE PROPERTY STRINGS cpu gpu)

# This is basically a copy of the CMakeLists.txt file in swm_c/c_OpenACC/CMakeLists.txt just adapted for Fortran. It sets the OpenACC flags based on the compiler and device type. Maybe just do this in a common file later.
# Set OpenACC flags based on compiler and device
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if(ACC_DEVICE STREQUAL "cpu")
        set(OPENACC_COMPILE_FLAGS -fopenacc -foffload=disable)
    elseif(ACC_DEVICE STREQUAL "gpu")
        set(OPENACC_COMPILE_FLAGS -fopenacc)
    endif()
    set(OPENACC_LINK_FLAGS -fopenacc)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    if(ACC_DEVICE STREQUAL "cpu")
        set(OPENACC_COMPILE_FLAGS -acc=multicore -Minfo=accel -Mnofma)
        set(OPENACC_LINK_FLAGS -acc=multicore)
    elseif(ACC_DEVICE STREQUAL "gpu")
        set(OPENACC_COMPILE_FLAGS -acc=gpu -Minfo=accel -Mnofma)
        set(OPENACC_LINK_FLAGS -acc=gpu)
    endif()
else()
    # Maybe add a warning here that we are not using OpenACC for whatever compiler you are using.
    set(OPENACC_COMPILE_FLAGS "")
    set(OPENACC_LINK_FLAGS "")
endif()

# Fortran with OpenACC on CPU 
add_executable(swm_fortran_acc swm_fortran_driver.F90 swm_fortran_kernels.F90 ../fortran/params.F90)
target_compile_options(swm_fortran_acc PRIVATE ${DEFAULT_REAL_FLAG} ${OPENACC_COMPILE_FLAGS})
target_link_options(swm_fortran_acc PRIVATE ${OPENACC_LINK_FLAGS})