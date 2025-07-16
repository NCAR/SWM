find_package(AMReX REQUIRED FORTRAN FINTERFACES) 

# These version are only for offloading onto the GPU with CUDA backend
if(SWM_DEVICE STREQUAL "gpu" AND SWM_CUDA)

  # For now only build the NVIDIA compiler... could update later to support other compilers if we test this.
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVIDIA")

    if(SWM_OPENACC)
        add_subdirectory(OpenACC)
    endif()

    if(SWM_OPENMP)
        add_subdirectory(OpenMP)
    endif()

  endif()

endif()

