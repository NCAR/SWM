
set(AMREX_COMPONENTS 2D DOUBLE)

if(SWM_CUDA)
  list(APPEND AMREX_COMPONENTS CUDA)

  enable_language(CUDA)
  if(CMAKE_VERSION VERSION_LESS 3.20)
     include(AMReX_SetupCUDA)
  endif()

endif()

if(SWM_MPI)
  list(APPEND AMREX_COMPONENTS MPI)
endif()

if(SWM_OPENMP)
  list(APPEND AMREX_COMPONENTS OMP)
endif()

find_package(AMReX REQUIRED ${AMREX_COMPONENTS})

add_subdirectory(swm_AMReX)
#add_subdirectory(swm_AMReX_Fkernels)
add_subdirectory(swm_AMReX_Fsubroutine)