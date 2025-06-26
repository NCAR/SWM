
if (SWM_CUDA)
  find_package(AMReX REQUIRED CUDA)
else()
  find_package(AMReX REQUIRED)
endif()

add_subdirectory(swm_AMReX)
#add_subdirectory(swm_AMReX_Fkernels)
#add_subdirectory(swm_AMReX_Fsubroutine)