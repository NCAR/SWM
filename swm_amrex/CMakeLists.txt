#find_package(AMReX)

#enable_language(CUDA) # <-- Add this line at the top if not already present
find_package(AMReX REQUIRED CUDA)

add_subdirectory(swm_AMReX)
#add_subdirectory(swm_AMReX_Fkernels)
#add_subdirectory(swm_AMReX_Fsubroutine)