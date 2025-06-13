find_package(AMReX REQUIRED FORTRAN) 

# ./funcF.H gets included in the mini app utils.cpp file here. Maybe build the utils object at this level and then just link to it in the lower level directories?

add_subdirectory(OpenACC)
add_subdirectory(OpenMP)