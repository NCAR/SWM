# Plain Fortran on CPU
add_executable(swm_fortran swm_fortran.F90 ../common/params.F90)
target_compile_options(swm_fortran PRIVATE ${DEFAULT_REAL_FLAG})

target_compile_definitions(swm_fortran PRIVATE MNUM=${M})
target_compile_definitions(swm_fortran PRIVATE NNUM=${N})