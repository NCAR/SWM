# Plain Fortran on CPU
add_executable(swm_fortran swm_fortran.F90 params.F90)
target_compile_options(swm_fortran PRIVATE ${DEFAULT_REAL_FLAG})
