# Fortran mock up of AMReX-style driver. But but still only using plain Fortran for CPU. The driver does not use AMReX library.
add_executable(swm_fortran_amrex_driver swm_fortran_amrex_driver.F90 swm_fortran_amrex_kernels.F90)
target_compile_options(swm_fortran_amrex_driver PRIVATE ${DEFAULT_REAL_FLAG})