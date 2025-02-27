# Define the compiler
FC = gfortran

# Common source files used by all executables
COMMON_SRCS =

FFLAGS = -O2 -fdefault-real-8# -cpp
#FFLAGS = -Wall -g

ifdef GPU
    FFLAGS += gpu=cc70 -Minfo=accel -Mnofma
endif

# Libraries to link to
LDLIBS = -lm

# Default target
all: swm_fortran swm_fortran_driver

nvfortran:
	make FC=nvfortran FFLAGS="-r8 -acc=gpu" swm_fortran swm_fortran_driver

nvfortran_noacc:
	make FC=nvfortran FFLAGS="-r8 -noacc" swm_fortran swm_fortran_driver

swm_fortran: swm_fortran.F90
	$(FC) $(FFLAGS) -o $@ swm_fortran.F90

swm_fortran_driver: swm_fortran_driver.F90 swm_fortran_kernels.o
	$(FC) $(FFLAGS) -o $@ swm_fortran_kernels.o swm_fortran_driver.F90

swm_fortran_kernels.o: swm_fortran_kernels.F90
	$(FC) $(FFLAGS) -c $@ swm_fortran_kernels.F90

# Modularized with AMREX-styled kernels (pass 3D array inside loop)
swm_fortran_amrex_driver: swm_fortran_amrex_driver.F90 swm_fortran_amrex_kernels.o
	$(FC) $(FFLAGS) -o $@ swm_fortran_amrex_kernels.o swm_fortran_amrex_driver.F90

swm_fortran_amrex_kernels.o: swm_fortran_amrex_kernels.F90
	$(FC) $(FFLAGS) -c $@ swm_fortran_amrex_kernels.F90

clean:
	rm -f swm_fortran swm_fortran_driver swm_fortran_amrex_driver *.o *.mod
