# Define the compiler
FC = gfortran

# Common source files used by all executables
COMMON_SRCS =

FFLAGS = -O2 -fdefault-real-8# -cpp
#FFLAGS = -Wall -g

GPUFFLAGS += gpu=cc80 -Minfo=accel -Mnofma

# Libraries to link to
LDLIBS = -lm

# Default target
#all: swm_fortran swm_fortran_driver
all:nvfortran_gpu nvfortran

nvfortran_gpu:
	touch swm_fortran.F90
	make FC=nvfortran FFLAGS="-O2 -r8 -acc=gpu" swm_fortran_gpu

nvfortran:
	touch swm_fortran.F90
	make FC=nvfortran FFLAGS="-O2 -r8 -noacc" swm_fortran_cpu

gfortran_gpu:
	touch swm_fortran.F90
	make FC=gfortran FFLAGS="-O2 -fdefault-real-8 -fopenacc -foffload=default" swm_fortran_gpu

swm_fortran_cpu: swm_fortran.F90 params.o
	$(FC) $(FFLAGS) -o $@ params.o swm_fortran.F90

swm_fortran_noacc: swm_fortran.F90
	$(FC) $(FFLAGS) -o $@ swm_fortran.F90

swm_fortran_gpu: swm_fortran_driver.F90 swm_fortran_kernels.o params.o
	$(FC) $(FFLAGS) -o $@ swm_fortran_kernels.o params.o swm_fortran_driver.F90

swm_fortran_kernels.o: swm_fortran_kernels.F90
	$(FC) $(FFLAGS) -c $@ swm_fortran_kernels.F90

params.o: params.F90
	$(FC) $(FFLAGS) -c $@ params.F90

# Modularized with AMREX-styled kernels (pass 3D array inside loop)
swm_fortran_amrex_driver: swm_fortran_amrex_driver.F90 swm_fortran_amrex_kernels.o
	$(FC) $(FFLAGS) -o $@ swm_fortran_amrex_kernels.o swm_fortran_amrex_driver.F90

swm_fortran_amrex_kernels.o: swm_fortran_amrex_kernels.F90
	$(FC) $(FFLAGS) -c $@ swm_fortran_amrex_kernels.F90

clean:
	rm -f swm_fortran swm_fortran_noacc *.o *.mod
