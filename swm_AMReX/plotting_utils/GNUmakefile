# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME ?= PROVIDE_YOUR_OWN_PATH_TO_AMREX_HERE
HDF5_HOME  ?= PROVIDE_YOUR_OWN_PATH_TO_HDF5_HERE

DEBUG        = FALSE
USE_MPI      = TRUE   # Dont really nedd MPI for the program but I am using an hdf5 that needs MPI
USE_OMP      = FALSE
COMP         = gnu
DIM          = 2

USE_HDF5     = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

## Set the target name explicitly
#TARGET = plotfile2hdf5
## Set the target name explicitly
#EXEC = plotfile2hdf5
#
#executable = plotfile2hdf5.exe

include $(AMREX_HOME)/Src/Base/Make.package

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
