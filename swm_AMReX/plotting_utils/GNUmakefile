# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME ?= PROVIDE_YOUR_OWN_PATH_TO_AMREX_HERE
HDF5_HOME  ?= PROVIDE_YOUR_OWN_PATH_TO_HDF5_HERE

DEBUG        = FALSE
USE_MPI      = TRUE   # Dont really need MPI for the program but I am using an hdf5 that needs MPI so I am setting this to TRUE
USE_OMP      = FALSE
COMP         = gnu
DIM          = 2

USE_HDF5     = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
