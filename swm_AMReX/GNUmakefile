# AMREX_HOME defines the directory in which we will find all the AMReX code.
# I recommend adding this to your environment before running make. However if it is not defined then it will get set here:
AMREX_HOME ?= PROVIDE_YOUR_OWN_PATH_TO_AMREX_HERE

DEBUG        = FALSE
USE_MPI      = FALSE
USE_OMP      = FALSE
COMP         = gnu
DIM          = 2

CXXSTD       = c++20

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
