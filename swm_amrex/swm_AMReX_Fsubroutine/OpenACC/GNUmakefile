# AMREX_HOME = the directory where you pulled the AMReX repo.
# I recommend adding this to your environment before running make. However if it is not defined then it will get set here:
AMREX_HOME ?= PROVIDE_YOUR_OWN_PATH_TO_AMREX_HERE

DEBUG        = FALSE
USE_MPI      = FALSE
USE_OMP      = FALSE
USE_CUDA     = TRUE
USE_ACC         = TRUE
USE_OMP_OFFLOAD = FALSE
COMP         = nvhpc
DIM          = 2

TINY_PROFILE = TRUE
PROFILE = FALSE

LIBRARIES += -lstdc++ -lm

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
