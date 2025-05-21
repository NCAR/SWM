# AMREX_HOME = the directory where you pulled the AMReX repo.
# I recommend adding this to your environment before running make. However if it is not defined then it will get set here:
AMREX_HOME ?= PROVIDE_YOUR_OWN_PATH_TO_AMREX_HERE

DEBUG        = FALSE
USE_MPI      = FALSE
USE_OMP      = FALSE
USE_CUDA     = TRUE
USE_ACC         = FALSE
USE_OMP_OFFLOAD = TRUE
COMP         = nvhpc
DIM          = 2

TINY_PROFILE = TRUE
PROFILE = FALSE

LIBRARIES += -lstdc++ -lm

#USE_HDF5=TRUE
#HDF5_HOME=/home/lalo/spack/opt/spack/linux-ubuntu24.04-skylake/gcc-13.3.0/hdf5-1.14.5-vm6rifl2savzhynhyx56ayutfhmycs3s/

#CXXSTD       = c++20

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
