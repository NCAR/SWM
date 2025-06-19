#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status
set -u  # Treat expanding empty variables as an error
set -x  # Print each command before executing it

###############################################################################
# User Input
###############################################################################


# Directory where you pulled the AMREX repositry
export AMREX_HOME=/glade/u/home/htorres/amrex

# Directory where you pulled the SWM repository
export SWM_ROOT=/glade/u/home/htorres/SWM

# Directory where we are going to create the SWM build
export SWM_BUILD_DIR=/glade/u/home/htorres/SWM_build

# Options... example to pick cpu vs gpu?

###############################################################################
# Module Setup
###############################################################################

module purge

# Modules we always use
module load cmake

## Serial
module load gcc

# MPI Only
#module load gcc cray-mpich

## MPI with w/ nvhpc so we can profile the code with nsys
#module load nvhpc cray-mpich

## CUDA Only
#module load nvhpc cuda

## MPI and CUDA
#module load nvhpc cuda cray-mpich

# MPI and CUDA on multi-gpu
#module load gcc cuda cray-mpich ncarcompilers

# HDF5
module load hdf5
# Loading the hdf5 module on derecheo should define a variable NCAR_ROOT_HDF5. We need to set HDF5_HOME
export HDF5_HOME="${NCAR_ROOT_HDF5}"

# Profile or Debugg with linary forge
#module load linaro-forge

## Python for plotting
#module load conda
##conda activate npl
#conda activate swm_amrex # A conda environment with numpy matplotlib and h5py

module list

###############################################################################
# Build the version of amrex that we are asking for
###############################################################################

export AMREX_BUILD_DIR=/glade/u/home/htorres/amrex_build
export AMREX_INSTALL_DIR=$AMREX_BUILD_DIR/install

# Initialize an array for CMake options
amrex_cmake_opts=()

# Always-used options
amrex_cmake_opts+=("-DAMReX_SPACEDIM=2")
amrex_cmake_opts+=("-DAMReX_PRECISION=DOUBLE")
amrex_cmake_opts+=("-DAMReX_FORTRAN=YES")
# These were on by default but I don't think we are using any of these features so I am turning them off
amrex_cmake_opts+=("-DAMReX_GPU_RDC=NO")
amrex_cmake_opts+=("-DAMReX_LINEAR_SOLVERS=NO")
amrex_cmake_opts+=("-DAMReX_LINEAR_SOLVERS_INCFLO=NO")
amrex_cmake_opts+=("-DAMReX_LINEAR_SOLVERS_EM=NO")
amrex_cmake_opts+=("-DAMReX_AMRLEVEL=NO")
amrex_cmake_opts+=("-DAMReX_PARTICLES=NO")
amrex_cmake_opts+=("-DAMReX_TINY_PROFILE=NO")

# These take on differnent values depending on the options chosen
amrex_cmake_opts+=("-DAMReX_MPI=NO")
amrex_cmake_opts+=("-DAMReX_OMP=NO")
amrex_cmake_opts+=("-DAMReX_GPU_BACKEND=NONE")

#amrex_cmake_opts+=("-D=")

## Example: add options based on logic
#if [[ "$USE_GPU" == "YES" ]]; then
#    amrex_cmake_opts+=("-DAMReX_GPU=YES")
#fi

# Add more options as needed...

mkdir -p $AMREX_BUILD_DIR
cd    $AMREX_BUILD_DIR

cmake "${amrex_cmake_opts[@]}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$AMREX_INSTALL_DIR" \
    -S "$AMREX_HOME" \
    -B "$AMREX_BUILD_DIR"

make  install
make  test_install  # optional step to test if the installation is working

###############################################################################
# Build SWM Using the version of AMReX that we just built
###############################################################################
cmake -DAMReX_ROOT=$AMREX_INSTALL_DIR/lib/cmake/AMReX -S $SWM_ROOT -B $SWM_BUILD_DIR

cd $SWM_BUILD_DIR
make

###############################################################################
# Run all the versions of SWM mini-app
###############################################################################
$SWM_BUILD_DIR/swm_c/c/swm_c
$SWM_BUILD_DIR/swm_c/c_OpenACC/swm_c_acc
$SWM_BUILD_DIR/swm_c/c_OpenACC/swm_c_acc_tile
$SWM_BUILD_DIR/swm_fortran/fortran/swm_fortran
$SWM_BUILD_DIR/swm_fortran/fortran_AMReX_proxy/swm_fortran_amrex_driver
$SWM_BUILD_DIR/swm_fortran/fortran_OpenACC/swm_fortran_acc
$SWM_BUILD_DIR/swm_amrex/swm_AMReX/swm_amrex $SWM_ROOT/swm_amrex/inputs
$SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fkernels/swm_amrex_fkernels $SWM_ROOT/swm_amrex/inputs
$SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenMP/swm_amrex_fsubroutine_omp $SWM_ROOT/swm_amrex/inputs
$SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenACC/swm_amrex_fsubroutine_acc $SWM_ROOT/swm_amrex/inputs

