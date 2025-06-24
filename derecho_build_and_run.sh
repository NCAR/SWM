#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status
set -u  # Treat expanding empty variables as an error
set -x  # Print each command before executing it

###############################################################################
# User Input
###############################################################################

# Directory where you pulled the AMReX repositry
export AMREX_HOME=/glade/u/home/htorres/amrex

# Directory where the amrex AMReX build will be created
# This is used as a base name and will appended to based on the options you pick
export AMREX_BUILD_DIR=/glade/u/home/htorres/amrex_build

# Directory where you pulled the SWM repository
export SWM_ROOT=/glade/u/home/htorres/SWM

# Directory where the SWM build will be created
# This is used as a base name and will appended to based on the options you pick
export SWM_BUILD_DIR=/glade/u/home/htorres/SWM_build

# Options... example to pick cpu vs gpu?
export AMREX_USE_MPI=NO   # Set to YES or NO

export AMREX_USE_CUDA=YES   # Set to YES or NO

###############################################################################
# Module Setup
###############################################################################

module purge

# Modules we always use
module load cmake

# Compiler
#module load gcc
#export SWM_BUILD_DIR="${SWM_BUILD_DIR}_GNU"
#export AMREX_BUILD_DIR="${AMREX_BUILD_DIR}_GNU"

module load nvhpc
export SWM_BUILD_DIR="${SWM_BUILD_DIR}_NVHPC"
export AMREX_BUILD_DIR="${AMREX_BUILD_DIR}_NVHPC"

if [[ "${AMREX_USE_MPI}" == "YES" ]]; then
    module load cray-mpich
    export SWM_BUILD_DIR="${SWM_BUILD_DIR}_MPI"
    export AMREX_BUILD_DIR="${AMREX_BUILD_DIR}_MPI"
fi

if [[ "${AMREX_USE_CUDA}" == "YES" ]]; then
    module load cuda
    export SWM_BUILD_DIR="${SWM_BUILD_DIR}_CUDA"
    export AMREX_BUILD_DIR="${AMREX_BUILD_DIR}_CUDA"
fi

module load ncarcompilers

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

module list

###############################################################################
# Build the version of AMReX that we are asking for
###############################################################################

# Initialize an array for CMake options
amrex_cmake_opts=()

# Always-used options
amrex_cmake_opts+=("-DAMReX_SPACEDIM=2")
amrex_cmake_opts+=("-DAMReX_PRECISION=DOUBLE")
amrex_cmake_opts+=("-DAMReX_FORTRAN=YES")
# These were on by default but I don't think we are using any of these features so I am turning them off
amrex_cmake_opts+=("-DAMReX_LINEAR_SOLVERS=NO")
amrex_cmake_opts+=("-DAMReX_LINEAR_SOLVERS_INCFLO=NO")
amrex_cmake_opts+=("-DAMReX_LINEAR_SOLVERS_EM=NO")
amrex_cmake_opts+=("-DAMReX_AMRLEVEL=NO")
amrex_cmake_opts+=("-DAMReX_PARTICLES=NO")
amrex_cmake_opts+=("-DAMReX_TINY_PROFILE=NO")

amrex_cmake_opts+=("-DAMReX_BUILD_SHARED_LIBS=YES")

#if [[ "${AMREX_USE_CUDA}" == "YES" ]]; then
#
#    #CMAKE_C_COMPILER=$(which nvc)
#    #CMAKE_CXX_COMPILER=$(which nvc++)
#    #CMAKE_Fortran_COMPILER=$(which nvfortran)
#    #amrex_cmake_opts+=("-DAMReX_DIFFERENT_COMPILER=ON")
#
#    amrex_cmake_opts+=("-DCMAKE_C_COMPILER=$(which nvc)")
#    amrex_cmake_opts+=("-DCMAKE_CXX_COMPILER=$(which nvc++)")
#    amrex_cmake_opts+=("-DCMAKE_Fortran_COMPILER=$(which nvfortran)")
#fi

amrex_cmake_opts+=("-DCMAKE_C_COMPILER=$CC")
amrex_cmake_opts+=("-DCMAKE_CXX_COMPILER=$CXX")
amrex_cmake_opts+=("-DCMAKE_Fortran_COMPILER=$FC")
amrex_cmake_opts+=("-DAMReX_DIFFERENT_COMPILER=ON")

# These take on differnent values depending on the options chosen
if [[ "${AMREX_USE_MPI}" == "YES" ]]; then
    amrex_cmake_opts+=("-DAMReX_MPI=YES")
else
    amrex_cmake_opts+=("-DAMReX_MPI=NO")
fi

amrex_cmake_opts+=("-DAMReX_OMP=NO")

if [[ "${AMREX_USE_CUDA}" == "YES" ]]; then
    amrex_cmake_opts+=("-DAMReX_GPU_BACKEND=CUDA")
    amrex_cmake_opts+=("-DAMReX_GPU_RDC=YES")
    amrex_cmake_opts+=("-DAMReX_CUDA_ARCH=8.0") # Set the CUDA architecture version, adjust as needed
else
    amrex_cmake_opts+=("-DAMReX_GPU_BACKEND=NONE")
fi

export AMREX_INSTALL_DIR=$AMREX_BUILD_DIR/install

mkdir -p $AMREX_BUILD_DIR
cd $AMREX_BUILD_DIR

cmake "${amrex_cmake_opts[@]}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$AMREX_INSTALL_DIR" \
    -S "$AMREX_HOME" \
    -B "$AMREX_BUILD_DIR"
#    --trace \

make -j 32 install 
#make -j 32 install 2>&1 | less -R
#make test_install  # optional step to test if the installation is working

#exit 0 # Exit early for testing purposes

###############################################################################
# Build SWM Using the version of AMReX that we just built
###############################################################################
#cmake -DAMReX_ROOT=$AMREX_INSTALL_DIR/lib/cmake/AMReX \
#      -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_Fortran_COMPILER=$FC \
#      -S $SWM_ROOT -B $SWM_BUILD_DIR

#cmake -DAMReX_ROOT=$AMREX_INSTALL_DIR/lib/cmake/AMReX \
#      -DAMReX_GPU_BACKEND=CUDA \
#      -DCMAKE_C_COMPILER=nvc -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_Fortran_COMPILER=nvfortran \
#      -DCMAKE_CUDA_FLAGS="-O3 -DNDEBUG  --expt-relaxed-constexpr --expt-extended-lambda -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored -Xcudafe --diag_suppress=implicit_return_from_non_void_function -maxrregcount=255 -Xcudafe --display_error_number --Wext-lambda-captures-this --use_fast_math --generate-line-info" \
#      -S $SWM_ROOT -B $SWM_BUILD_DIR

cmake -DAMReX_ROOT=$AMREX_INSTALL_DIR/lib/cmake/AMReX \
      -DAMReX_GPU_BACKEND=CUDA \
      -DCMAKE_C_COMPILER=nvc -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_Fortran_COMPILER=nvfortran \
      -S $SWM_ROOT -B $SWM_BUILD_DIR

#cmake -DCMAKE_C_COMPILER=nvc -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_Fortran_COMPILER=nvfortran \
#      -S $SWM_ROOT -B $SWM_BUILD_DIR


cd $SWM_BUILD_DIR
make 
#make 2>&1 | log.txt

###############################################################################
# Run all the versions of SWM mini-app
###############################################################################
#$SWM_BUILD_DIR/swm_c/c/swm_c
#$SWM_BUILD_DIR/swm_c/c_OpenACC/swm_c_acc
#$SWM_BUILD_DIR/swm_c/c_OpenACC/swm_c_acc_tile
#$SWM_BUILD_DIR/swm_fortran/fortran/swm_fortran
#$SWM_BUILD_DIR/swm_fortran/fortran_AMReX_proxy/swm_fortran_amrex_driver
#$SWM_BUILD_DIR/swm_fortran/fortran_OpenACC/swm_fortran_acc

$SWM_BUILD_DIR/swm_amrex/swm_AMReX/swm_amrex $SWM_ROOT/swm_amrex/inputs
$SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fkernels/swm_amrex_fkernels $SWM_ROOT/swm_amrex/inputs
$SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenMP/swm_amrex_fsubroutine_omp $SWM_ROOT/swm_amrex/inputs
$SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenACC/swm_amrex_fsubroutine_acc $SWM_ROOT/swm_amrex/inputs

if [[ "${AMREX_USE_MPI}" == "YES" ]]; then
    mpirun -np 2 $SWM_BUILD_DIR/swm_amrex/swm_AMReX/swm_amrex $SWM_ROOT/swm_amrex/inputs
    mpirun -np 2 $SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fkernels/swm_amrex_fkernels $SWM_ROOT/swm_amrex/inputs
    mpirun -np 2 $SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenMP/swm_amrex_fsubroutine_omp $SWM_ROOT/swm_amrex/inputs
    mpirun -np 2 $SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenACC/swm_amrex_fsubroutine_acc $SWM_ROOT/swm_amrex/inputs
else
    $SWM_BUILD_DIR/swm_amrex/swm_AMReX/swm_amrex $SWM_ROOT/swm_amrex/inputs
    $SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fkernels/swm_amrex_fkernels $SWM_ROOT/swm_amrex/inputs
    $SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenMP/swm_amrex_fsubroutine_omp $SWM_ROOT/swm_amrex/inputs
    $SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenACC/swm_amrex_fsubroutine_acc $SWM_ROOT/swm_amrex/inputs
fi

#echo "C Compiler:   $CC"
#echo "C++ Compiler: $CXX"
#echo "Fortran Compiler: $FC"