#!/bin/bash
#PBS -A NTDD0005
#PBS -N cf_gpu
#PBS -q develop@desched1
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=64:ngpus=1:gpu_type=a100
#PBS -o swm_cf_gpu.log

# Load modules to match compile-time environment
module --force purge
module load ncarenv-basic/25.10
module load nvhpc/26.1

# Set which compilers to use. These should be default with nvhpc loaded.
# export CC=nvc
# export CXX=nvc++
# export FC=nvfortran

# Set the SWM root and build directory 
export SWM_ROOT=~/SWM
export SWM_BUILD_DIR=$SWM_ROOT/../SWM_build

for num in 64 128 256 512 1024 2048 4096 8192 16384; do
    echo "========================================"
    echo "${num} x ${num}"
    echo "========================================"

    # Generate the build directory. 

    # To use compile time arguments: 
    # add these two lines to $SWM_ROOT/swm_c/c_OpenACC/CMakeLists.txt
    # target_compile_definitions(swm_c PRIVATE M=${M})
    # target_compile_definitions(swm_c PRIVATE N=${N})
    # and comment out lines 36 and 37 in $SWM_ROOT/swm_c/c_OpenACC/shallow_swap.acc.c
    # // #define M 256
    # // #define N 256
    # add these two lines to $SWM_ROOT/swm_fortran/fortran_OpenACC/CMakeLists.txt
    # target_compile_definitions(swm_fortran PRIVATE MNUM=${M})
    # target_compile_definitions(swm_fortran PRIVATE NNUM=${N})
    # and change out lines 5 and 6 in $SWM_ROOT/swm_fortran/common/params.F90 to
    # integer, parameter :: M = MNUM
    # integer, parameter :: N = NNUM
    cmake -DSWM_DEVICE=gpu -DSWM_C=ON -DSWM_FORTRAN=ON -DSWM_OPENACC=ON -S $SWM_ROOT -B $SWM_BUILD_DIR -DM=$num -DN=$num

    # or change hardcoded lines 36 and 37 in $SWM_ROOT/swm_c/c_OpenACC/shallow_swap.acc.c
    #                 and lines  5 and  6 in $SWM_ROOT/swm_fortran/common/params.F90
    # cmake -S $SWM_ROOT -B $SWM_BUILD_DIR

    # Build the code. 
    cmake --build $SWM_BUILD_DIR

    # Run the code. 
    echo
    echo "C"
    echo
    $SWM_BUILD_DIR/swm_c/c_OpenACC/swm_c_acc
    $SWM_BUILD_DIR/swm_c/c_OpenACC/swm_c_acc
    $SWM_BUILD_DIR/swm_c/c_OpenACC/swm_c_acc

    # Run the code. 
    echo
    echo "Fortran"
    echo
    $SWM_BUILD_DIR/swm_fortran/fortran_OpenACC/swm_fortran_acc
    $SWM_BUILD_DIR/swm_fortran/fortran_OpenACC/swm_fortran_acc
    $SWM_BUILD_DIR/swm_fortran/fortran_OpenACC/swm_fortran_acc

    echo
done
