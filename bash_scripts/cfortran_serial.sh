#!/bin/bash
#PBS -A NTDD0005
#PBS -N cfortran_serial
#PBS -q develop@desched1
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=128
#PBS -o cfortran_serial.log

# Load modules to match compile-time environment
module --force purge
module load ncarenv/25.10 intel/2025.2.1 cray-mpich/8.1.32

# Set the SWM root and build directory 
export SWM_ROOT=~/SWM
export SWM_BUILD_DIR=$SWM_ROOT/../SWM_build

# Set which compilers to use.
export CC=gcc
export CXX=g++
export FC=gfortran

# Build and run the code.

for num in 64 128 256 512 1024; do
    echo "========================================"
    echo "${num} x ${num}"
    echo "========================================"

    # Generate the build directory.

    # To use compile time arguments: 
    # add these two lines to $SWM_ROOT/swm_c/c/CMakeLists.txt
    # target_compile_definitions(swm_c PRIVATE M=${M})
    # target_compile_definitions(swm_c PRIVATE N=${N})
    # and comment out lines 34 and 35 in $SWM_ROOT/swm_c/c/shallow_swap.c
    # // #define M 256
    # // #define N 256
    # add these two lines to $SWM_ROOT/swm_fortran/fortran/CMakeLists.txt
    # target_compile_definitions(swm_fortran PRIVATE MNUM=${M})
    # target_compile_definitions(swm_fortran PRIVATE NNUM=${N})
    # and change out lines 5 and 6 in $SWM_ROOT/swm_fortran/common/params.F90 to
    # integer, parameter :: M = MNUM
    # integer, parameter :: N = NNUM
    cmake -S $SWM_ROOT -B $SWM_BUILD_DIR -DM=$num -DN=$num

    # or change hardcoded lines 34 and 35 in $SWM_ROOT/swm_c/c/shallow_swap.c
    #                 and lines  5 and  6 in $SWM_ROOT/swm_fortran/common/params.F90
    # cmake -S $SWM_ROOT -B $SWM_BUILD_DIR

    # Build the code. 
    cmake --build $SWM_BUILD_DIR

    # Run the C code. 
    $SWM_BUILD_DIR/swm_c/c/swm_c
    $SWM_BUILD_DIR/swm_c/c/swm_c
    $SWM_BUILD_DIR/swm_c/c/swm_c
    $SWM_BUILD_DIR/swm_c/c/swm_c
    $SWM_BUILD_DIR/swm_c/c/swm_c

    # Run the Fortran code. 
    $SWM_BUILD_DIR/swm_fortran/fortran/swm_fortran
    $SWM_BUILD_DIR/swm_fortran/fortran/swm_fortran
    $SWM_BUILD_DIR/swm_fortran/fortran/swm_fortran
    $SWM_BUILD_DIR/swm_fortran/fortran/swm_fortran
    $SWM_BUILD_DIR/swm_fortran/fortran/swm_fortran

    echo
done