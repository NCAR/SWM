#!/bin/bash

# SWM_ROOT should be set and point to the directory where you cloned the SWM repository.
if [[ -z "$SWM_ROOT" || ! -d "$SWM_ROOT" ]]; then
    echo "Error: SWM_ROOT environment variable is not set or does not point to a valid directory. It should point to where you cloned the SWM repository." >&2
    exit 1
fi

# Set the where the build directory will be created... You can change this to your desired location 
export SWM_BUILD_DIR=$SWM_ROOT/../swm_build

# Set which compilers to use. 
export CC=gcc
export CXX=g++
export FC=gfortran

# Generate the build directory. 
cmake -S $SWM_ROOT -B $SWM_BUILD_DIR

# Build the code. 
cmake --build $SWM_BUILD_DIR

# Run the code. 
$SWM_BUILD_DIR/swm_c/c/swm_c
$SWM_BUILD_DIR/swm_fortran/fortran/swm_fortran