#!/bin/bash


###############################################################################
# Generated different builds for different compilers {gcc and nvidia} for different targets {cpu and gpu} in their own directory
# Build and run all the executables in each build directory
# This script assumes that the environment variable SWM_ROOT is set to the root directory of the SWM project
# It also assumes that the Spack package manager is installed and you have the gcc and nvhpc compilers installed via Spack
###############################################################################
spack unload --all
spack load gcc
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran -DDEVICE=cpu -S $SWM_ROOT -B $SWM_ROOT/build_gcc_cpu
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran -DDEVICE=gpu -S $SWM_ROOT -B $SWM_ROOT/build_gcc_gpu

spack unload --all
spack load nvhpc
cmake -DCMAKE_C_COMPILER=nvc -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_Fortran_COMPILER=nvfortran -DDEVICE=cpu -S $SWM_ROOT -B $SWM_ROOT/build_nvhpc_cpu
cmake -DCMAKE_C_COMPILER=nvc -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_Fortran_COMPILER=nvfortran -DDEVICE=gpu -S $SWM_ROOT -B $SWM_ROOT/build_nvhpc_gpu

# Loop over the four build directories and run make in each
for build_dir in $SWM_ROOT/build_gcc_cpu $SWM_ROOT/build_gcc_gpu $SWM_ROOT/build_nvhpc_cpu $SWM_ROOT/build_nvhpc_gpu; do

    spack unload --all
    if [[ "$build_dir" == *gcc* ]]; then
        spack load gcc
    fi
    if [[ "$build_dir" == *nvhpc* ]]; then
        spack load nvhpc
    fi

    echo "Building in $build_dir"
    cmake --build $build_dir
    $build_dir/swm_c/c/swm_c
    $build_dir/swm_c/c_OpenACC/swm_c_acc
    $build_dir/swm_c/c_OpenACC/swm_c_acc_tile
    $build_dir/swm_fortran/fortran/swm_fortran
    $build_dir/swm_fortran/fortran_OpenACC/swm_fortran_acc
    $build_dir/swm_fortran/fortran_AMReX_proxy/swm_fortran_amrex_driver

done


###############################################################################
# Create a one off build directory with taking the compiler from the environment variables
###############################################################################
#export CC=gcc
#export CXX=g++
#export FC=gfortran

#export CC=nvc
#export CXX=nvc++
#export FC=nvfortran
#
#cmake -S $SWM_ROOT -B $SWM_ROOT/build
