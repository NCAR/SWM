#!/bin/bash

# Get original working directory
orig_dir=$PWD

# Load required modules
module load gnu

# Change to directory of desire kokkos profiling tool
profiler_dir=/glade/work/$USER/kokkos-tools/profiling/simple-kernel-timer
cd $profiler_dir

# Build kokkos profiling tool
[ -f *.so ] && make clean
make

# Set profiler environment variable so Kokkos can find the dynamic libraries on initialize
export KOKKOS_PROFILE_LIBRARY=${profiler_dir}/kp_kernel_timer.so

# Add profiler directory to path so kp_reader command can be used from anywhere
export PATH=${PATH}:${profiler_dir}

# Return to original working directory
cd $orig_dir