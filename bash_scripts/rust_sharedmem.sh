#!/bin/bash
#PBS -A NTDD0005
#PBS -N rust_sharedmem
#PBS -q develop@desched1
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=128:ompthreads=128:mem=235GB
#PBS -o rust_sharedmem.log

# Load modules to match compile-time environment
module --force purge

# Set the rust directory 
export SWM_RUST=~/SWM/swm_rust/zip_ndarray_rayon
cd $SWM_RUST

# Build and run the code.

# choose grid size
M=4096 N=4096 RUSTFLAGS="-C target-cpu=native" cargo build --release

for num in 2 4 8 16 32 64 128; do
    echo "========================================"
    echo "${num} cores"
    echo "========================================"

    export RAYON_NUM_THREADS=$num

    ./target/release/swm_rust
    ./target/release/swm_rust
    ./target/release/swm_rust

    echo
done
