#!/bin/bash
#PBS -A NTDD0005
#PBS -N rust_gpu
#PBS -q develop@desched1
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=64:ngpus=1:gpu_type=a100
#PBS -o rust_gpu.log

# Load modules to match compile-time environment
module --force purge
module load ncarenv-basic/25.10
module load cuda/12.9.0

# Set the rust directory 
export SWM_RUST=~/SWM/swm_rust/gpu_cudarc
cd $SWM_RUST

# Build and run the code.

for num in 64 128 256 512 1024 2048 4096 8192 16834; do
    echo "========================================"
    echo "${num} x ${num}"
    echo "========================================"

    M=$num N=$num RUSTFLAGS="-C target-cpu=native" cargo build --release

    ./target/release/swm_rust
    ./target/release/swm_rust
    ./target/release/swm_rust

    echo
done
