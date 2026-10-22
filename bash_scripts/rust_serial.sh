#!/bin/bash
#PBS -A NTDD0005
#PBS -N rust_serial
#PBS -q develop@desched1
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=128
#PBS -o rust_serial.log

# Load modules to match compile-time environment
module --force purge

# Set the rust directory 
export SWM_RUST=~/SWM/swm_rust/zip_ndarray
cd $SWM_RUST

# Build and run the code.

for num in 64 128 256 512 1024; do
    echo "========================================"
    echo "${num} x ${num}"
    echo "========================================"

    M=$num N=$num RUSTFLAGS="-C target-cpu=native" cargo build --release

    ./target/release/swm_rust
    ./target/release/swm_rust
    ./target/release/swm_rust

    echo
done