#!/bin/bash
#PBS -A NTDD0005
#PBS -N python_serial
#PBS -q develop@desched1
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=128
#PBS -o python_serial.log

# Load modules to match compile-time environment
module --force purge
module load ncarenv-basic/25.10
module load conda
conda activate npl

# Set python directory
export SWM_PYTHON=~/SWM/swm_python
cd $SWM_PYTHON

# Build and run the code.

for num in 64 128 256 512 1024; do
    echo "========================================"
    echo "${num} x ${num}"
    echo "========================================"

    python swm_numpy.py --M=$num --N=$num
    python swm_numpy.py --M=$num --N=$num
    python swm_numpy.py --M=$num --N=$num

    echo
done