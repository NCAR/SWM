#!/bin/bash -l

#PBS -N SWM_GPU_CPP
#PBS -A NTDD0002
#PBS -l walltime=00:05:00
#PBS -q casper
### Merge output and error files
#PBS -o SWM_GPU_CPP.out
#PBS -e SWM_GPU_CPP.err
#PBS -l select=1:ncpus=1:ngpus=1
#PBS -l gpu_type=v100

module purge
module load nvhpc/21.3
module load cuda
module list
nvidia-smi

[ -f SWM_gpu_cpp ] && rm SWM_gpu_cpp

nvcc -O3 -arch=compute_70 -std=c++11 shallow_unroll.acc.omp.cpp wtime.cpp -o SWM_gpu_cpp
./SWM_gpu_cpp > results.gpu.cpp.$(date +%m%d%H%M%S).txt
