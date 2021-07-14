#!/bin/bash -l

#PBS -N SWM_CPU_MP
#PBS -A NTDD0002
#PBS -l walltime=00:05:00
#PBS -q casper
### Merge output and error files
#PBS -o SWM_GPU.out
#PBS -e SWM_GPU.err
#PBS -l select=1:ncpus=1:ngpus=1
#PBS -l gpu_type=v100

module purge
module load nvhpc/21.3
module load cuda
module list
nvidia-smi

rm SWM_gpu

nvc -O2 -acc -ta=tesla:cc70 -Minfo -Mnofma shallow_unroll.acc.omp.c wtime.c -o SWM_gpu
./SWM_gpu > results.gpu.$(date +%m%d%H%M%S).txt

