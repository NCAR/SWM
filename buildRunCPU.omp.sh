#!/bin/bash -l
## to prepare a test case labeled by the name of this folder: define module path, load module, set library path, set environment variables
## the following template is translated from /glade/scratch/ssuresh/muram/pgi1910 bash script

#PBS -N SWM_CPU_MP
#PBS -A NTDD0002
#PBS -l walltime=00:05:00
#PBS -q casper
### Merge output and error files
#PBS -o SWM_CPU_MP.out
#PBS -e SWM_CPU_MP.err
#PBS -l select=1:ncpus=1:ompthreads=18


#source ./loadEnv.sh

module purge
module load nvhpc/21.3
module list

export OMP_NUM_THREADS=18

nvc -mp shallow_unroll.acc.omp.c wtime.c -o SWM_cpu_mp
./SWM_cpu_mp > results.cpu.mp..$(date +%m%d%H%M%S).txt
