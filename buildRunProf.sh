#!/bin/bash -l
## to prepare a test case labeled by the name of this folder: define module path, load module, set library path, set environment variables
## the following template is translated from /glade/scratch/ssuresh/muram/pgi1910 bash script

#SBATCH -J MUR_B
#SBATCH -n 1
#SBATCH --ntasks-per-node=18
#SBATCH -t 00:10:00
#SBATCH -A NTDD0002
#SBATCH -p dav
##SBATCH --reservation=NTDD0002
##SBATCH --reservation=$caspResv
#SBATCH --reservation=TDD_4xV100
#SBATCH -e MURaMB.err
#SBATCH -o MURaMB.out

#source ./loadEnv.sh

module purge
module load pgi/20.4
module load cuda
module list
nvidia-smi

pgcc -O0 -g -acc -ta=tesla:cc70 -Minfo -Mnofma shallow_swap.acc.Tile.c wtime.c -o SWM_nvprof
nvprof --print-gpu-trace --print-openacc-trace ./SWM_nvprof > results.gpuprof.txt
