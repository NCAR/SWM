#!/bin/bash -l

# Example PBS commands for Casper
#PBS -N SWM
#PBS -l nodes=1:ncpus=1:ngpus=1
#PBS -l walltime=00:05:00
#PBS -A NTDD0004
#PBS -q casper
#PBS -e SWM.err
#PBS -o SWM.out

# Example interactive job request for Casper
#execcasper -l gpu_type=v100 -l select=1:ncpus=18:ngpus=1 -A NTDD0004 -l walltime=06:00:00

module purge
module load nvhpc
module load cuda
module list
nvidia-smi

#make clean
#make    # set GPU=1 to enable GPU execution
#./SWM_acc > results_acc.txt

rm -f SWM_ori SWM_cpu SWM_gpu *.o  

# running the original code
#nvc shallow_swap.c wtime.c -o SWM_ori
#./SWM_ori > results_ref.64x64.IT4000.txt

# running the accelerated code on CPU
nvc shallow_swap.acc.c wtime.c -o SWM_cpu
./SWM_cpu > results.cpu.txt

# running the accelerated code on GPU
nvc -O2 -gpu=cc70 -Minfo -Mnofma shallow_swap.acc.c wtime.c -o SWM_gpu
./SWM_gpu > results.gpu.txt



