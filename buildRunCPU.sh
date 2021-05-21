#!/bin/bash -l

module purge
module load gnu/9.1.0
module list

gcc -O2 -lm  shallow_unroll.c wtime.c -o SWM_cpu
./SWM_cpu > results.cpu.$(date +%m%d%H%M%S).txt
