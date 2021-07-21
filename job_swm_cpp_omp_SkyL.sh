#!/bin/bash
source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1

echo
echo start: $(date "+%y%m%d.%H%M%S.%3N")
echo

export OMP_NUM_THREADS=36
./buildRun.cpp_omp_SkyL.sh -s -a
 
echo
echo stop:  $(date "+%y%m%d.%H%M%S.%3N")
echo

# qsub -l nodes=1:gold6128:ppn=2 -d . job_swm_cpp_omp_SkyL.sh