#!/bin/bash
source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1

echo
echo start: $(date "+%y%m%d.%H%M%S.%3N")
echo

./buildRun.dpcpp_shared_P630.sh -s -a
 
echo
echo stop:  $(date "+%y%m%d.%H%M%S.%3N")
echo

# qsub -l nodes=1:gpu:ppn=2 -d . job_swm_dpcpp_shared_P630.sh