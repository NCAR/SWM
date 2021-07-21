#!/bin/bash
source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1

echo
echo start: $(date "+%y%m%d.%H%M%S.%3N")
echo

export OverrideDefaultFP64Settings=1
export IGC_EnableDPEmulation=1

./buildRun.dpcpp_shared_IrisQuad.sh -s -a
 
echo
echo stop:  $(date "+%y%m%d.%H%M%S.%3N")
echo

# qsub -l nodes=1:iris_xe_max:quad_gpu:ppn=2 -d . job_swm_dpcpp_shared_IrisQuad.sh
