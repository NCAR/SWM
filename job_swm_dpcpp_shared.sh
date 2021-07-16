#!/bin/bash
source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1
 
echo
echo start: $(date "+%y%m%d.%H%M%S.%3N")
echo

#export SYCL_DEVICE_FILTER=gpu
./buildRun.dpcpp_usm_shared.sh -s -a
 
echo
echo stop:  $(date "+%y%m%d.%H%M%S.%3N")
echo