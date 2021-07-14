#!/bin/bash
source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1
 
echo
echo start: $(date "+%y%m%d.%H%M%S.%3N")
echo

export SYCL_DEVICE_FILTER=gpu
./swm_dpcpp_usm >$PWD/results/results.dpcpp_usm.$(date +%m%d%H%M%S).txt
 
echo
echo stop:  $(date "+%y%m%d.%H%M%S.%3N")
echo