#!/bin/bash
source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1
 
echo
echo start: $(date "+%y%m%d.%H%M%S.%3N")
echo

export SYCL_DEVICE_FILTER=gpu
./swm_dpcpp
 
echo
echo stop:  $(date "+%y%m%d.%H%M%S.%3N")
echo