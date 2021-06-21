#!/bin/bash -l

diff swm_init.csv swm_init_orig.csv >error_init

diff swm_h100.csv swm_h100_orig.csv >error_h100
