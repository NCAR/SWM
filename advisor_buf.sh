#!/bin/bash

advixe-cl –collect=survey --enable-gpu-profiling --project-dir=./adv -- ./swm_dpcpp_buf

advixe-cl -–collect=tripcounts --stacks --flop --enable-gpu-profiling --project-dir=./adv -- ./swm_dpcpp_buf

advixe-cl --report=roofline --gpu  --project-dir=./adv
