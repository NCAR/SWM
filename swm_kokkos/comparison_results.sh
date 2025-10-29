#!/bin/bash

# This script is used to compare the results of SWM_C and SWM_KOKKOS implementations.
# Works for Derecho only.

# 1. Enable the conda environment on Derecho
module purge
module load conda
conda activate npl

# 2. Define the paths to the output files from both implementations
path1="/glade/derecho/scratch/sunjian/SWM_KOKKOS_BUILD/swm_c/c/"
path2="/glade/derecho/scratch/sunjian/SWM_KOKKOS_BUILD/swm_kokkos/"

# 3. Define output file names to compare
file_lists=("p.txt" "u.txt" "v.txt")

# 4. Loop through each file and compare
for file in "${file_lists[@]}"; do
    echo "...Comparing $file between SWM_C and SWM_KOKKOS..."
    python comparison_results.py "${path1}${file}" "${path2}${file}"
    echo ""
done