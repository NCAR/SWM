#!/bin/bash

# Prerequisites:
#    SWM_AMREX_ROOT = the directory containing the AMReX version of the SWM mini-app.
#    You must provide an executable directory.

##############################################################################
# User Input
##############################################################################

#n_processors=(1 2)
#n_processors=(1 2 4 8 16)
#n_processors=(1 2 4 8 16 32 64 128)
n_processors=(16 32 64 128)
#n_processors=(1 2 4 8 16 32 64 128 256)
#n_processors=(32 64 128 256)
#n_processors=($(seq 1 16))

n_samples_per_case=5

#main_exe="${SWM_AMREX_ROOT}"/main2d.gnu.TPROF.MPI.ex
main_exe="${SWM_AMREX_ROOT}"/main2d.gnu.x86-milan.TPROF.MPI.ex
#main_exe="${SWM_AMREX_ROOT}"/main2d.intel-llvm.x86-milan.TPROF.MPI.ex

inputs_file="${SWM_AMREX_ROOT}"/inputs

output_dir="${SWM_AMREX_ROOT}"/strong_scaling_runs

timer_names=("main" "UpdateIntermediateVariables" "UpdateNewVariables" "UpdateOldVariables" "FabArray::FillBoundary")

##############################################################################
# Setup
##############################################################################

set -e
set -u
#set -x

mkdir -p "${output_dir}"

# Save a copy of the inputs file we are using can be usefull for reference later when you copy this whole folder to look at the scaling output.
cd "${output_dir}"
cp $inputs_file inputs.txt

# Print headers - Note this overwrites the file if it already exists
echo "n_proc, run_idx, runtime_min [s], runtime_avg [s], runtime_max [s]" > "${output_dir}/total_runtime.txt"
for timer_type in "inclusive" "exclusive"; do
  for timer_name in "${timer_names[@]}"; do
      echo "n_proc, run_idx, runtime_min [s], runtime_avg [s], runtime_max [s], max_percent" > "${output_dir}/${timer_type}_runtime_${timer_name}.txt"
  done
done

##############################################################################
# Run all the Cases and Extract the Timer Data into Summary CSV Files
##############################################################################

# Looping through the array
for n_proc in "${n_processors[@]}"; do
    for ((run_idx=1; run_idx<=n_samples_per_case; run_idx++)); do

        #######################################################################
        # Run the case
        #######################################################################
        echo " Running case with n_proc $n_proc sample $run_idx"

        run_dir="${output_dir}/n_proc_$(printf "%04d" ${n_proc})/run_number_$(printf "%04d" ${run_idx})"
        mkdir -p "${run_dir}"
        cd "${run_dir}"

        #if ! mpiexec -n ${n_proc} ${main_exe} ${inputs_file} > output.txt 2>&1; then
        #if ! mpiexec -n ${n_proc} --cpu-bind none --mem-bind none ${main_exe} ${inputs_file} > output.txt 2>&1; then
	if ! map --profile --report=txt,summary mpiexec -n ${n_proc} -ppn ${n_proc} --cpu-bind none -env OMP_NUM_THREADS=1 ${main_exe} ${inputs_file} > output.txt 2>&1; then
            echo "Error: mpiexec failed for n_proc=${n_proc}, run_idx=${run_idx}"
            echo "Output Saved Saved to file: ${run_dir}"
            cat output.txt
            exit 1
        fi

        #######################################################################
        # Extract Total Runtime from Output File
        #######################################################################
        runtime_line=$(grep "TinyProfiler total time across processes" output.txt)

        runtime_min=$(echo "$runtime_line" | awk -F': ' '{print $2}' | awk -F' ... ' '{print $1}')
        runtime_avg=$(echo "$runtime_line" | awk -F': ' '{print $2}' | awk -F' ... ' '{print $2}')
        runtime_max=$(echo "$runtime_line" | awk -F': ' '{print $2}' | awk -F' ... ' '{print $3}')

        echo "${n_proc}, ${run_idx}, ${runtime_min}, ${runtime_avg}, ${runtime_max}" >> "${output_dir}/total_runtime.txt"

        #######################################################################
        # Extract Timers from Output File
        #######################################################################
        # Extract the inclusive and exclusive runtimes  for the selected timers from the output file
        for timer_type in "exclusive" "inclusive" ; do
            for timer_name in "${timer_names[@]}"; do

                if [ "$timer_type" == "exclusive" ]; then
                    grep -A 10 "^Name\s\+NCalls\s\+Excl\. Min\s\+Excl\. Avg\s\+Excl\. Max\s\+Max %$" output.txt > table.txt
                else
                   grep -A 10 "^Name\s\+NCalls\s\+Incl\. Min\s\+Incl\. Avg\s\+Incl\. Max\s\+Max %$" output.txt > table.txt
                fi

                runtime_min=$(awk '/^'"${timer_name}"'/ {print $3}' table.txt)
                runtime_avg=$(awk '/^'"${timer_name}"'/ {print $4}' table.txt)
                runtime_max=$(awk '/^'"${timer_name}"'/ {print $5}' table.txt)
                percent_max=$(awk '/^'"${timer_name}"'/ {print $6}' table.txt | sed 's/%//') # Remove the percent sign

		rm table.txt

                echo "${n_proc}, ${run_idx}, ${runtime_min}, ${runtime_avg}, ${runtime_max}, ${percent_max}" >> "${output_dir}/${timer_type}_runtime_${timer_name}.txt"

            done
        done

    done
done
