#!/bin/bash

# Prerequisite 
#    AMREX_HOME = the directory where you pulled the AMReX repo.
#    SWM_AMREX_ROOT = the directory containing the AMReX version of the SWM mini-app.

##############################################################################
# User Input
##############################################################################
# The exact executable name may vary depending on the compiler and build options you used.
# You may need to change this to match what is produced by your build.

main_exe="$SWM_AMREX_ROOT"/main2d.gnu.ex

#main_exe="$SWM_AMREX_ROOT"/main2d.gnu.TPROF.MPI.ex
#num_procs=4 # Number of MPI ranks to run with. Only used if the executable name contains "MPI"

fcompare_exe="$AMREX_HOME"/Tools/Plotfile/fcompare.gnu.ex

##############################################################################
# Setup 
##############################################################################

set -e
set -u

# Convenience function to print 80 character wide banners with centered text
print_banner() {
    local message="$1"
    local border_char="#"
    local width=80
    local padding=$(( (width - ${#message} - 2) / 2 ))
    local border=$(printf "%-${width}s" | tr ' ' "$border_char")

    echo -e "\n$border"
    printf "%*s %s %*s\n" $padding "" "$message" $padding ""
    echo -e "$border\n"
}

# Will put you back in this directory after script is finished
dir_backup=$(pwd)

##############################################################################
# Build
##############################################################################
print_banner "Build"
cd "$SWM_AMREX_ROOT"
make

##############################################################################
# Run
##############################################################################
print_banner "Run"

run_dir="$SWM_AMREX_ROOT"/run_dir
mkdir -p "$run_dir"
cd "$run_dir"

# TODO: Make the executable name a variable so that we can run the MPI, OpenMP, and Cuda versions with the same script
# Check if the executable name contains "MPI"
if [[ "$main_exe" == *"MPI"* ]]; then
    mpiexec -n $num_procs "$main_exe" "$SWM_AMREX_ROOT"/inputs
else
    "$main_exe" "$SWM_AMREX_ROOT"/inputs
fi

##############################################################################
# Solution Verification
##############################################################################
print_banner "Solution Verification"

# Turn off exit on error for this one command since we want to check the return value
set +e 

# TODO: Make the plotfile name a so it automatically finds the plotfile for the last time step
# TODO: Do a check to see if the corresponding reference plotfile exists step
$fcompare_exe plt04000 "$SWM_AMREX_ROOT"/plt04000_reference

if [ $? -eq 0 ]; then
    echo -e "\nSolution Verification: PASS"
else
    echo -e "\nSolution Verification: FAIL"
    exit 1
fi

# Turn exit on error back on
set -e 

##############################################################################
# Plotting 
##############################################################################
#print_banner "Plotting"
#python "$SWM_AMREX_ROOT"/plot_with_yt.py

##############################################################################
# Clean Up
##############################################################################
cd "$dir_backup"
