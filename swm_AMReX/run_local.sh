#!/bin/bash

# Prerequisite 
#    SWM_AMREX_ROOT must be set to the top level of the AMReX version of the mini-app.

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
"$SWM_AMREX_ROOT"/main2d.gnu.ex "$SWM_AMREX_ROOT"/inputs
#"$SWM_AMREX_ROOT"/main2d.gnu.MPI.ex "$SWM_AMREX_ROOT"/inputs

##############################################################################
# Solution Verification
##############################################################################
print_banner "Solution Verification"

# Turn off exit on error for this one command since we want to check the return value
set +e 

# TODO: Make the plotfile name a so it automatically finds the plotfile for the last time step
# TODO: Do a check to see if the corresponding reference plotfile exists step
$AMREX_HOME/Tools/Plotfile/fcompare.gnu.ex plt04000 "$SWM_AMREX_ROOT"/plt04000_reference

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
