#!/bin/bash

# Prerequisite 
#    SWM_AMREX_ROOT must be set to the top level of the AMReX version of the mini-app.

set -e
set -u

# Will put you back in this directory after script is finished
dir_backup=$(pwd)

##############################################################################
# Build
##############################################################################
cd "$SWM_AMREX_ROOT"
make

##############################################################################
# Run
##############################################################################
run_dir="$SWM_AMREX_ROOT"/run_dir
mkdir -p "$run_dir"
cd "$run_dir"

"$SWM_AMREX_ROOT"/main2d.gnu.ex "$SWM_AMREX_ROOT"/inputs

##############################################################################
# Postprocess
##############################################################################
python "$SWM_AMREX_ROOT"/plot_with_yt.py


##############################################################################
# Clean Up
##############################################################################
cd "$dir_backup"