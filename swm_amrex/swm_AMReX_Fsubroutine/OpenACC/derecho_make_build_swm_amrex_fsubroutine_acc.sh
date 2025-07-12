#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status
set -u  # Treat expanding empty variables as an error
set -x  # Print each command before executing it

###############################################################################
# User Input
###############################################################################

# Directory where you pulled the SWM repository
export SWM_ROOT=/glade/u/home/htorres/SWM

# Directory where you pulled the AMReX repositry
export AMREX_HOME=/glade/u/home/htorres/amrex

fresh_build=1  # Set to 1 to force a fresh build... A.K.A run make clean, 0 to use existing build

###############################################################################
# Module Setup
###############################################################################

module purge

# Modules we always use
module load nvhpc cuda

###############################################################################
# Build SWM AMReX with Fortran Subroutines that use OpenACC
###############################################################################
cd "$SWM_ROOT"/swm_amrex/swm_AMReX_Fsubroutine/OpenACC

if [[ $fresh_build -eq 1 ]]; then
  make clean
fi

make 2>&1 | tee build.log