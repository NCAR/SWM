#!/bin/bash -l
## to prepare a test case labeled by the name of this folder: define module path, load module, set library path, set environment variables
## the following template is translated from /glade/scratch/ssuresh/muram/pgi1910 bash script

module purge
module load cmake
module load gnu/8.3.0

export OMP_NUM_THREADS=18

# Set install_name to name of particular Kokkos install directory you would
# like to use to compile the source code
install_name=gcc_omp_serial_skx
# Set project_dir to the directory containing the source code you would like
# to compile
project_dir=/glade/work/$USER/SWM
# Set Kokkos_DIR environment variable to the full path to the particular Kokkos
# install directory you would like to use to compile the source code
export Kokkos_DIR=/glade/work/$USER/kokkos/installs/$install_name

# Check to see if Kokkos install exists
[ ! -d "$Kokkos_DIR" ] && echo "Build canceled - Requested Kokkos install doesn't exist" 
[ ! -d "$Kokkos_DIR" ] && exit 1

# Check to see if a build directory exists for this Kokkos install in the
# project directory
[ -d "$project_dir/$install_name" ] && rm -r $project_dir/$install_name

# Create a build directory, build, and generate an executable using chosen
# Kokkos install
mkdir $project_dir/$install_name
cd $project_dir/$install_name
cmake ..
make

./swm_kokkos > $project_dir/results.cpu.omp.$(date +%m%d%H%M%S).txt

