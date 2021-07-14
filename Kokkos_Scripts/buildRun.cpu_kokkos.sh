#!/bin/bash -l

# Get original working directory
orig_dir=$PWD

module purge
module load cmake
module load gnu/8.3.0

# Set install_name to name of particular Kokkos install directory you would
# like to use to compile the source code
install_name=gcc_serial_skx
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

./swm_kokkos
#[ ! -d "${project_dir}/results" ] && mkdir ${project_dir/results}
#[ ! -d "${project_dir}/results/cpu_kokkos" ] && mkdir ${project_dir}/results/cpu_kokkos
#./swm_kokkos > $project_dir/results/cpu_kokkos/results.cpu.kokkos.48.txt
#./swm_kokkos > $project_dir/results/cpu_kokkos/results.cpu.kokkos.$(date +%m%d%H%M%S).txt

# Return to original working directory
cd $orig_dir