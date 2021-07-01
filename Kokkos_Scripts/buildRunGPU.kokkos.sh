#!/bin/bash -l

#PBS -N SWM_CPU_MP
#PBS -A NTDD0002
#PBS -l walltime=00:05:00
#PBS -q casper
### Merge output and error files
#PBS -o SWM_GPU.out
#PBS -e SWM_GPU.err
#PBS -l select=1:ncpus=1:ngpus=1
#PBS -l gpu_type=v100

# Get original working directory
orig_dir=$PWD

#module purge
module load cmake
module load cuda/10.1
module load gnu/8.3.0

# Set install_name to name of particular Kokkos install directory you would
# like to use to compile the source code
install_name=nvcc_cuda_omp_serial_volta70_skx
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

# Set environment variables to squash warning
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

./swm_kokkos > $project_dir/results.gpu.$(date +%m%d%H%M%S).txt

# Return to original working directory
cd $orig_dir