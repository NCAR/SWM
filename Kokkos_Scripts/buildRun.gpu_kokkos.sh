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

ID="date"

for arg in "$@"; do
    case $arg in
        -s|--save) 
            SAVE="true";
            shift;;
        -i|--id)
            ID=$2;
            shift;
            shift;;
        -m|--hsize)
            M=$2;
            shift;
            shift;;
        -n|--vsize)
            N=$2;
            shift;
            shift;;
        -a|--all)
            ALL="true";
            shift;;
        -h|--help)
            echo "-s|--save        : Save console output";
            echo "-i|--id <name>   : Set identifier to problem size with <name> in saved output filename";
            echo "                   If -f option not included, identifier defaults to problem size with time stamp";
            echo "-m|--hsize <int> : Set the horizontal dimension of the problem size to <int>";
            echo "                   If -n option not included, set both dimensions to <int>";
            echo "                   If neither -m or -n options included, run on default problem size 64x128";
            echo "-n|--vsize <int> : Set the vertical dimension of the problem size to <int>";
            echo "                   If -m option not included, set both dimensions to <int>";
            echo "                   If neither -m or -n options included, run on default problem size 64x128";
            echo "-a|--all         : Run all of the problem sizes in the following predefined set -"
            echo "                   (48x48 64x64 96x96 128x128 192x192 256x256 384x384 512x512 768x768 1024x1024 1536x1536 2048x2048 3072x3072";
            echo "-h|--help        : Display this help output"
            shift;
            exit 0;;
    esac
done

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
build_name=gpu_kokkos

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

# Create results directories if they don't already exist
[ ! -d "${project_dir}/results" ] && mkdir ${project_dir}/results
[ ! -d "${project_dir}/results/${build_name}" ] && mkdir ${project_dir}/results/${build_name}

# Set problem size based on command line arguments or lack there of
if [ "$ALL" == "true" ]; then
    Ms=(48 64 96 128 192 256 384 512 768 1024 1536 2048 3072);
    Ns=(48 64 96 128 192 256 384 512 768 1024 1536 2048 3072);
elif [ ! -z $M ]; then
    Ms=($M);
    if [ -z $N ]; then
        N=$M;
    fi
    Ns=($N);
elif [ ! -z $N ]; then
    M=$N;
    Ms=($M);
    Ns=($N);
else
    Ms=(64);
    Ns=(128);
fi

# Execute and save results if selected
NUM_PROBS=${#Ms[@]}
for (( i=0; i<$NUM_PROBS; i++ )) do
    M=${Ms[$i]};
    N=${Ns[$i]};

    if [ "$ID" == "date" ]; then
        ID=$(date +%m%d%H%M%S);
    fi

    if [ "$SAVE" == "true" ]; then
        ./swm_kokkos $M $N $build_name.$ID > $project_dir/results/$build_name/results.$M.$N.$build_name.$ID.txt;
    else
        ./swm_kokkos $M $N $build_name.$ID;
    fi
done

# Move csv file to kokkos gpu results directory in the working directory
mv *.csv $project_dir/results/$build_name

cd $orig_dir
