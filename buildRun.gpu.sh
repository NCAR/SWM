#!/bin/bash -l

#PBS -N SWM_GPU
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

orig_dir=$PWD

module purge
module load nvhpc/21.3
module load cuda
module list
nvidia-smi

project_dir=/glade/work/$USER/SWM
build_name=gpu

# Check to see if a build directory exists for this Kokkos install in the
# project directory
[ -d "$project_dir/$build_name" ] && rm -r $project_dir/$build_name

# Create a build directory, build, and generate an executable using chosen
# Kokkos install
mkdir $project_dir/$build_name
cd $project_dir/$build_name
nvc -O3 -acc -ta=tesla:cc70 -Minfo -Mnofma $project_dir/shallow_unroll.acc.omp.c $project_dir/wtime.c -o SWM_$build_name

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
        ./SWM_$build_name $M $N $build_name.$ID > $project_dir/results/$build_name/results.$M.$N.$build_name.$ID.txt;
    else
        ./SWM_$build_name $M $N $build_name.$ID;
    fi
done

# Move csv file to kokkos gpu results directory in the working directory
mv $project_dir/$build_name/*.csv $project_dir/results/$build_name

cd $orig_dir
