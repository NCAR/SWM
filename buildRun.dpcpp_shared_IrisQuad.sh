#!/bin/bash -l

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
project_dir=/home/u76219/Leila/SWM
build_name=dpcpp_usm_IrisQuad
exec_name=swm_dpcpp_usm_shared

# Create build directories
[ ! -d "${project_dir}/build" ] && mkdir ${project_dir}/build
[ ! -d "${project_dir}/build/${build_name}" ] && mkdir ${project_dir}/build/${build_name}

# Make
[ ! -d "$exec_name" ] && make -B $exec_name
mv $exec_name ${project_dir}/build/${build_name}/.
cd ${project_dir}/build/${build_name}

# Create results directories if they don't already exist
[ ! -d "${project_dir}/results" ] && mkdir ${project_dir}/results
[ ! -d "${project_dir}/results/${build_name}" ] && mkdir ${project_dir}/results/${build_name}
[ ! -d "${project_dir}/results/${build_name}/csv" ] && mkdir ${project_dir}/results/${build_name}/csv

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
        ${project_dir}/build/${build_name}/$exec_name $M $N $build_name.$ID > $project_dir/results/$build_name/results.$M.$N.$build_name.$ID.txt;
    else
        ${project_dir}/build/${build_name}/$exec_name $M $N $build_name.$ID;
    fi
done

# Move csv file to results directory
mv $project_dir/build/$build_name/*.csv $project_dir/results/$build_name/csv

cd $orig_dir
