#!/bin/bash

# Prerequisites:
#    AMREX_HOME = the directory where you pulled the AMReX repo.
#    SWM_AMREX_ROOT = the directory containing the AMReX version of the SWM mini-app.
#    HDF5_HOME = the directory where you installed HDF5. Only need if you are running with the "hdf5" solution verification method, or set create_plot="true" or create_movie="true".

##############################################################################
# User Input
##############################################################################

# Number of MPI ranks to run with. Only used if the executable name contains "MPI"
n_mpi_ranks=16
# Number of MPI ranks to run with. Only used if the executable name contains "OMP"
# Note that this script ignores the value of OMP_NUM_THREADS in the users environment and uses this hard coded value instead. Maybe not the best choice but just doing it fow now to know for sure how many threads I am asking for.
n_omp_threads=8

# Set the solution verification method: "none", "hdf5", or "plotfile"
#    You should also change the reference_file and file_to_compare variables in corresponding the Solution Verification section below.
solution_verification_method="none"       # Do not perform solution verification
#solution_verification_method="hdf5"       # If you use the hdf5 option you must have the HDF5_HOME environment variable set.
#solution_verification_method="plotfile"   # If you want to use the plotfile option you must first build the fcomprare executable in the AMREX_HOME/Tools/Plotfile directory.

# Convert all plotfiles to HDF5 format: "true" or "false"
convert_all_plotfiles_to_hdf5="false"

# Plot all results using python and matplotlib: "true" or "false"
#     If this option is set to "true" then the convert_all_plotfiles_to_hdf5 option will also be treated as "true". Since we need the hdf5 files to make the plots.
create_plots="false"

# Create an MP4 movie from the PNG files created by the plotting script: "true" or "false"
#     If this option is set to "true" then the create_plots option will also be treated as "true". Since we need the PNG files of the plots to make the movie.
create_movie="false"

##############################################################################
# Setup 
##############################################################################

set -e
set -u
#set -x

# Source common functions
source "$SWM_AMREX_ROOT"/common_functions.sh

##############################################################################
# Build
##############################################################################
print_banner "Build"

cd "$SWM_AMREX_ROOT"

# Run make and capture the output
make_output=$(mktemp) # Create a temporary file to store the output of make
make | tee "$make_output"
make_exit_status=${PIPESTATUS[0]} # Capture the exit status of make

if [ $make_exit_status -ne 0 ]; then
    echo "Make failed with exit status $make_exit_status"
    exit $make_exit_status
fi

# Parse the output to find the executable name
main_exe="${SWM_AMREX_ROOT}"/$(grep "executable is" "$make_output" | awk '{print $3}')
rm "$make_output" # Remove the temporary file used to store the output of make

##############################################################################
# Run
##############################################################################
print_banner "Run"

run_dir="$SWM_AMREX_ROOT"/run_dir
mkdir -p "$run_dir"
cd "$run_dir"

launch_prefix=""

if [[ "$main_exe" == *"OMP"* ]]; then
    launch_prefix+=" OMP_NUM_THREADS=$n_omp_threads "

    # This is an alternative way to set the number of threads but it is not nice to change the user's environment
    #export OMP_NUM_THREADS=$n_omp_threads 
fi

if [[ "$main_exe" == *"MPI"* ]]; then
    launch_prefix+=" mpiexec -n $n_mpi_ranks "
fi

swe_amrex_command="$launch_prefix $main_exe ${SWM_AMREX_ROOT}/inputs"
eval $swe_amrex_command

##############################################################################
# Solution Verification
##############################################################################
if [ "$solution_verification_method" != "none" ]; then

    print_banner "Solution Verification"
    echo " "

    if [ "$solution_verification_method" == "plotfile" ]; then
        
        # Look for files that match the pattern fcompare*.ex in the directory $AMREX_HOME/Tools/Plotfile/
        fcompare_files=($(find "$AMREX_HOME/Tools/Plotfile/" -name 'fcompare*.ex'))
        
        # Check if there is exactly one matching file
        if [ ${#fcompare_files[@]} -eq 1 ]; then
            fcompare_exe="${fcompare_files[0]}"
        elif [ ${#fcompare_files[@]} -gt 1 ]; then
            echo "Error: Multiple fcompare executables found:"
            for file in "${fcompare_files[@]}"; do
                echo "  $file"
            done
            echo "I am not sure which one to use. Aborting script."
            exit 1
        else
            echo "Error: No fcompare executable found."
            exit 1
        fi
  
        compare_exe="$fcompare_exe"
  
        # TODO: Make the plotfile name a so it automatically finds the plotfile for the last time step
        # TODO: Do a check to see if the corresponding reference plotfile exists step
        reference_file="$SWM_AMREX_ROOT"/plt00100_reference
        file_to_compare="$run_dir"/plt00100
  
    elif [ "$solution_verification_method" == "hdf5" ]; then

        compare_exe="h5diff" # exact comparison
        #compare_exe="h5diff --use-system-epsilon" # compare to machine epsilon
        #compare_exe="h5diff --relative=1.0e-2" # compare to 1% relative error
  
        #reference_file="$SWM_AMREX_ROOT"/plt00100_64_reference.h5
        reference_file="$SWM_AMREX_ROOT"/plt00100_4096_reference.h5
  
        ## Call the function and get the path to the executable
        #      This function sets the variable plotfile_2_hdf5_exe to the path of the executable.
        build_plotfile_to_hdf5_exe
        
        # Make sure that plotfile_2_hdf5_exe is an executable
        if [ ! -x "$plotfile_2_hdf5_exe" ]; then
            echo "Error: $plotfile_2_hdf5_exe is not an executable file."
            exit 1
        fi

        plt_file="$run_dir"/plt00100
        hdf5_file=${plt_file}.h5 # put the hdf5 file in the same directory as the plotfile and have the same base name as the plotfile, just append the .h5 extension
        "${plotfile_2_hdf5_exe}" infile="$plt_file" outfile="$hdf5_file"

        echo "Converted $plt_file to $hdf5_file"

        file_to_compare="$hdf5_file"
  
    else
        echo "Unknown solution verification method: $solution_verification_method"
        echo "The only supported options are: hdf5 and plotfile"
        exit 1
    fi
  
    echo "Comparing $reference_file to $file_to_compare"
  
    # Turn off exit on error for this one command since we want to check the return value
    set +e 
  
    $compare_exe $reference_file $file_to_compare
  
    if [ $? -eq 0 ]; then
        echo -e "\nSolution Verification: PASS"
    else
        echo -e "\nSolution Verification: FAIL"
        exit 1
    fi
  
    # Turn exit on error back on
    set -e 

fi

##############################################################################
# Convert all AMReX plotfiles to HDF5 
##############################################################################

if [ "$convert_all_plotfiles_to_hdf5" == "true" ] || [ "$create_plots" == "true" ] || [ "$create_movie" == "true" ]; then
    # For now this does all three steps: convert to hdf5, plot, and create movie
    # TODO break this up into separate functions to do each step
    $SWM_AMREX_ROOT/plotting_utils/convert_to_hdf5_and_plot.sh "$run_dir"
fi

##############################################################################
# Plotting 
##############################################################################

##############################################################################
# Create Movie 
##############################################################################
