# Define common functions here

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

# Function to build Plotfile to HDF5 Conversion Executable.
#     Sets the variable plotfile_2_hdf5_exe to the path of the executable.
build_plotfile_to_hdf5_exe() {

    local plotting_utils_dir="${SWM_AMREX_ROOT}"/plotting_utils

    cd "${plotting_utils_dir}"

    # Run make and capture the output
    make_output=$(mktemp) # Create a temporary file to store the output of make
    make | tee "$make_output"
    make_exit_status=${PIPESTATUS[0]} # Capture the exit status of make
    
    if [ $make_exit_status -ne 0 ]; then
        echo "Make failed with exit status $make_exit_status"
        exit $make_exit_status
    fi

    # Parse the output to find the executable name
    local plotfile_2_hdf5_exe_base=$(grep "executable is" "$make_output" | awk '{print $3}')
    rm "$make_output" # Remove the temporary file used to store the output of make

    # This variable is the output of the function
    plotfile_2_hdf5_exe="${plotting_utils_dir}/${plotfile_2_hdf5_exe_base}"

    # Make sure that plotfile_2_hdf5_exe is an executable
    if [ ! -x "$plotfile_2_hdf5_exe" ]; then
        echo "Error: $plotfile_2_hdf5_exe is not an executable file."
        exit 1
    fi
}

