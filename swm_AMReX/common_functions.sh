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

## Function to build Plotfile to HDF5 Conversion Executable.
## Sets the variable plotfile_2_hdf5_exe to the path of the executable.
#build_plotfile_to_hdf5_exe() {
#
#    local plotting_utils_dir="${SWM_AMREX_ROOT}"/plotting_utils
#
#    cd "${plotting_utils_dir}"
#
#    # Run make and capture the output
#    local make_output=$(mktemp) # Create a temporary file to store the output of make
#    #make | tee "$make_output"
#    #make > "$make_output" 2>&1
#    make > "$make_output" 
#
#    # Parse the output to find the executable name
#    local plotfile_2_hdf5_exe_base=$(grep "executable is" "$make_output" | awk '{print $3}')
#    rm "$make_output" # Remove the temporary file used to store the output of make
#
#    local plotfile_2_hdf5_exe="${plotting_utils_dir}/${plotfile_2_hdf5_exe_base}"
#
#    echo "$plotfile_2_hdf5_exe"
#}