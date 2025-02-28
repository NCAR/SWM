#!/bin/bash

# This script is used to convert all AMReX plotfiles in a specified directory to HDF5 format
# and then plot the hdf5 files using matplotlib. The script will also create an MP4 movie from the
# PNG files created by the plotting script.

set -u
set -e
#set -x

###############################################################################
# Command Line Input
###############################################################################
# Input directory argument must be provided
if [ -z "$1" ]; then
    echo "Usage: $0 <directory with amrex plotfiles of the form plt_*>"
    exit 1
fi

# Set the directory to the provided argument
input_directory="$1"

# Remove trailing slash from input_directory if it exists... so the pattern in the for loop below works
input_directory="${input_directory%/}"

###############################################################################
# Setup
###############################################################################

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

###############################################################################
# Build Plotfile to HDF5 Conversion Executable
###############################################################################

print_banner "Build Plotfile to HDF5 Conversion Executable"

plotting_utils_dir="${SWM_AMREX_ROOT}"/plotting_utils

cd "${plotting_utils_dir}"

# Run make and capture the output
make_output=$(mktemp) # Create a temporary file to store the output of make
make | tee "$make_output"

# Parse the output to find the executable name
plotfile_2_hdf5_exe_base=$(grep "executable is" "$make_output" | awk '{print $3}')
rm "$make_output" # Remove the temporary file used to store the output of make

plotfile_2_hdf5_exe="${plotting_utils_dir}/${plotfile_2_hdf5_exe_base}"

###############################################################################
# Convert AMReX Plotfiles to HDF5 files
###############################################################################

# The directory where the HDF5 files will be saved... just use the same directory that contains the plotfiles
hdf5_output_dir="$input_directory"

# Initialize an array to store the HDF5 filenames
hdf5_files=()

# Loop over all files that start with plt_ followed by a series of numbers in the specified directory
for plt_file in "$input_directory"/plt[0-9]*; do

    # Skip files that contain .old or .h5
    if [[ "$plt_file" == *".old"* || "$plt_file" == *".h5"* ]]; then
        continue
    fi

    print_banner "Converting $plt_file to HDF5"

    # Plotfiles are actually directories. Check that plt_file is a directory
    if [[ -d "$plt_file" ]]; then
        # Set the outfile name to the same as plt_file but with .h5 extension
        hdf5_file="${hdf5_output_dir}/$(basename "${plt_file}").h5"
        
        "$plotfile_2_hdf5_exe" infile="$plt_file" outfile="$hdf5_file"
        
        # Add the HDF5 filename to the array
        hdf5_files+=("$hdf5_file")
    fi

done

##############################################################################
# Loop over each HDF5 files and plot. Creates a series of PNG files.
##############################################################################
print_banner "Plotting HDF5 files"

# The directory where the png files will be saved... just use the same directory that contains the plotfiles
png_output_dir="$input_directory"

for hdf5_file in "${hdf5_files[@]}"; do
    python "$SWM_AMREX_ROOT"/plotting_utils/plot_hdf5.py "$hdf5_file" --output_dir "$png_output_dir"
done

##############################################################################
# Create MP4 Movie from PNG Files
##############################################################################
print_banner "Creating MP4 Movies"

# Ensure ffmpeg is available
if ! command -v ffmpeg &> /dev/null
then
    echo "ffmpeg could not be found, please install it or add it to your path to create the MP4 movies."
    exit 1
fi

# The directory where the mp4 files will be saved... just use the same directory that contains the plotfiles
mp4_output_dir="$input_directory"

# Create the MP4 movie from PNG files
for prefix in u v p; do
    ffmpeg -y -framerate 4 -pattern_type glob -i "${png_output_dir}/${prefix}_*.png" -c:v libx264 -pix_fmt yuv420p "${mp4_output_dir}/${prefix}.mp4"
done
