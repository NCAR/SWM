#!/bin/bash

# This script is used to convert all AMReX plotfiles in a specified directory to HDF5 format
# and then plot the hdf5 files using matplotlib. The script will also create an MP4 movie from the
# PNG files created by the plotting script.

###############################################################################
# User Input
###############################################################################
plotfile2hdf5_exe="$SWM_AMREX_ROOT"/plotting_utils/main2d.gnu.MPI.ex

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

set -u
set -e

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
# Convert to HDF5 and Plot to create PNG files
###############################################################################

# Loop over all files that start with plt_ followed by a series of numbers in the specified directory
for infile in "$input_directory"/plt[0-9]*; do

    # Skip files that contain .old or .h5
    if [[ "$infile" == *".old"* || "$infile" == *".h5"* ]]; then
        continue
    fi

    print_banner "Converting $infile to HDF5 and Plotting"

    # Plotfiles are actually directories. Check that infile is a directory
    if [[ -d "$infile" ]]; then
        # Set the outfile name to the same as infile but with .h5 extension
        hdf5_file="${infile}.h5"
        
        # Run the program with the infile and outfile arguments
        "$plotfile2hdf5_exe" infile="$infile" outfile="$hdf5_file"
    fi

    python "$SWM_AMREX_ROOT"/plotting_utils/plot_hdf5.py "$hdf5_file" --output_dir "$input_directory"

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

# Create the MP4 movie from PNG files
for prefix in u v p; do
    ffmpeg -y -framerate 4 -pattern_type glob -i "${prefix}_*.png" -c:v libx264 -pix_fmt yuv420p "${prefix}.mp4"
done