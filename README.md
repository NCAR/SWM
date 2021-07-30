SWM is a simplified kernel representing the nonlinear PDE's governing geophysical fluid flow.
Please contact the Point of Contact below before distributing it further, or if any questions or issues arise.

# COPY vs SWAP:

There is an optimization that avoids copies in loop 300 by swapping pointers. COPY can be turned on and off by defining _COPY_ or not, as desired.

# Building and Running Executables

The buildRun*.sh scripts can be used to generate and run an executable for a particular version of the code. They can be executed from any directory.

    ./buildRun*.sh [-s [-i id [-a] [-m 64 [-n 128]]]]

Use the following command to see a description of the available build and run options.

    ./buildRun*.sh -h
    
The variable **project_dir** should be the path to the directory containing the source code that you would like to compile.

The variable **build_name** should be the desired name for the directory in which results will be stored, and it will also be included in the name of the generated executable.

# Validation: 

Compare the results.[].txt file(s) with the corresponding results.[].txt generated using the buildRun.cpu.sh file, particularly the last values printed, namely those below-

    diagonal elements of p
    diagonal elements of u
    diagonal elements of v
 
Additionally you can use one of the following command to print the L infinity norm of the difference between two csv output files to the console.

    python compare_results.py -r path/to/ref/results.csv -t path/to/test/results.csv
    
Or you can use the following command to print to the console and save to a csv file the L infinity norm of the difference between all corresponding csv output files of the predefined set of problem sizes from the -a option of the buildRun.[].sh files for two particular directories and IDs.

    python compare_results_all.py -r path/to/ref/dir -i id_of_ref_files -t path/to/test/dir -d id_of_test_files

Point of Contact:

Dr. Richard Loft
National Center for Atmospheric Research
loft@ucar.edu

for Kokkos versions
Zephaniah Connell
zephconnell@gmail.com
