# Installing Kokkos

The install_kokkos-*.sh scripts can be used to install a unique build of kokkos

They can be executed from any directory.

The line defining SRC_DIR must be modified to the path to where you have cloned the kokkos repository (https://github.com/kokkos/kokkos.git).


# Building and Running Executables

The buildRun*.sh scripts can be used to generate an executable using a particular kokkos install. They can be executed from any directory. And they should be executed using **source**, so the exported environment variables persist.

    source buildRun*.sh [-s [-i id [-a] [-m 64 [-n 128]]]]

Use the following command to see a description of the available build and run options.

    ./buildRun*.sh -h

The variable **install_name** should match the name of the directory containing the kokkos install.

The variable **project_dir** should be the path to the directory containing the source code that you would like to compile.

The environment variable **Kokkos_DIR** should be the path to the kokkos install.

The variable **build_name** should be the desired name for the directory in which results will be stored.


# Setting up Kokkos Profiler

The setup_kokkos_profiler.sh script can be used to prepare the Kokkos simple_kernel_timer profiler, so whenever you run kokkos program, a *.dat file will be generated that can be viewed using the following command.

    kp_reader *.dat

The setup_kokkos_profiler.sh script should be executed using **source**, so the modification to your PATH environment variable will persist.

    source setup_kokkos_profiler.sh

The variable **profiler_dir** must be modified to be the path to the simple_kernel_timer tool in your cloned kokkos-tools repository.
