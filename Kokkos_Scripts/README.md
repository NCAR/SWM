# Installing Kokkos

The install_kokkos-*.sh scripts can be used to install a unique build of kokkos

They can be executed from any directory.

The line defining SRC_DIR must be modified to the path to where you have cloned the kokkos repository (https://github.com/kokkos/kokkos.git).


# Generating Executables

The buildRun*.sh scripts can be used to generate an executable using a particular kokkos install. They can be executed from any directory. And they should be executed using **source**, so the exported environment variables persist.

    source buildRun*.sh

You have to use the script that contains in its filename the same compiler command that was used to build the particular kokkos install that you would like to be used to compile your source code (hereafter referred to as "the kokkos install").

The variable **install_name** must be modified match the name of the directory containing the kokkos install.

The variable **project_dir** must be modified to be the path to the directory containing the source code that you would like to compile.

The environment variable **Kokkos_DIR** must be modified to be the path to the kokkos install.


# Setting up Kokkos Profiler

The setup_kokkos_profiler.sh script can be used to prepare the Kokkos simple_kernel_timer profiler, so whenever you run kokkos program, a *.dat file will be generated that can be viewed using the following command.

    kp_reader *.dat

The setup_kokkos_profiler.sh script should be executed using **source**, so the modification to your PATH environment variable will persist.

    source setup_kokkos_profiler.sh

The variable **profiler_dir** must be modified to be the path to the simple_kernel_timer tool in your cloned kokkos-tools repository.