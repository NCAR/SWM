# Installing Kokkos

The install_kokkos-*.sh scripts can be used to install a unique build of kokkos

They can be executed from any directory.

The line defining SRC_DIR must be modified to the path to where you have cloned the kokkos repository (https://github.com/kokkos/kokkos.git).

# Generating Executables

The generate_executable-*.sh scripts can be used to generate an executable using a particular kokkos install.

They can be executed from any directory.

You have to use the script that contains in its filename the same compiler command that was used to build the particular kokkos install that you would like to be used to compile your source code (hereafter referred to as "the kokkos install").

The environment variable **install_name** must be modified match the name of the directory containing the kokkos install.

    export install_name=<install directory name>

The environment variable **project_dir** must be modified to be the path to the directory containing the source code that you would like to compile.

    export project_dir=<path to code directory>

The environment variable Kokkos_DIR must be modified to be the path to the kokkos install.

    export Kokkos_DIR=<path to install directory>
