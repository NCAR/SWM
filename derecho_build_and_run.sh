#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status
set -u  # Treat expanding empty variables as an error
set -x  # Print each command before executing it

###############################################################################
# User Input
###############################################################################

# Directory where you pulled the SWM repository
export SWM_ROOT=/glade/u/home/htorres/SWM

# Directory where the SWM build will be created
# This is used as a base name and will appended to based on the options you pick
export SWM_BUILD_DIR_BASE=/glade/u/home/htorres/SWM_build

# Directory where you pulled the AMReX repositry
export AMREX_HOME=/glade/u/home/htorres/amrex

# Directory where the amrex AMReX build will be created
# This is used as a base name and will appended to based on the options you pick
export AMREX_BUILD_DIR_BASE=/glade/u/home/htorres/amrex_build

# Set to GNU or NVHPC
#export COMPILER=GNU   
#export COMPILER=NVHPC  

# Set to "cpu", "gpu"
#export SWM_DEVICE=cpu  

# Save the base directory for the SWM build
#for SWM_DEVICE in cpu gpu; do
for SWM_DEVICE in gpu; do
  #for COMPILER in GNU NVHPC; do
  #for COMPILER in GNU; do
  for COMPILER in NVHPC; do


# All other options are ON or OFF
export SWM_C=OFF        
export SWM_FORTRAN=OFF
export SWM_AMREX=ON

export SWM_ACC=ON        
export SWM_MPI=OFF        
export SWM_OMP=OFF        
export SWM_CUDA=ON


###############################################################################
# Setup based on user input
###############################################################################

# If you are building for the cpu then turn off the CUDA option... even if the user set it to ON
if [[ "${SWM_DEVICE}" == "cpu" ]] && [[ "${SWM_CUDA}" == "ON" ]]; then
    export SWM_CUDA=OFF
fi

# These build directories will be appended to based on the user input
export SWM_BUILD_DIR="${SWM_BUILD_DIR_BASE}_${SWM_DEVICE}_${COMPILER}"
export AMREX_BUILD_DIR="${AMREX_BUILD_DIR_BASE}_${COMPILER}"

###############################################################################
# Module Setup
###############################################################################

module purge

# Modules we always use
module load cmake

if [[ "${COMPILER}" == "GNU" ]]; then
    module load gcc
elif [[ "${COMPILER}" == "NVHPC" ]]; then
    module load nvhpc
else
    echo "Unsupported compiler option: ${COMPILER}"
    exit 1
fi

if [[ "${SWM_AMREX}" == "ON" ]]; then

    if [[ "${SWM_MPI}" == "ON" ]]; then
        module load cray-mpich
        #export SWM_BUILD_DIR="${SWM_BUILD_DIR}_MPI"
        export AMREX_BUILD_DIR="${AMREX_BUILD_DIR}_MPI"
    fi

    if [[ "${SWM_OMP}" == "ON" ]]; then
        #export SWM_BUILD_DIR="${SWM_BUILD_DIR}_OMP"
        export AMREX_BUILD_DIR="${AMREX_BUILD_DIR}_OMP"
    fi

    if [[ "${SWM_DEVICE}" == "gpu" ]] && [[ "${SWM_CUDA}" == "ON" ]]; then
        module load cuda                                  # If we are using the nvhpc module then that already load the cuda module, but this is good to be explicit
        #export SWM_BUILD_DIR="${SWM_BUILD_DIR}_CUDA"
        export AMREX_BUILD_DIR="${AMREX_BUILD_DIR}_CUDA"
    fi

fi

if [[ "${COMPILER}" == "GNU" ]]; then
    module load ncarcompilers
fi

# HDF5 is only used by the AMReX mini-app version of SWM, still loading it by default for now but will eventually move to netcdf output
module load hdf5
# Loading the hdf5 module on derecheo should define a variable NCAR_ROOT_HDF5. We need to set HDF5_HOME.
export HDF5_HOME="${NCAR_ROOT_HDF5}"

# Profile or Debugg with linary forge
#module load linaro-forge

module list


# WARNING: YES will delete the build directory if it exists. Make sure this wont delete anything important!
fresh_build_amrex=NO
if [[ "${fresh_build_amrex}" == "YES" ]]; then
    if [[ -d "${AMREX_BUILD_DIR}" ]]; then
        echo "Deleting existing AMReX build directory: ${AMREX_BUILD_DIR}"
        rm -rf "${AMREX_BUILD_DIR}"
    fi
fi

# WARNING: YES will delete the build directory if it exists. Make sure this wont delete anything important!
fresh_build_swm=YES
if [[ "${fresh_build_swm}" == "YES" ]]; then
    if [[ -d "${SWM_BUILD_DIR}" ]]; then
        echo "Deleting existing SWM build directory: ${SWM_BUILD_DIR}"
        rm -rf "${SWM_BUILD_DIR}"
    fi
fi

###############################################################################
# Build the version of AMReX that we are asking fori

###############################################################################
# Build the version of AMReX that we are asking for
###############################################################################

if [[ "${SWM_AMREX}" == "ON" ]]; then

    # Initialize an array for CMake AMReX build options
    amrex_cmake_opts=()
    
    # AMReX options we always-use
    amrex_cmake_opts+=("-DAMReX_SPACEDIM=2")
    amrex_cmake_opts+=("-DAMReX_PRECISION=DOUBLE")
    amrex_cmake_opts+=("-DAMReX_FORTRAN=YES")
    amrex_cmake_opts+=("-DCMAKE_C_COMPILER=$CC")
    amrex_cmake_opts+=("-DCMAKE_CXX_COMPILER=$CXX")
    amrex_cmake_opts+=("-DCMAKE_Fortran_COMPILER=$FC")
    
    #if [[ "${COMPILER}" == "NVHPC" ]] && [[ "${AMREX_USE_MPI}" == "YES" ]]; then
    #    #echo "Building AMReX with NVHPC, MPI, and CUDA is not working in derecho yet. Please use GNU compiler instead."
    #    #exit 1
    #    amrex_cmake_opts+=("-DCMAKE_C_COMPILER=mpicc")
    #    amrex_cmake_opts+=("-DCMAKE_CXX_COMPILER=mpicxx")
    #    amrex_cmake_opts+=("-DCMAKE_Fortran_COMPILER=mpifort")
    #    #amrex_cmake_opts+=("-DCMAKE_CUDA_COMPILER=nvcc")
    #else
    #    amrex_cmake_opts+=("-DCMAKE_C_COMPILER=$CC")
    #    amrex_cmake_opts+=("-DCMAKE_CXX_COMPILER=$CXX")
    #    amrex_cmake_opts+=("-DCMAKE_Fortran_COMPILER=$FC")
    #fi
    
    # These AMReX options were on by default, however I don't think we are using any of these features so I am turning them off.
    amrex_cmake_opts+=("-DAMReX_LINEAR_SOLVERS=NO")
    amrex_cmake_opts+=("-DAMReX_LINEAR_SOLVERS_INCFLO=NO")
    amrex_cmake_opts+=("-DAMReX_LINEAR_SOLVERS_EM=NO")
    amrex_cmake_opts+=("-DAMReX_AMRLEVEL=NO")
    amrex_cmake_opts+=("-DAMReX_PARTICLES=NO")
    amrex_cmake_opts+=("-DAMReX_TINY_PROFILE=NO")
    
    # TODO: Decide if we want to use these options or just base the logic on the SWM ones
    #if [[ "${SWM_DEVICE}" == "cpu" ]]; then
    #  # If we are building for the CPU then just take whatever the user set for the SWM options
    #  export AMREX_USE_MPI=$SWM_MPI   # Set to YES or NO
    #  export AMREX_USE_OMP=YES   # Set to YES or NO
    #  export AMREX_USE_CUDA=YES  # Set to YES or NO
    #elif [[ "${SWM_DEVICE}" == "gpu" ]]; then
    #  # If we are building for the GPU then do somethinge else?
    #  export AMREX_USE_MPI=YES   # Set to YES or NO
    #  export AMREX_USE_OMP=YES   # Set to YES or NO
    #  export AMREX_USE_CUDA=YES  # Set to YES or NO
    #fi

    ## These AMReX options take on differnent values depending on the user input
    #if [[ "${AMREX_USE_MPI}" == "YES" ]]; then
    #    amrex_cmake_opts+=("-DAMReX_MPI=YES")
    #else
    #    amrex_cmake_opts+=("-DAMReX_MPI=NO")
    #fi
    
    #if [[ "${AMREX_USE_OMP}" == "YES" ]]; then
    #    amrex_cmake_opts+=("-DAMReX_OMP=YES")
    #else
    #    amrex_cmake_opts+=("-DAMReX_OMP=NO")
    #fi
    
    #if [[ "${AMREX_USE_CUDA}" == "YES" ]]; then
    #    amrex_cmake_opts+=("-DAMReX_GPU_BACKEND=CUDA")
    #    amrex_cmake_opts+=("-DAMReX_GPU_RDC=YES")
    #    amrex_cmake_opts+=("-DAMReX_CUDA_ARCH=8.0") # Set the CUDA architecture version, adjust as needed
    #    amrex_cmake_opts+=("-DAMReX_DIFFERENT_COMPILER=ON")
    #else
    #    amrex_cmake_opts+=("-DAMReX_GPU_BACKEND=NONE")
    #fi

    # AMReX options based on the user input... had to translate from our ON and OFF to what AMReX uses, YES and NO
    if [[ "${SWM_MPI}" == "ON" ]]; then
        amrex_cmake_opts+=("-DAMReX_MPI=YES")
    else
        amrex_cmake_opts+=("-DAMReX_MPI=NO")
    fi

    if [[ "${SWM_OMP}" == "ON" ]]; then
        amrex_cmake_opts+=("-DAMReX_OMP=YES")
    else
        amrex_cmake_opts+=("-DAMReX_OMP=NO")
    fi

    if [[ "${SWM_CUDA}" == "ON" ]]; then
        amrex_cmake_opts+=("-DAMReX_GPU_BACKEND=CUDA")
        amrex_cmake_opts+=("-DAMReX_GPU_RDC=YES")
        amrex_cmake_opts+=("-DAMReX_CUDA_ARCH=8.0") 
        amrex_cmake_opts+=("-DAMReX_DIFFERENT_COMPILER=ON")
    else
        amrex_cmake_opts+=("-DAMReX_GPU_BACKEND=NONE")
    fi
    
    # Installing AMReX as a subdirectory of the AMReX build directory. Could change this to a different location if needed.
    export AMREX_INSTALL_DIR=$AMREX_BUILD_DIR/install
    
    mkdir -p $AMREX_BUILD_DIR
    cd $AMREX_BUILD_DIR
    
    cmake "${amrex_cmake_opts[@]}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$AMREX_INSTALL_DIR" \
        -S "$AMREX_HOME" \
        -B "$AMREX_BUILD_DIR"
    #    --trace \
    
    make -j 32 install 
    #make -j 32 VERBOSE=1 install 
    
    #make test_install  # optional step to test if the installation is working
    #exit 0 # Exit early for testing purposes

fi

###############################################################################
# Build SWM Using the version of AMReX that we just built
###############################################################################

# Initialize an array for SWM CMake build options
swm_cmake_opts=()

swm_cmake_opts+=("-DSWM_DEVICE=${SWM_DEVICE}")

swm_cmake_opts+=("-DSWM_C=${SWM_C}")
swm_cmake_opts+=("-DSWM_FORTRAN=${SWM_FORTRAN}")
swm_cmake_opts+=("-DSWM_AMREX=${SWM_AMREX}")

if [[ "${SWM_AMREX}" == "ON" ]]; then
    swm_cmake_opts+=("-DAMReX_ROOT=$AMREX_INSTALL_DIR/lib/cmake/AMReX")
fi

swm_cmake_opts+=("-DSWM_OPENACC=${SWM_ACC}")
swm_cmake_opts+=("-DSWM_OPENMP=${SWM_OMP}")
swm_cmake_opts+=("-DSWM_MPI=${SWM_MPI}")
swm_cmake_opts+=("-DSWM_CUDA=${SWM_CUDA}")

swm_cmake_opts+=("-DCMAKE_C_COMPILER=$CC")
swm_cmake_opts+=("-DCMAKE_CXX_COMPILER=$CXX")
swm_cmake_opts+=("-DCMAKE_Fortran_COMPILER=$FC")

# Hard coded options for the SWM build
#swm_cmake_opts+=("-DSWM_DEVICE=cpu;gpu")
#
#swm_cmake_opts+=("-DSWM_C=ON")
#swm_cmake_opts+=("-DSWM_FORTRAN=ON")
#swm_cmake_opts+=("-DSWM_AMREX=ON")
#swm_cmake_opts+=("-DAMReX_ROOT=$AMREX_INSTALL_DIR/lib/cmake/AMReX")
#
#swm_cmake_opts+=("-DSWM_OPENACC=ON")
#swm_cmake_opts+=("-DSWM_OPENMP=ON")
#swm_cmake_opts+=("-DSWM_MPI=ON")
#swm_cmake_opts+=("-DSWM_CUDA=ON")
#
#swm_cmake_opts+=("-DCMAKE_C_COMPILER=$CC")
#swm_cmake_opts+=("-DCMAKE_CXX_COMPILER=$CXX")
#swm_cmake_opts+=("-DCMAKE_Fortran_COMPILER=$FC")


cmake "${swm_cmake_opts[@]}" -S $SWM_ROOT -B $SWM_BUILD_DIR
#cmake "${swm_cmake_opts[@]}" -S $SWM_ROOT -B $SWM_BUILD_DIR --trace-expand

#cmake -DAMReX_ROOT=$AMREX_INSTALL_DIR/lib/cmake/AMReX \
#      -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_Fortran_COMPILER=$FC \
#      -S $SWM_ROOT -B $SWM_BUILD_DIR

cd $SWM_BUILD_DIR

#make 
#make 2>&1 | tee build.log
make VERBOSE=1 2>&1 | tee build.log
if [ "${PIPESTATUS[0]}" -ne 0 ]; then
    exit 1
fi

###############################################################################
# Run all the versions of SWM mini-app
###############################################################################

echo "Running SWM mini-apps for ${SWM_DEVICE} that were built with ${COMPILER} compiler"
#SWM_C=OFF
#SWM_FORTRAN=OFF
#SWM_AMREX=ON

if [[ "${SWM_C}" == "ON" ]]; then
    echo "Running SWM C mini-apps"

    if [[ "${SWM_DEVICE}" == "cpu" ]]; then
        echo "Running SWM plain C mini-app "
        $SWM_BUILD_DIR/swm_c/c/swm_c
    fi

    if [[ "${SWM_ACC}" == "ON" ]]; then
        echo "Running SWM C mini-apps with OpenACC"
        $SWM_BUILD_DIR/swm_c/c_OpenACC/swm_c_acc
        $SWM_BUILD_DIR/swm_c/c_OpenACC/swm_c_acc_tile
    fi
fi

if [[ "${SWM_FORTRAN}" == "ON" ]]; then
    echo "Running SWM Fortran mini-apps"

    if [[ "${SWM_DEVICE}" == "cpu" ]]; then
        echo "Running SWM plain Fortran mini-apps"
        $SWM_BUILD_DIR/swm_fortran/fortran/swm_fortran
        #$SWM_BUILD_DIR/swm_fortran/fortran_AMReX_proxy/swm_fortran_amrex_driver # This runs very slowly particulary with the NVHPC compiler, so I am commenting it out for now
    fi

    if [[ "${SWM_ACC}" == "ON" ]]; then
        echo "Running SWM Fortran mini-apps with OpenACC"
        $SWM_BUILD_DIR/swm_fortran/fortran_OpenACC/swm_fortran_acc
    fi
fi

if [[ "${SWM_AMREX}" == "ON" ]]; then

  # Order matters for the suffix of the executable name
  exe_suffix=""
  #if [[ "${AMREX_USE_OMP}" == "YES" ]]; then # Maybe comeback to this if we decide to have seperate otions for openMP code we wrote vs that in amrex?
  if [[ "${SWM_OMP}" == "ON" ]]; then
    exe_suffix="${exe_suffix}_OMP"
  fi
  if [[ "${SWM_MPI}" == "ON" ]]; then
    exe_suffix="_MPI"
  fi
  if [[ "${SWM_CUDA}" == "ON" ]]; then
    exe_suffix="${exe_suffix}_CUDA"
  fi

  mpi_launcher=""
  if [[ "${SWM_MPI}" == "ON" ]]; then
    #mpi_launcher="mpirun -np 2"
    mpi_launcher="mpirun -np 1"
  fi

  input_file=$SWM_ROOT/swm_amrex/inputs

  if [[ "${SWM_OMP}" == "ON" ]]; then
    #env OMP_NUM_THREADS=2 ${mpi_launcher} $SWM_BUILD_DIR/swm_amrex/swm_AMReX/swm_amrex${exe_suffix} ${input_file}
    env OMP_NUM_THREADS=1 ${mpi_launcher} $SWM_BUILD_DIR/swm_amrex/swm_AMReX/swm_amrex${exe_suffix} ${input_file}
  else
    ${mpi_launcher} $SWM_BUILD_DIR/swm_amrex/swm_AMReX/swm_amrex${exe_suffix} ${input_file}
  fi
  
  #$SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenMP/swm_amrex_fsubroutine_omp ${input_file}
  $SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenACC/swm_amrex_fsubroutine_acc ${input_file}

fi

#$SWM_BUILD_DIR/swm_c/c/swm_c
#$SWM_BUILD_DIR/swm_c/c_OpenACC/swm_c_acc
#$SWM_BUILD_DIR/swm_c/c_OpenACC/swm_c_acc_tile
#$SWM_BUILD_DIR/swm_fortran/fortran/swm_fortran
#$SWM_BUILD_DIR/swm_fortran/fortran_AMReX_proxy/swm_fortran_amrex_driver
#$SWM_BUILD_DIR/swm_fortran/fortran_OpenACC/swm_fortran_acc

#$SWM_BUILD_DIR/swm_amrex/swm_AMReX/swm_amrex $SWM_ROOT/swm_amrex/inputs
#$SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fkernels/swm_amrex_fkernels $SWM_ROOT/swm_amrex/inputs
#$SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenMP/swm_amrex_fsubroutine_omp $SWM_ROOT/swm_amrex/inputs
#$SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenACC/swm_amrex_fsubroutine_acc $SWM_ROOT/swm_amrex/inputs

#if [[ "${AMREX_USE_MPI}" == "YES" ]]; then
#    mpirun -np 2 $SWM_BUILD_DIR/swm_amrex/swm_AMReX/swm_amrex $SWM_ROOT/swm_amrex/inputs
#    mpirun -np 2 $SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fkernels/swm_amrex_fkernels $SWM_ROOT/swm_amrex/inputs
#    mpirun -np 2 $SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenMP/swm_amrex_fsubroutine_omp $SWM_ROOT/swm_amrex/inputs
#    mpirun -np 2 $SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenACC/swm_amrex_fsubroutine_acc $SWM_ROOT/swm_amrex/inputs
#else
#    $SWM_BUILD_DIR/swm_amrex/swm_AMReX/swm_amrex $SWM_ROOT/swm_amrex/inputs
#    $SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fkernels/swm_amrex_fkernels $SWM_ROOT/swm_amrex/inputs
#    $SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenMP/swm_amrex_fsubroutine_omp $SWM_ROOT/swm_amrex/inputs
#    $SWM_BUILD_DIR/swm_amrex/swm_AMReX_Fsubroutine/OpenACC/swm_amrex_fsubroutine_acc $SWM_ROOT/swm_amrex/inputs
#fi

  done # SWM_DEVICE loop
done # COMPILER loop