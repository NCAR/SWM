#!/bin/bash

module load cmake
module load cuda/10.1
module load gnu/7.4.0
module list

SRC_DIR=/glade/work/$USER/kokkos

[ ! -d "$SRC_DIR" ] && exit 1

BLDS_DIR=$SRC_DIR/blds
BLD_DIR=$SRC_DIR/blds/nvcc_cuda_omp_serial_volta70_skx
INSTALLS_DIR=$SRC_DIR/installs
INSTALL_DIR=$SRC_DIR/installs/nvcc_cuda_omp_serial_volta70_skx

[ ! -d "$BLDS_DIR" ] && mkdir $BLDS_DIR
[ ! -d "$BLD_DIR" ] && mkdir $BLD_DIR
[ ! -d "$INSTALLS_DIR" ] && mkdir $INSTALLS_DIR
[ ! -d "$INSTALL_DIR" ] && mkdir $INSTALL_DIR

cd $BLD_DIR

cmake $SRC_DIR \
  -DCMAKE_CXX_COMPILER=$SRC_DIR/bin/nvcc_wrapper \
  -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
  -DCMAKE_CXX_EXTENSIONS=On \
  -DKokkos_ARCH_VOLTA70=On \
  -DKokkos_ARCH_SKX=On \
  -DKokkos_ENABLE_CUDA=On \
  -DKokkos_ENABLE_CUDA_LAMBDA=On \
  -DKokkos_ENABLE_OPENMP=On \
  -DKokkos_ENABLE_SERIAL=On

make install