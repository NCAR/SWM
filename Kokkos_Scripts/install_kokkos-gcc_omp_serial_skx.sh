#!/bin/bash

module load cmake
module load gnu/8.3.0
module list

SRC_DIR=/glade/work/$USER/kokkos

[ ! -d "$SRC_DIR" ] && exit 1

BLDS_DIR=$SRC_DIR/blds
BLD_DIR=$SRC_DIR/blds/gcc_omp_serial_skx
INSTALLS_DIR=$SRC_DIR/installs
INSTALL_DIR=$SRC_DIR/installs/gcc_omp_serial_skx

[ ! -d "$BLDS_DIR" ] && mkdir $BLDS_DIR
[ -d "$BLD_DIR" ] && rm -r $BLD_DIR
mkdir $BLD_DIR

[ ! -d "$INSTALLS_DIR" ] && mkdir $INSTALLS_DIR
[ -d "$INSTALL_DIR" ] && rm -r $INSTALL_DIR
mkdir $INSTALL_DIR

cd $BLD_DIR

cmake $SRC_DIR \
  -DCMAKE_CXX_COMPILER=gcc \
  -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
  -DCMAKE_CXX_EXTENSIONS=On \
  -DKokkos_ARCH_SKX=On \
  -DKokkos_ENABLE_OPENMP=On \
  -DKokkos_ENABLE_SERIAL=On

make install VERBOSE=1