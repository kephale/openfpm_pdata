#! /bin/bash

source script/detect_gcc
source script/discover_os

detect_gcc_or_clang g++
discover_os

# check if the directory $1/SUITESPARSE exist

if [ -d "$1/SUITESPARSE" ]; then
  echo "SUITESPARSE already installed"
  exit 0
fi

wget http://ppmcore.mpi-cbg.de/upload/SuiteSparse-4.4.5.tar.gz
rm -rf SuiteSparse
tar -xf SuiteSparse-4.4.5.tar.gz
if [ $? != 0 ]; then
  echo "Fail to download SuiteSparse"
  exit 1
fi
cd SuiteSparse

# configuration

if [ x"$platform" = x"osx"  ]; then
    # installation for OSX

    sed -i "" -e "s|INSTALL_LIB = \/usr\/local\/lib|INSTALL_LIB = "$1"\/SUITESPARSE\/lib|" SuiteSparse_config/SuiteSparse_config_Mac.mk
    sed -i "" -e "s|INSTALL_INCLUDE = \/usr\/local\/include|INSTALL_INCLUDE = "$1"\/SUITESPARSE\/include|" SuiteSparse_config/SuiteSparse_config_Mac.mk
    sed -i "" -e "s| LAPACK = -llapack|LAPACK = |" SuiteSparse_config/SuiteSparse_config_Mac.mk
    sed -i "" -e "s| BLAS = -lopenblas|BLAS = -L"$1"/OPENBLAS/lib -lopenblas|" SuiteSparse_config/SuiteSparse_config_Mac.mk

    ### Overwrite SuiteSparse_config.mk

    rm SuiteSparse_config/SuiteSparse_config.mk
    mv SuiteSparse_config/SuiteSparse_config_Mac.mk SuiteSparse_config/SuiteSparse_config.mk

else
    # Installation for linux

    sed -i "/INSTALL_LIB\s=\s\/usr\/local\/lib/c\INSTALL_LIB = $1\/SUITESPARSE\/lib" SuiteSparse_config/SuiteSparse_config.mk
    sed -i "/INSTALL_INCLUDE\s=\s\/usr\/local\/include/c\INSTALL_INCLUDE = $1\/SUITESPARSE\/include" SuiteSparse_config/SuiteSparse_config.mk
    sed -i "/\sLAPACK\s=\s-llapack/c\LAPACK = " SuiteSparse_config/SuiteSparse_config.mk
    sed -i "/\sBLAS\s=\s\-lopenblas/c\BLAS = -L$1/OPENBLAS/lib -lopenblas" SuiteSparse_config/SuiteSparse_config.mk

fi

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$1/OPENBLAS/lib"

make
if [ $? != 0 ]; then
  echo "Fail to compile SuiteSparse"
  exit 1
fi
mkdir $1/SUITESPARSE
mkdir $1/SUITESPARSE/lib
mkdir $1/SUITESPARSE/include
make install
rm -rf SuiteSparse
rm SuiteSparse-4.4.5.tar.gz
