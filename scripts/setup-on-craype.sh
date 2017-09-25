#!/bin/sh
INSTALL_LOCATION=$HOME/local

# Download submodules from github.com
git submodule init
git submodule update

# load modules
module load cray-hdf5-parallel

mkdir build
cd build
cmake .. -DLOCAL_LIBRARYDIR=${INSTALL_LOCATION} -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment -DCMAKE_PREFIX_PATH=""
