#!/bin/bash
INSTALL_LOCATION_A=$HOME/local_a
INSTALL_LOCATION_B=$HOME/local_b
hostname_new="camphor.*"
hostname_old="laurel.*"

# Download submodules from github.com
git submodule init
git submodule update

mkdir build
cd build

if [[ `hostname` =~ $hostname_new ]]; then
    module load cray-hdf5-parallel
    cmake .. -DLOCAL_LIBRARYDIR=$INSTALL_LOCATION_A -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment -DCMAKE_PREFIX_PATH=""
elif [[ `hostname` =~ $hostname_old ]]; then
    module load hdf5-parallel
    cmake .. -DLOCAL_LIBRARYDIR=$INSTALL_LOCATION_B -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment -DCMAKE_INSTALL_PREFIX="" -DOLD_CRAY_ENV=ON -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc
else
    echo 'No target hostname was found!.'
fi
