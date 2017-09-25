#!/bin/sh
INSTALL_LOCATION=$HOME/local
MPI_COMPILER=/usr/local/bin/mpicc
HDF5_VER=1.10.1
HDF5_URL=https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-${HDF5_VER}.tar.gz

# Download submodules from github.com
git submodule init
git submodule update

# -- install libraries --

if [ ! -d ./libraries ]; then
    mkdir libraries
fi

cd libraries

# install zlib
if [ ! -f ./zlib-1.2.11.tar.gz ]; then
    wget http://zlib.net/zlib-1.2.11.tar.gz
    tar zxf zlib-1.2.11.tar.gz
fi

cd zlib-1.2.11
./configure --prefix=${INSTALL_LOCATION}
make
make install
cd ..

# install hdf5
if [ ! -f ./hdf5-{HDF5_VER}.tar.gz ]; then
    wget ${HDF5_URL}
    tar zxf hdf5-${HDF5_VER}.tar.gz
fi

cd hdf5-${HDF5_VER}
CC=${MPI_COMPILER} ./configure --prefix=${INSTALL_LOCATION} --disable-shared --enable-build-mode=production --enable-parallel --with-zlib=${INSTALL_LOCATION}/include,${INSTALL_LOCATION}/lib
make
make install
cd ..

cd ..

mkdir build
cd build
cmake .. -DLOCAL_LIBRARYDIR=${INSTALL_LOCATION} -DCMAKE_PREFIX_PATH=""
