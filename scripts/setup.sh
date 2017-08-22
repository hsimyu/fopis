#!/bin/sh
INSTALL_LOCATION=$HOME/local
HDF5_VER=1.10.1
HDF5_URL=https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-${HDF5_VER}.tar.bz2

# Download submodules from github.com
git submodule init
git submodule update

# -- install libraries --
mkdir libraries
cd libraries

# install zlib
wget http://zlib.net/zlib-1.2.11.tar.gz
tar zxvf zlib-1.2.11.tar.gz
cd zlib-1.2.11
./configure --prefix=${INSTALL_LOCATION}
make
make install
cd ..

# install hdf5
wget ${HDF5_URL}
tar zxvf hdf5-${HDF5_VER}.tar.bz2
cd hdf5-${HDF5_VER}
./configure --enable-build-mode=production --enable-cxx --prefix=${INSTALL_LOCATION} --with-zlib=${INSTALL_LOCATION}/include,${INSTALL_LOCATION}/lib
make
make install
cd ..

cd ..
rm -rf ./libraries

mkdir build
cd build
cmake ..
