#!/bin/sh
#
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
./configure --prefix=$HOME/local
make
make install
cd ..

# install hdf5
wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.0-patch1.tar.bz2
tar zxvf hdf5-1.10.0-patch1.tar.bz2
cd hdf5-1.10.0-patch1
./configure --enable-build-mode=production --prefix=$HOME/local
make
make install
cd ..

# install silo
wget https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2.tar.gz
tar zxvf silo-4.10.2.tar.gz
cd silo-4.10.2
./configure --prefix=$HOME/local --with-hdf5=$HOME/local/include,$HOME/local/lib --with-zlib=$HOME/local/include,$HOME/local/lib
make
make install
cd ..

cd ..
rm -rf ./libraries

mkdir build
cd build
cmake ..
