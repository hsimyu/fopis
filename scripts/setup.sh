#!/bin/sh

# Download submodules from github.com
git submodule init
git submodule update
#wget hdf5

mkdir libraries
cd libraries

# install silo
wget https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2.tar.gz
tar zxvf silo-4.10.2.tar.gz
cd silo-4.10.2
./configure --prefix=$HOME/local
make
make install

cd ..
rm -rf ./libraries

mkdir build
cd build
cmake ..
