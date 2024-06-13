#!/bin/bash

# this installer script assumes you module loaded petsc 
# so it doesn't install it

mkdir packages
cd packages

# build and install spherepack
git clone https://github.com/comp-physics/spherepack3.2.git
cd spherepack3.2
make -j 8

# build and install makedepf90
cd ..
git clone https://github.com/comp-physics/makedepf90.git
cd ../install/scripts
echo "PWD $(pwd)"
python3 mdf90_replace.py
cd ../../packages/makedepf90
make -j 8
make install

echo "Done installing RBC3D!"