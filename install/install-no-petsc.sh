#!/bin/bash

# this installer script assumes you module loaded petsc 
# so it doesn't install it

# build and install spherepack
cd ..
git clone https://github.com/comp-physics/spherepack3.2.git
cd spherepack3.2
make

# build and install makedepf90
cd ..
git clone https://github.com/comp-physics/makedepf90.git
cd makedepf90
# it works without setting prefix
make
make install

echo "Done installing RBC3D!"