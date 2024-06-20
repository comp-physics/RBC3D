#!/bin/bash

# assume packages directory is already created
cd packages

# build and install makedepf90
cd ..
git clone https://github.com/comp-physics/makedepf90.git
cd ../install/scripts
python3 mdf90_replace.py
cd ../../packages/makedepf90
make -j 8
make install

echo "Done installing makedepf90"