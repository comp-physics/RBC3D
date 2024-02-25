#!/bin/bash

# on pace cluster, run these commands to install fprettify
# ml python/3.9.12-rkxvr6
# make sure python bin in path, 1 way to do this is to add to your bashrc
# export PATH="$PATH:$HOME/.local/bin"
# pip3 install --upgrade fprettify

# run bash ./format.sh in RBC3D directory

echo "formatting example cases and common directory"

echo "$PWD"
fprettify ./examples -r --indent 2
fprettify ./common -r --indent 2

echo "done"