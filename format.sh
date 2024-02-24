#!/bin/bash

# on pace cluster, run these 2 commands to install fprettify
# ml python/3.9.12-rkxvr6
# pip3 install --upgrade fprettify

# run bash ./format.sh in RBC3D directory

echo "formatting example cases and common directory"

echo "$PWD"
fprettify ./examples -r --indent 2
fprettify ./common -r --indent 2

echo "done"