#!/bin/bash

# run this in RBC3D directory

ml python/3.9.12-rkxvr6
# make sure fprettify installed via this command
# pip3 install --upgrade fprettify

echo "formatting example cases and common directory"

echo "$PWD"
fprettify ./examples -r --indent 2
fprettify ./common -r --indent 2

echo "done"