#!/bin/bash

echo "formatting example cases and common directory in $PWD"

pip3 install --upgrade fprettify

fprettify ./examples -r --indent 2
fprettify ./mycases -r --indent 2
fprettify ./testcases -r --indent 2
fprettify ./common -r --indent 2

echo "done"