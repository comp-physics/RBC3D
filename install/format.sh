#!/bin/bash

echo "[rbc.sh] Formatting example cases and common directory in $PWD"

pip3 install --upgrade fprettify

fprettify ./examples -r --indent 2
fprettify ./common -r --indent 2

echo "[rbc.sh] Done formatting RBC3D!"