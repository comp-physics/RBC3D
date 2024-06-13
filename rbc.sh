#!/bin/bash

# change this to the correct module load for python3 on your cluster
ml python/3.9.12-rkxvr6

if [[ ":$PATH:" != *":$HOME/.local/bin:"* ]]; then
    export PATH="$PATH:$HOME/.local/bin"
fi

# check python3 works
command -v python3 > /dev/null 2>&1
if (($?)); then
    echo "[rbc.sh] Error: Couldn't find Python3. Please ensure it is discoverable."
    exit 1
fi

python3 -c 'print("")' > /dev/null
if (($?)); then
    echo "[rbc.sh] Error: Python3 is present but can't execute a simple program. Please ensure that python3 is working."
    exit 1
fi

if [ "$1" == 'format' ]; then
    . "$(pwd)/install/format.sh" $@; exit
fi

if [ "$1" == 'install-with-mkl' ]; then
    . "$(pwd)/install/install-with-mkl.sh" $@; exit
fi

if [ "$1" == 'install-with-lapack' ]; then
    echo "HERE"
    . "$(pwd)/install/install-with-lapack.sh" $@; exit
fi

if [ "$1" == 'install-no-petsc' ]; then
    . "$(pwd)/install/install-no-petsc.sh" $@; exit
fi
