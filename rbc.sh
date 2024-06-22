#!/bin/bash

# module load the latest version of python3 on your cluster

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

if [ "$1" == 'install' ]; then
    . "$(pwd)/install/install.sh" $@; exit
fi

if [ "$1" == 'install-phoenix' ]; then
    . "$(pwd)/install/install-phoenix.sh" $@; exit
fi

if [ "$1" == 'install-makedepf90' ]; then
    . "$(pwd)/install/install-makedepf90.sh" $@; exit
fi

if [ "$1" == 'cmake' ]; then
    . "$(pwd)/install/cmake.sh" $@; exit
fi
