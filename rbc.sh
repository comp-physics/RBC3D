#!/bin/bash

load_modules() {
    # A list with all the necessary modules
    local -a modules=('gcc'
                        'mvapich2'
                        'mkl'
                        'netcdf-c'
                        'netcdf-cxx'
                        'netcdf-fortran'
                        'fftw')
    local mod

    echo "Running function"
    echo

    # Check if modules are loaded. If not, load them
    for mod in "${modules[@]}"
    do
        if module is-loaded "$mod"
        then
            echo "$mod is already loaded"
        else
            echo "Loading $mod…"
            module load $mod
        fi
    done

    echo
    echo "Finished loading modules on Phoenix!"
}

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
    . "$(pwd)/install/install-with-lapack.sh" $@; exit
fi

if [ "$1" == 'phoenix-modules' ]; then
    . "$(pwd)/install/phoenix-modules.sh" $@
    bash; exit
    # . "$(pwd)/install/phoenix-modules.sh" $@; exit
fi