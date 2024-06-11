#!/bin/bash

# ml gcc mvapich2 mkl python/3.9.12-rkxvr6 netcdf-c netcdf-cxx netcdf-fortran fftw

load_modules() {
    # A list with all the necessary modules
    local -a modules=('gcc'
                        'mvapich2'
                        'mkl'
                        'netcdf-c'
                        'netcdf-cxx'
                        'netcdf-fortran'
                        'fftw/3.3.10-mva2-dgx5sz')
    local mod

    echo "Running function"
    echo

    # Check if modules are loaded. If not, load them
    for mod in "${modules[@]}"; do
        if module is-loaded "$mod"
        then
            echo "$mod is already loaded"
        else
            echo "Loading $mod…"
            module load $mod > /dev/null 2>&1

            code=$?
            if [ "$code" != '0' ]; then
                error "Failed to load module $M$element$CR:"

                # Run load again to show error message
                module load "$element"

                return
            fi
        fi
    done

    echo
    echo "Finished loading modules on Phoenix!"
}

shopt -s expand_aliases
load_modules