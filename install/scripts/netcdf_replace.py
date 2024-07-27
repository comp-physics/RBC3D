#!/usr/bin/env python

import os
import fileinput
import sys

def addnewline(file, old, new):
    for line in fileinput.input(file, inplace=True):
        if old in line:
            sys.stdout.write(new + '\n')
        sys.stdout.write(line)

def replace(file, old, new):
    for line in fileinput.input(file, inplace=1):
        if old in line:
            line = line.replace(old, new)
        sys.stdout.write(line)

curdir = os.path.dirname(os.path.dirname(os.getcwd()))
targetdir = f'{curdir}/packages/srcNETCDF/netcdf-fortran-4.6.1/fortran'
print(f"target: {targetdir}")
replace(f"{targetdir}/netcdf.F90", "f90", "F90")
replace(f"{targetdir}/netcdf.F90", "_subset", "")
replace(f"{targetdir}/netcdf.F90", "#include \"netcdf_get_nd_expanded.F90\"", "")

addnewline(f"{targetdir}/module_netcdf_nc_data.F90", "#ifdef USE_NETCDF4", "#define USE_NETCDF4")