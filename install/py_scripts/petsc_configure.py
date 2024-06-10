#!/usr/bin/env python
""" Configuration helper for PETSc

Execute from PETSC_DIR

Note that older versions of PETSc require Python 2

Run with -h to see arguments

Adds a thin extra layer around PETSc's configuration script,
to help combine commonly-used combinations of configuration options and
construct meaningful PETSC_ARCH values.

Proceeds by collecting the set of passed arguments and processing them
before sending them to PETSc's configure script.
This logic will likely be somewhat brittle. Always do a sanity check
and look at the options that are actually being sent. This script should
be simple enough to figure out what's going on.
"""

from __future__ import print_function
import sys
import os
import argparse
import re


def main():
    """ Main script logic """
    args, options_in = get_args()
    options = process_args(options_in, args)
    petsc_configure(options, args)


def get_args():
    """ Retrieve custom arguments and remaining arguments"""
    parser = argparse.ArgumentParser(
        description='Compute arguments to pass to PETSc\'s configure script')
    parser.add_argument(
        '--dryrun',
        action="store_true",
        help="don't actually configure")
    parser.add_argument(
        '--mkl-only',
        action="store_true",
        help="Configure with MKL")
    parser.add_argument(
        '--blas-lapack',
        action="store_true",
        help="Configure with blas/lapack")
    args, unknown = parser.parse_known_args()
    return args, unknown

def process_args(options_in, args):
    """ Main logic to create a set of options for PETSc's configure script,
    along with a corresponding PETSC_ARCH string, if required """

    if args.mkl_only:
        return options_for_mkl()

    if args.blas_lapack:
        return options_for_blaslapack()


def options_for_mkl():
    """ Return a custom set of arguments to simply download and build MPICH """
    
    MKL_DIR = os.getenv('MKLROOT')

    options = []
    options.append('--with-cc=mpicc')
    options.append('--with-cxx=mpicxx')
    options.append('--with-fc=mpif90')
    options.append('--with-fortran-datatypes')
    options.append('--with-debugging=0')
    options.append("--COPTFLAGS=-g -O3 -march=native")
    options.append("--CXXOPTFLAGS=-g -O3 -march=native")
    options.append("--FOPTFLAGS=-g -O3 -march=native")
    options.append("--with-blaslapack-dir=" + MKL_DIR)
    options.append("--with-mpiexec=srun")
    # option from original configure but i think this is unneccessary
    # options.append("--with-shared-libraries=0")
    options.append("--with-x11=0")
    options.append("--with-x=0")
    options.append("--with-windows-graphics=0")
    return options

def options_for_blaslapack():
    """ Return a custom set of arguments to simply download and build MPICH """
    
    PACKAGES_DIR = os.path.dirname(os.getcwd())
    LAPACK_LIB = PACKAGES_DIR + "/lapack-3.11/liblapack.a"
    BLAS_LIB = PACKAGES_DIR + "/lapack-3.11/librefblas.a"

    print("PACKAGES_DIR: " + PACKAGES_DIR)
    print("BLAS_LIB: " + BLAS_LIB)
    print("LAPACK_LIB: " + LAPACK_LIB)

    options = []
    options.append('--with-cc=mpicc')
    options.append('--with-cxx=mpicxx')
    options.append('--with-fc=mpif90')
    options.append('--with-fortran-datatypes')
    options.append('--with-debugging=0')
    options.append("--COPTFLAGS=-g -O3 -march=native")
    options.append("--CXXOPTFLAGS=-g -O3 -march=native")
    options.append("--FOPTFLAGS=-g -O3 -march=native")
    options.append("--with-blas-lib=" + BLAS_LIB)
    options.append("--with-lapack-lib=" + LAPACK_LIB)
    options.append("--with-mpiexec=srun")
    # options.append("--with-shared-libraries=0")
    options.append("--with-x11=0")
    options.append("--with-x=0")
    options.append("--with-windows-graphics=0")
    return options

def petsc_configure(options, args):
    """ Standard PETSc configuration script logic (from config/examples) """
    if args.dryrun:
        print("Dry Run. Would configure with these options:")
        print("\n".join(options))
    else:
        sys.path.insert(0, os.path.abspath('config'))
        # print(sys.version)
        try:
            import configure
        except ImportError:
            print(
                'PETSc configure module not found. Make sure you are executing from PETSC_DIR'
            )
            sys.exit(1)
        print('Configuring with these options:')
        print("\n".join(options))

        # Since petsc_configure looks directly at sys.argv, remove and replace arguments
        argv_temp = sys.argv
        sys.argv = sys.argv[:1]
        configure.petsc_configure(options)
        sys.argv = argv_temp


if __name__ == '__main__':
    main()