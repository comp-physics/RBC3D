#!/usr/bin/env python

import os
import fileinput
import sys

def replace(file, old, new):
    for line in fileinput.input(file, inplace=1):
        if old in line:
            line = line.replace(old, new)
        sys.stdout.write(line)

curdir = os.path.dirname(os.path.dirname(os.getcwd()))
targetdir = f'{curdir}/packages/makedepf90'
print(f"target: {targetdir}")
replace(f"{targetdir}/Makefile", "PREFIX		?= /usr/local", f"PREFIX		?= {targetdir}")