#!/usr/bin/env python

import os
import fileinput
import sys
import subprocess


def replace(file, old, new):
    for line in fileinput.input(file, inplace=1):
        if old in line:
            line = line.replace(old, new)
        sys.stdout.write(line)


curdir = os.getcwd()
target = f'{curdir}/Input/tube.in'
old = "'D/restart.LATEST.dat'"

for i in range(21500, 25700, 100):
    numZeros = 9 - len(str(i))
    print(f'i: {i}, numZeros: {numZeros}')
    numStr = ("0"*numZeros) + str(i)
    filePath = "'R/" + "restart" + numStr + ".dat'"
    print(filePath)
    replace(target, old, filePath)
    old = filePath
    
    subprocess.run(["srun", "-n", "1", "../../build/field_visual/traction"])
    
    src = "./D/wallvectors.csv"
    dst = f"./D/csvs/wallvectors{numStr}.csv"
    os.rename(src, dst)
