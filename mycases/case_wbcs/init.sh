#!/bin/bash

echo "This is a message."
cd D
rm -rf *
cd ../

make clean
make .depend
make
srun -n 2 ./initcond