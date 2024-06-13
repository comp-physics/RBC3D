#!/bin/bash

set -e -x

base_dir="examples"

cases=$(find "$base_dir" -maxdepth 1 -type d ! -path "$base_dir")

for case in $cases; do
    make -C "$case" .depend
    make -C "$case"
done

echo "`pwd`/lapack-3.11/liblapack.a"