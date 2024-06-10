#!/bin/bash

if [ "$1" == 'format' ]; then
    . "$(pwd)/scripts/format.sh" $@; exit
fi