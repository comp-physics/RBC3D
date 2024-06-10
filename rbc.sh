#!/bin/bash

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