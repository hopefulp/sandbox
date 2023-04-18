#!/bin/bash

if [ $# -lt 2 ]; then
    echo "input two numbers: \$1 for start process \$2 for the last process"
    exit 1
fi

i=$1
n=$(expr $2 + 1)


list=()
host=$(hostname)

while [ $i -lt $n ]; do
    echo $qdel $i
    qdel $i
    i=$(expr $i + 1)
    done

