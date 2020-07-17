#!/bin/bash

if [ $# -lt 2 ]; then
    echo "input two numbers"
    exit 1
fi

i=$1
n=$2

while [ $i -ne $n ]; do
    echo qdel $i
    qdel $i
    i=$(expr $i + 1)
    done

