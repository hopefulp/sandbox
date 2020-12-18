#!/bin/bash

if [ $# -lt 2 ]; then
    echo "input two numbers: \$1 for start process \$2 for the last process"
    exit 1
fi

i=$1
f=$2

for ps in `ps aux | grep python | awk '{ print $2 }'`; do
    if [ $i -le $ps -a $ps -le $f ]; then
        echo "kill $ps"
    fi
    done

#n=$(expr $2 + 1)

#while [ $i -lt $n ]; do
#    echo qdel $i
#    qdel $i
#    i=$(expr $i + 1)
#    done

