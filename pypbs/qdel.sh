#!/bin/bash

function cut_dot {
    pre=`cut $1 -d. -f1`
    echo $pre
    }


if [ $# -lt 2 ]; then
    echo "input two numbers: \$1 for start process \$2 for the last process"
    exit 1
fi

list=()

for x in $@; do
    if [[ $x =~ "." ]]; then
        x=`echo $x | cut -d. -f1`
    fi
    list+=($x)
done

i=${list[0]}
n=${list[1]}

echo $i, $n

while [ $i -le $n ]; do
    echo qdel $i
    qdel $i
    i=$(expr $i + 1)
    done

