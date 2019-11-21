#!/bin/bash

if [ $# -lt 3 ]; then
    echo "Usage:: $0 inode iprocess nprocess[default=16]"
    exit 1
fi
inode=$1
iprocess=$2
npro=16
nprocess=${3:-$npro}
i=0

while [ $i -lt $nprocess ]; do
    ipro=$(($iprocess + $i))
    echo "ssh  $inode kill $ipro"
    if [ ! -z $4 ]; then
        ssh  $inode kill $ipro
    fi
    i=$(($i + 1))
    done
if [ -z $4 ]; then
    echo "this will run if \$4"
fi
