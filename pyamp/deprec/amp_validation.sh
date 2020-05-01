#!/bin/bash

if [ $1 == "-h" ] ; then
    echo "Usage:: amp_validation.sh Ethylene.extxyz '3 3 3' 0.001 -g "
    echo "Usage for run:: amp_validation.sh Ethylene.extxyz '3 3 3' 0.001 -g | sh"
    exit
fi

set PYTHON = "$HOME/anaconda3/bin/python"

# ethylene validataion check
for var in 0 1 2 3
do
    echo "$PYTHON amp_ene.py $1 val -n 5  -hl $2 -el $3  -i $var $4"
done

