#!/bin/bash

fn=$1
nmolA=$2
nmolB=$3

if [ $# -lt 3 ]; then
    echo "Usage:: $0 fname[-suffix] nmol_A nmol_B"
    echo "fname == dirname"
    exit
fi
if [ ! -f "$fn.in" ]; then
    echo "Should have whole input file of $fn.in"
    exit
fi
if [ -d $fn ]; then
    echo "Error: there is the same directory already"
    exit
else
    mkdir $fn
fi

cp $fn.in $fn
cd $fn
bsse_2mol.pl $fn $2 $3
