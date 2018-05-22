#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Error: need 1 arg for suffix of rem"
    exit
fi

rem=$1  #rpbed-opt	# b3lyp


for x in `ls *.mol`
  do
    f_pre=$( basename $x | cut -d'.' -f1 )
    cat $x rem.$rem > $f_pre.in
  done

