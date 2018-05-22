#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Error: need 2 arg for outfile name and opt-methodm"
    exit
fi

f_pre=$( basename $1 | cut -d'.' -f1 )

base=$2  #rpbed-opt	# b3lyp


qgeo_mo.pl $1 > $f_pre-$base.xyz
xyz22mol.pl $f_pre-$base.xyz

