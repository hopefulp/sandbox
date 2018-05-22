#!/bin/bash

mol=$1
mtmp=rimp2
jtmp=opt

if [ $# -lt 1 ]; then
    echo "Usage:: $0 fname[-suffix]"
    echo "    default:: $mtmp $jtmp"
    echo "if want to run other job"
    echo "        $0 fname[-suffix] method job_type"
    exit
fi

method=${2:-$mtmp}
job=${3:-$jtmp}

rem=${method}_$job

cat $mol.mol > $mol.in
cat rem.$rem >> $mol.in

qchem $mol.in $mol.out &

