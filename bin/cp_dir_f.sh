#!/bin/bash

job_list=( qchem )

d_old=$1
d_new=$2

for f in `ls $d_old/*.out`
    do
	f_pre=$( basename $f | cut -d'.' -f1 )
	diff $d_old/$f_pre.out $d_new/$f_pre.out
	if [ $? -eq 0 ]; then
	    echo $f_pre.out
	else
	    echo "cp $d_old/$f_pre.out $d_new"
	fi
	
    done
