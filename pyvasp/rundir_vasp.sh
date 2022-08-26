#!/bin/bash

if [ $# -gt 0 ]; then
    job=$1
else
    echo "input job"
    echo "Usage:: $0 job"
    exit 1
fi

dirfile='dirs.txt'

dirs=${2:-$dirfile}

joblist=( 'pos2cif.pl' )

run=0

for d in $(cat $dirs); do
    dir=${d}dos
    #echo ${dir}
    case $job in
        'pos2cif.pl')               # to plot figure of POSCAR to MS
            s1="$job $dir/POSCAR"
            s2="cp $dir/POSCAR.cif pos_${d}.cif"
            echo $s1 
            echo $s2 
            #if [ -n "$run" ]; then
            #    echo $s1 | sh
            #    echo $s2 | sh
            #fi
            ;;
        *)
            echo "input $job is not in job list ($joblist)"
            ;;
    esac
done
