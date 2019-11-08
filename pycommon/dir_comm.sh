#!/usr/bin/bash

f=$1

for d in $(ls -d p*); do
    #echo "cp $f $d"             # copy incar to all directories
    echo "cd $d"
    cd $d
    dir_clean_p2.py -w vasp -j rm
    echo "cd .."
    cd ..
    done

