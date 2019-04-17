#!/bin/bash

sandbox="$HOME/sandbox_gl"

PYTHON="$HOME/anaconda3/bin/python"
EXE="$sandbox/py_ai/amp_ene.py"

if [ $# -lt 4 ]; then
    echo "Error:: It requires input argument"
    echo "Usage:: amp_loop.bash scan val water128.extxyz 6"
    echo "1 : scan or not"
    echo "2 : [tr|val]"
    echo "  3         4 "
    echo " file-name  num-of-core"
    exit 1
else
    scan=$1
    py_job=$2
    fname=$3
    nc=$4
fi    

if [  $1 == 'noscan'  ]; then
    echo "pass"
    #$PYTHON $EXE $fname $py_job -n 5 -hl 4 4 4 -el 0.001 -g
else
    for i in $(seq 2 2 10); do
        echo $i >> a.log
        ### repeat 5 times iteration
        for j in  1 2 3 4 5 ; do
            echo $i $j >> a.log
            ### validation check 1~4 parts among 5 parts
            if [ $2 == 'val' ]; then
                for k in  0 1 2 3 ; do
                    echo $i $j $k >> a.log
                    echo "$PYTHON $EXE $fname $py_job -n 5 -hl $i $i $i -el 0.001 -i $k -g -nc $nc"
                    $PYTHON $EXE $fname $py_job -n 5 -hl $i $i $i -el 0.001 -i $k -g -nc $nc
                    done
            else
                #$PYTHON $EXE $fname $py_job -n 5 -hl $i $i $i -el 0.001 -g
                echo "pass"
            fi
            done
        done
fi

