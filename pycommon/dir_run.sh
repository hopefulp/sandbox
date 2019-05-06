#!/bin/bash

cmd=$1


for f in $(pwd)/*; do
    fname=$(basename -- "$f")
    ext="${fname##*.}"
    fhead="${fname%.*}"
    if [ $ext == "in" ]; then
        #echo "$cmd i=$fname"
        echo "qsub -N SP -v qcjob=$fhead -v nc=8 -pe numa 8 /gpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.csh"
    fi
done    
