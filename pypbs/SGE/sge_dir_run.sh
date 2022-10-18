#!/bin/bash


#cmd=$1
input=$1

if [ $# -lt 1 ]; then
    echo "Argument required: qchem input type"
    exit
fi    

for f in $(pwd)/*; do
    fname=$(basename -- "$f")
    ext="${fname##*.}"
    fhead="${fname%.*}"
    if [ $ext == $input ]; then
        # for qcout_mol_in.pl m=$f r=$f
        echo "$cmd m=$fname r=$fname "
        #echo "xyz22mol.pl $fname"
        ### submit all .in in SGE
        #echo "qsub -N SP -v qcjob=$fhead -v nc=8 -pe numa 8 /gpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.csh"
    fi
done    
