#!/bin/bash


#cmd=$1
dir=$1
fname=$2
#input=$1

if [ $# -lt 1 ]; then
    echo "Argument required: qchem input type"
    exit
fi    

#for f in $(pwd)/*; do
for f in $dir/*.n; do
    #fname=$(basename -- "$f")
    #ext="${f##*.}"
    #fhead="${f%.*}"
    echo $f
    #qcin_mol_rem.pl m=$fhead.out r=rem.b3lyp_vdz_opt i=$fhead.in # i is not working
    qcrun.sh 
    #if [ $ext == $input ]; then
        # for qcout_mol_in.pl m=$f r=$f
     #   echo "$cmd m=$fname r=$fname "
        #echo "xyz22mol.pl $fname"
        ### submit all .in in SGE
        #echo "qsub -N SP -v qcjob=$fhead -v nc=8 -pe numa 8 /gpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.csh"
    #fi
done    
