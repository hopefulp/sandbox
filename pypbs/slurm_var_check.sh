#!/usr/bin/bash

npar=$1
dir=$2

SLURM_JOB_PARTITION=X1
SLURM_JOB_NUM_NODES=4
if [ -z "$npar" ]; then
    if [ $SLURM_JOB_PARTITION == 'X1' -o $SLURM_JOB_PARTITION == 'X2' ]; then
        par=2
    else
        par=4
    fi
    npar=$(expr $SLURM_JOB_NUM_NODES \* $par )
fi
echo "npar = $npar" 

natom=$(python -m myvasp -j ak -s $dir/POSCAR)
echo "natom kinds = $natom"
ucorr=$(python -m myvasp -j Ucorr -s $dir/POSCAR -ind 0)
#echo ${ucorr[@]}
echo $ucorr
#echo ${ucorr[0]} # this returns only first element split by space
