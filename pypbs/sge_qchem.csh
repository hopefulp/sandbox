#!/usr/bin/csh
#$ -cwd
#$ -V

## bash is not working
## PBS_VARIABLE is not working


if ( ! $?qcjob ) then
    echo "variable job is not defined"
    echo "Usage:: qsub -N jobname -v ver=[4.3s|5.1p]  -v job=a(.in)  sge_qchem.sh"
    exit(2)
endif

#### parallel version 5.1
/gpfs/home/joonho/sciware/qchem5.1p/bin/qchem -np $nc $qcjob.in $qcjob.out




