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

######## version 4.3
#### serial
if ( $ver == "4.3s" ) then
    /gpfs/opt/qchem_trunk_4.3/bin/qchem $qcjob.in $qcjob.out
#### parallel version 5.1
else if ( $ver == "5.1p" ) then
    /gpfs/home/joonho/sciware/qchem/bin/qchem -np $nc $qcjob.in $qcjob.out
endif    




