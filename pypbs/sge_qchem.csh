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

set iqc = '5.1'

#### parallel version 5.1 (keyword) && 3.2 (Jmol)
if ( $iqc == '5.1' ) then
    if ( $opt == 'save') then
        /gpfs/home/joonho/sciware/qchem5.1p/bin/qchem -save -np $nc $qcjob.in $qcjob.out $qcjob
    else
        /gpfs/home/joonho/sciware/qchem5.1p/bin/qchem -np $nc $qcjob.in $qcjob.out
    endif
else if ( $iqc == '3.2' ) then
    /gpfs/opt/qchem/bin/qchem -np $nc $qcjob.in $qcjob.out
endif    




