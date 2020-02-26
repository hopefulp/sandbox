#!/usr/bin/csh
#$ -cwd
#$ -V

## bash is not working
## PBS_VARIABLE is not working


if ( ! $?qcjob ) then
    echo "variable job is not defined"
    echo "Usage:: qsub -N jobname -pe numa np -v np=np -v np=np -v qcjob=a(.in) $SB/sge_qchem.csh"
    echo "Usage:: qsub -N jobname -pe numa np  -v np=np -v qcjob=a(.in) [-v iqc=5.1p -v save=ok]  $SB/sge_qchem.csh"
    exit(2)
endif

set iqc = '5.1p'
set save = no
#### parallel version 5.1 (keyword) && 3.2 (Jmol)
if ( $iqc == '5.1p' ) then
    if ( $save == 'ok') then
        /gpfs/home/joonho/sciwares/qchem5.1p/bin/qchem -save -np $np $qcjob.in $qcjob.out $qcjob
    else
        /gpfs/home/joonho/sciwares/qchem5.1p/bin/qchem -np $np $qcjob.in $qcjob.out
    endif
else if ( $iqc == '3.2' ) then
    /gpfs/opt/qchem/bin/qchem -np $np $qcjob.in $qcjob.out
endif    




