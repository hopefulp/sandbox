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

# '5.1p' '3.2'
set iqc = '3.2p'    
set save = no

#### parallel version 5.1 (keyword) && 3.2 (Jmol)
if ( $iqc == '5.1p' ) then
    set QCHEM = /gpfs/home/joonho/sciwares/qchem5.1p/bin/qchem
    if ( $save == 'ok') then
        mpirun -np $np $QCHEM -save $qcjob.in $qcjob.out $qcjob
    else
        #mpirun -np $np $QCHEM       $qcjob.in $qcjob.out    # looks like running the same process at each process
        mpirun -np $np $QCHEM       $qcjob.in $qcjob.out
        #$QCHEM -np $np $qcjob.in $qcjob.out         # run only 1 process
    endif
else if ( $iqc == '3.2p' ) then
    set QCHEM = /gpfs/opt/qchem/bin/qchem
    $QCHEM -np $np $qcjob.in $qcjob.out             # this is working
endif    




