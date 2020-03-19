#!/usr/bin/csh
#$ -cwd
#$ -V

## bash is not working
## PBS_VARIABLE is not working


if ( ! $?qcjob ) then
    echo "variable job is not defined"
    echo "Usage:: qsub -N jobname -pe numa np -l mem=2G -v np=np -v qcjob=a(.in) $SB/sge_qchem.csh"
    echo "Usage:: qsub -N jobname -pe numa np -l mem=2G -v np=np -v qcjob=a(.in) [-v iqc=5.1p -v save=ok]  $SB/sge_qchem.csh"
    exit(2)
endif

# depricate below to get from CLI
set iqc = '5.1p'    
set save = no

#if ( $qcjob =~ 'in' ) then
#    set qcjob = `cut -d. -f1 <<< $qcjob`
#endif

#### parallel version 5.1 (keyword) && 3.2 (Jmol)
if ( $iqc == '5.1p' ) then
    set QCHEM = /gpfs/home/joonho/sciwares/qchem5.1p/bin/qchem
else if ( $iqc == '5.1pt' ) then    
    set QCHEM = /gpfs/home/joonho/sciwares/qchem5.1pt/bin/qchem
else if ( $iqc == '3.2p' ) then
    set QCHEM = /gpfs/opt/qchem/bin/qchem
endif

if ( $save == 'ok') then
    #mpirun -np $np $QCHEM -save $qcjob.in $qcjob.out $qcjob
    $QCHEM -np $np -save $qcjob.in $qcjob.out
else
    #mpirun -np $np $QCHEM  $qcjob.in $QCSCRATCH/$qcjob > $qcjob.out    # for direct run
    $QCHEM -np $np $qcjob.in $qcjob.out         # run only 1 process
endif




