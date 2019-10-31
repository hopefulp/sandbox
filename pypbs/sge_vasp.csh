#!/usr/bin/csh
#$ -cwd
#$ -V

## bash is not working
## PBS_VARIABLE is not working


if ( ! $?np ) then
    echo "variable job is not defined"
    echo 'Usage:: qsub -N jobname -v np=np $SB/pypbs/sge_vasp.csh'
    exit(2)
endif

#### parallel version 5.4.4 compiled by Intel

set log_file = "job.log"

echo $jobname > $log_file
echo start >> $log_file
date >> $log_file
mpirun -np $np /gpfs/home/joonho/vasp.5.4.4/bin/vasp 
echo end >> $log_file
data >> $log_file


