#!/usr/bin/csh
#$ -cwd
#$ -V

## bash is not working
## PBS_VARIABLE is not working


if ( ! $?np ) then
    echo "variable job is not defined"
    echo 'Usage:: qsub -N jobname(queue) -pe numa np -v np=np -v dir=dirname(log) $SB/pypbs/sge_vasp.csh'
    exit(2)
endif

#### parallel version 5.4.4 compiled by Inte
set log_file = $dir.log
set queue_file = job.$dir
echo $dir > $queue_file
echo start >> $queue_file
date >> $queue_file

set EXE = /gpfs/home/joonho/vasp.5.4.4/bin/vasp_std

cd $dir
mpirun -np $np $EXE  > ../$log_file
cd ..
mv $log_file $dir.out
echo end >> $queue_file
date >> $queue_file

