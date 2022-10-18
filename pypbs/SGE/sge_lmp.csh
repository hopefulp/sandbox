#!/usr/bin/csh
#$ -cwd
#$ -V

## bash is not working
## PBS_VARIABLE is not working


if ( ! $?np ) then
    echo "variable job is not defined"
    echo 'Usage:: qsub -N jobname(queue) -pe numa np -l mem=2.5G -v np=np -v sys=(in.)sys $SB/pypbs/sge_lmp.csh'
    echo 'Lammps:: data.system in in.lammps'
    exit(2)
endif

#### parallel version 5.4.4 compiled by Inte
set log_file = $sys.log
set queue_file = job.$sys
echo $sys > $queue_file
echo start >> $queue_file
date >> $queue_file

#cd $dir
mpirun -np $np /gpfs/home/joonho/.local/bin/lmp_mpi_user -in in.$sys
#cd ..
echo end >> $queue_file
date >> $queue_file

