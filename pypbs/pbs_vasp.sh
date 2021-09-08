#!/bin/sh
#PBS -V
#PBS -A vasp
##PBS -N pe500
#PBS -q normal
#PBS -l select=1:ncpus=64:mpiprocs=64:ompthreads=1
#PBS -l walltime=04:00:00

if [ -z $PBS_JOBNAME ]; then
    echo "Usage:: qsub -N dirname $SB/pypbs/pbs_vasp.py"
    exit 1
fi

log_dir=$PBS_O_WORKDIR
jobname=$PBS_JOBNAME
wdir=$jobname
log_file=$log_dir/${PBS_JOBID}_$jobname

echo $jobname > $log_file
NPROC=`wc -l < $PBS_NODEFILE`
echo "NPROC = $NPROC" >> $log_file
echo start >> $log_file
date >> $log_file

EXEC="/home01/x1813a01/vasp.5.4.4/bin/vasp"

cd $log_dir/$wdir
#mpirun -np 64 $EXEC > $log_dir/$jobname.log
mpirun  $EXEC > $log_dir/$jobname.log
mv $log_dir/$jobname.log $log_dir/$jobname.out 
echo end >> $log_file
date >> $log_file


