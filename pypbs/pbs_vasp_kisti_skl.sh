#!/bin/sh
#PBS -V
#PBS -A vasp
#PBS -q normal
#PBS -l select=20:ncpus=40:mpiprocs=40:ompthreads=1
#PBS -l walltime=48:00:00

if [ -z $PBS_JOBNAME ]; then
    echo "Usage:: qsub -N dirname $SB/pypbs/pbs_vasp.sh"
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

if [ $exe ]; then
    EXEC="$HOME/bin/vasp_gam"
else
    EXEC="$HOME/bin/vasp_std"
fi

cd $log_dir/$wdir

echo "mpirun $EXEC" >> $log_file

mpirun $EXEC > $log_dir/$jobname.log

mv $log_dir/$jobname.log $log_dir/$jobname.out 

echo end >> $log_file
date >> $log_file


