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

#EXEC="/home01/g129a89/vasp.5.4.4/bin/vasp"
EXEC="/home01/x2232a02/bin/vasp_5.4.4_GRP7_NORMAL_20170903.x"
cd $log_dir/$wdir

mpirun $EXEC > $log_dir/$jobname.log

mv $log_dir/$jobname.log $log_dir/$jobname.out 

echo end >> $log_file
date >> $log_file


