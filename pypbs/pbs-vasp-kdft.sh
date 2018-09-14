#!/bin/bash
#PBS -N 
#PBS -l nodes=1:ppn=12
#PBS -q dque
#PBS -l walltime=480:00:00
#PBS -e ./pbsvasp.err
#PBS -o ./pbsvasp.out

NPROC=`wc -l < $PBS_NODEFILE`
jobname=$PBS_JOBNAME
wdir=$jobname
log_dir=$PBS_O_WORKDIR
log_file=$log_dir/${PBS_JOBID}_$jobname

MPI="/opt/intel/impi/5.0.3.049/intel64/bin/mpirun -n $NPROC -machinefile $PBS_NODEFILE"
EXEC="/opt/applic/vasp/bin/vasp-std"
export LD_LIBRARY_PATH=/opt/intel/composer_xe_2013_sp1.2.144/compiler/lib/intel64:/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64:/opt/intel/impi/5.0.3.049/intel64/lib:$LD_LIBRARY_PATH

if [ ! -d "$PBS_O_WORKDIR/$wdir"]; then
    echo "There is not $jobname directory" >> $log_file
    exit
fi

echo $jobname > $log_file
cat $PBS_NODEFILE >> $log_file
echo start >> $log_file
date >> $log_file

cd $PBS_O_WORKDIR/$wdir
$MPI $EXEC > $log_dir/$jobname.log
mv $PBS_O_WORKDIR/$PBS_JOBNAME.log $PBS_O_WORKDIR/$PBS_JOBNAME.out
echo end >> $log_file
date >> $log_file

