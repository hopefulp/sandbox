#PBS -N  hpp6-1W-B3D
#PBS -l nodes=1:ppn=16
#PBS -l walltime=120:00:00
#!/bin/tcsh

set jobname = $PBS_JOBNAME
set curr_dir = $PBS_O_WORKDIR
set log_dir = $curr_dir
set nprocs = `cat $PBS_NODEFILE | wc -l`

# for kdft
#setenv LD_LIBRARY_PATH /opt/intel/Compiler/9.1/042/lib:/opt/intel/Compiler/9.1/042/mkl/lib/em64t:/opt/mpi/intel9/mpich-1.2.7p1/lib:$LD_LIBRARY_PATH
#set path = (. /opt/mpi/intel9/mpich-1.2.7p1/bin $path)

#logging
echo $jobname > $curr_dir/$PBS_JOBID
cat $PBS_NODEFILE >> $curr_dir/$PBS_JOBID

date >> $log_dir/$PBS_JOBID
echo $curr_dir >> $log_dir/$PBS_JOBID

##### run
cd $PBS_O_WORKDIR

qchem -np $nprocs $jobname.in $jobname.out

date >> $log_dir/$PBS_JOBID

