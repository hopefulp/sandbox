#PBS -N C202-onec-ts
#PBS -l nodes=1:ppn=12
#PBS -l walltime=120:00:00
#!/bin/tcsh

set jobname = $PBS_JOBNAME
set curr_dir = $PBS_O_WORKDIR
set log_dir = $curr_dir
set nprocs = `cat $PBS_NODEFILE | wc -l`

setenv LD_LIBRARY_PATH /opt/intel/Compiler/9.1/042/lib:/opt/intel/Compiler/9.1/042/mkl/lib/em64t:/opt/mpi/intel9/mpich-1.2.7p1/lib:$LD_LIBRARY_PATH
#set path = (. /opt/mpi/intel9/mpich-1.2.7p1/bin $path)
set path = (. /opt/intel/Compiler/composer_xe_2011_sp1.13.367/mpi/openmpi-1.6.3/bin/ $path)

#logging
echo $jobname > $curr_dir/$PBS_JOBID
cat $PBS_NODEFILE >> $curr_dir/$PBS_JOBID

date >> $log_dir/$PBS_JOBID
echo $curr_dir >> $log_dir/$PBS_JOBID

##### run
cd $PBS_O_WORKDIR
#set rem=b3lyp_G631sD3_freqopt
#set rem=b3lyp_G6311ss_3copt
#set rem=b3lyp_G6311ss_2copt
#set rem=b3lyp_G6311ss_1copt
#set rem=b3lyp_G6311ss_opt
#set rem=b3lyp_G6311ss_ts
#set rem=b3lyp_G6311ss_spts
#set rem=b3lyp_lpath
#set rem=b3lyp_rpath

#set rem=M062X_opt

#cat $jobname.mol rem.$rem  > $jobname.in

qchem -np $nprocs $jobname.in $jobname.out

date >> $log_dir/$PBS_JOBID

