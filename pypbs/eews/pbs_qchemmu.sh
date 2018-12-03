#PBS -N C202Nt-TS1-irc
#PBS -q small
#PBS -l nodes=1:ppn=32
#PBS -l walltime=120:00:00
#!/bin/tcsh

#source /qcfs/joonho/scripts/mu.pbs.qchem
source /qcfs/joonho/scripts/.tcshrc.mu.env
#source /qcfs/joonho/scripts/.tcshrc.common

set jobname = $PBS_JOBNAME
set curr_dir = $PBS_O_WORKDIR
set log_dir = $curr_dir
set nprocs = `cat $PBS_NODEFILE | wc -l`


#logging
echo $jobname > $curr_dir/$PBS_JOBID
cat $PBS_NODEFILE >> $curr_dir/$PBS_JOBID

date >> $log_dir/$PBS_JOBID
#env >> $log_dir/$PBS_JOBID
echo $curr_dir >> $log_dir/$PBS_JOBID

##### run
cd $PBS_O_WORKDIR

#set rem=b3lyp_G6311ss_opt

#cat $jobname.mol rem.$rem > $jobname.in
#cat $jobname.mol rem.$rem  bas.def2-tzvp.NiFePNCOH > $jobname.in


/home/jackjack5/3party/qchem/20170414/bin/qchem -pbs -np $nprocs $jobname.in $jobname.out
#qchem -np $nprocs $jobname.in $jobname.out

date >> $log_dir/$PBS_JOBID

