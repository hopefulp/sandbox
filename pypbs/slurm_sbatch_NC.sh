#!/bin/bash
##SBATCH -J Grp              # job name, get jobname from dirname
##SBATCH -p X3               # X2(16), X3(22), X4(4)  
##SBATCH -N 6               # total number of nodesmpi tasks requested
##SBATCH -n 64               # total number of mpi tasks requested

## HPC ENVIRONMENT
. /etc/profile.d/SLURM.sh

#mpiexec.hydra -np $SLURM_NTASKS  /TGM/Apps/VASP/bin/5.4.4/NORMAL/vasp_5.4.4_GRP7_NORMAL_20170903.x
#mpirun -np $SLURM_NTASKS  /home/starnmj/bin/vasp_5.4.4_GRP7_NORMAL_20170903.x

pdir=$SLURM_SUBMIT_DIR
wdir=$pdir
jobname=$SLURM_JOB_NAME
partname=$SLURM_JOB_PARTITION

logfile=${pdir}/${SLURM_JOB_ID}.${jobname}.${partname}
outfile=$pdir/$jobname.log
date >> $logfile
#echo "pdir $pdir, jobname $jobname, wdir $wdir" > $logfile
echo "HOSTNAME    JOB_NAME   Nproc" >> $logfile
echo "$partname  $jobname $SLURM_NTASKS " >> $logfile
cd $wdir

python $main $SLURM_NTASKS 

date >> $logfile
