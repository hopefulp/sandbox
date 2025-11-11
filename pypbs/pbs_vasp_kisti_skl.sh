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

if [ $exe == 'gam' ]; then
    EXEC="$HOME/bin/vasp_gam"
elif [ $exe == 'xyrelax' ]; then
    EXEC="$HOME/bin/vasp_std-xy"
elif [ $exe == 'ncl' ]; then
    EXEC="$HOME/bin/vasp_ncl"
else
    EXEC="$HOME/bin/vasp_std"
fi

cd $log_dir/$wdir
### treat INCAR in wdir: remove NPAR, set NCORE
st="NNODE = $PBS_NNODES"
echo $st >> $log_file
if [[ $(grep -ic NCORE INCAR) -eq 0  &&  $(grep -ic NPAR INCAR) -ne 0 ]]  ; then
    sed -e "/NPAR/a NCORE = 20" -e "s/NPAR/\#NPAR/" -i INCAR
fi

echo "mpirun $EXEC" >> $log_file

mpirun $EXEC > $log_dir/$jobname.log

mv $log_dir/$jobname.log $log_dir/$jobname.out 

echo end >> $log_file
date >> $log_file


