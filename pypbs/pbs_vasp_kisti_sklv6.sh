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

### if not version, v6
if [ $v ]; then
version="5.4.4.skl"
else
version="6.3.1"
fi
vasp_dir="$HOME/bin"

### exe = [gam|ncl|std] for 5.4.4
if [ $exe ]; then
    EXEC="${vasp_dir}/vasp.${version}.{exe}.x"
elif [ $crelax ]; then
    EXEC="${vasp_dir}/vasp.${verison}.std-xy"
else
    EXEC="${vasp_dir}/vasp.${version}.std.x"
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


