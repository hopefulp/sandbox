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

if [ $exe ]; then
    EXEC="$HOME/bin/vasp_gam"
else
    EXEC="$HOME/bin/vasp_std"
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
cd $log_dir/$wdir

mpirun $EXEC > $log_dir/$jobname.log

mv $log_dir/$jobname.log $log_dir/$jobname.out 

echo end >> $log_file
date >> $log_file

sfiles=( POSCAR XDATCAR OUTCAR INCAR )
fail="error"
i=1
grep -q "$fail" $log_dir/$jobname.out
while [ $? -eq 0 ]; do
    for f in ${sfiles[@]}; do
        cp $f ${i}${f}
        done
    cp CONTCAR POSCAR
    if grep 'ISTART = 0' INCAR ; then
        sed -i 's/ISTART = 0/ISTART = 1/' INCAR
    fi
    if grep 'ICHARG = 2' INCAR ; then
        sed -i 's/ICHARG = 2/ICHARG = 0/' INCAR
    fi
    mv $log_dir/${jobname}.out ${i}${jobname}.out
    export i=`expr $i + 1`
    echo start >> $log_file
    date >> $log_file
    mpirun $EXEC > $log_dir/$jobname.log
    mv $log_dir/$jobname.log $log_dir/$jobname.out 
    echo end >> $log_file
    date >> $log_file
    grep -q "$fail" $log_dir/$jobname.out
    done
