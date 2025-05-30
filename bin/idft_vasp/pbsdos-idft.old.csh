#PBS -N cpo27Ti_crelax_a_d
#PBS -l nodes=4:ppn=8
#PBS -l walltime=480:00:00
#PBS -q batch
#PBS -e ./pbsvasp.err
#PBS -o ./pbsvasp.out
#### -N cpo27feRH cpo27ni6co2 cpo27znmd6co2 co2

#!/bin/tcsh

set nprocs = `wc -l < $PBS_NODEFILE`
set jobname = ${PBS_JOBNAME}os
set wdir = $jobname
set log_dir = $PBS_O_WORKDIR
set log_file = $log_dir/${PBS_JOBID}_$jobname

set DO_PARALLEL = "/home/qclab/common/ic11/mvapich2-1.6/bin/mpirun_rsh -np $nprocs -hostfile $PBS_NODEFILE"
set executable = "/home/qclab/vasp/bin/vasp-mvapich2"

#For openmpi to run vasp parallel
setenv LD_LIBRARY_PATH /opt/intel/11.1/cce/lib/intel64/:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /opt/intel/11.1/fce/lib/intel64/:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /opt/intel/11.1/mkl/lib/em64t/:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /home/qclab/common/ic11/mvapich2-1.6/lib:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /home/qclab/common/ic11/fftw3-mvapich2-1.6/lib:$LD_LIBRARY_PATH

if ( ! -d "$PBS_O_WORKDIR/$PBS_JOBNAME") then
    echo "There is not $PBS_JOBNAME directory" >> $log_file
    exit
else if ( -d "$PBS_O_WORKDIR/$wdir") then
    echo "There is $wdir directory" >> $log_file
    exit
else
    mkdir $PBS_O_WORKDIR/$wdir
endif

echo $jobname > $log_file
cat $PBS_NODEFILE >> $log_file
echo start >> $log_file
date >> $log_file

cp $PBS_O_WORKDIR/kp_dos		$PBS_O_WORKDIR/$wdir/KPOINTS
cp $PBS_O_WORKDIR/$PBS_JOBNAME/POTCAR   $PBS_O_WORKDIR/$wdir
cp $PBS_O_WORKDIR/$PBS_JOBNAME/CONTCAR  $PBS_O_WORKDIR/$wdir/POSCAR
cp $PBS_O_WORKDIR/incar.dos		$PBS_O_WORKDIR/$wdir/INCAR
cp $PBS_O_WORKDIR/$PBS_JOBNAME/CHGCAR   $PBS_O_WORKDIR/$wdir
cp $PBS_O_WORKDIR/$PBS_JOBNAME/WAVECAR  $PBS_O_WORKDIR/$wdir
#rm $PBS_O_WORKDIR/WAVECAR $PBS_O_WORKDIR/CHGCAR

#run
cd $PBS_O_WORKDIR/$wdir
$DO_PARALLEL $executable > $PBS_O_WORKDIR/$jobname.log
mv $PBS_O_WORKDIR/$jobname.log $PBS_O_WORKDIR/$jobname.out
echo end >> $log_file
date >> $log_file

rm $PBS_O_WORKDIR/$wdir/CHG
