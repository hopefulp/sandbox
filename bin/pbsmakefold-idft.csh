#PBS -N ggaCasp
#PBS -l nodes=1:ppn=8
#PBS -l walltime=480:00:00
#PBS -q batch
#PBS -e ./pbsvasp.err
#PBS -o ./pbsvasp.out
#### -N cpo27feRH cpo27ni6co2 cpo27znmd6co2 co2

#!/bin/tcsh

set Me = Ti
set nprocs = `wc -l < $PBS_NODEFILE`
set jobname = ${PBS_JOBNAME}
set wdir = $jobname
set log_dir = $PBS_O_WORKDIR
set log_file = $log_dir/${PBS_JOBID}_$jobname

#set poscar pos.${PBS_JOBNAME}                   # pos.Me or pos.Me.6co2
#set potcar pawpbe_${PBS_JOBNAME}OCH.pot
#set incar  incar.520.fine.${PBS_JOBNAME}

set DO_PARALLEL = "/home/qclab/common/ic11/mvapich2-1.6/bin/mpirun_rsh -np $nprocs -hostfile $PBS_NODEFILE"
set executable = "/home/qclab/vasp/bin/vasp-mvapich2"

#For openmpi to run vasp parallel
setenv LD_LIBRARY_PATH /opt/intel/11.1/cce/lib/intel64/:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /opt/intel/11.1/fce/lib/intel64/:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /opt/intel/11.1/mkl/lib/em64t/:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /home/qclab/common/ic11/mvapich2-1.6/lib:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /home/qclab/common/ic11/fftw3-mvapich2-1.6/lib:$LD_LIBRARY_PATH

if ( -d "$PBS_O_WORKDIR/$wdir") then
    echo "There exists $jobname directory" >> $log_file
    exit
else
    mkdir $PBS_O_WORKDIR/$wdir
endif

echo $jobname > $log_file
cat $PBS_NODEFILE >> $log_file
echo start >> $log_file
date >> $log_file

cp $PBS_O_WORKDIR/kp2.monk      		$PBS_O_WORKDIR/$wdir/KPOINTS
cp $PBS_O_WORKDIR/pos.${Me}       		$PBS_O_WORKDIR/$wdir/POSCAR
cp $PBS_O_WORKDIR/Pot/pawpbe_${Me}OCH.pot   	$PBS_O_WORKDIR/$wdir/POTCAR
cp $PBS_O_WORKDIR/incar.520.fine        	$PBS_O_WORKDIR/$wdir/INCAR

#rm $PBS_O_WORKDIR/WAVECAR $PBS_O_WORKDIR/CHGCAR

#run
cd $PBS_O_WORKDIR/$wdir
$DO_PARALLEL $executable > $PBS_O_WORKDIR/$jobname.log
mv $PBS_O_WORKDIR/$jobname.log $PBS_O_WORKDIR/$jobname.out
echo end >> $log_file
date >> $log_file

rm $PBS_O_WORKDIR/$wdir/CHG
