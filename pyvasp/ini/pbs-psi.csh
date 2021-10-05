#PBS -N OCBOO-dopedGrp-is2
#PBS -l nodes=1:ppn=16
#PBS -q dque
#PBS -l walltime=480:00:00
#PBS -e ./pbsvasp.err
#PBS -o ./pbsvasp.out

#!/bin/tcsh
set NPROC = `wc -l < $PBS_NODEFILE`
set jobname = $PBS_JOBNAME
set wdir = $jobname
set log_dir = $PBS_O_WORKDIR
set log_file = $log_dir/${PBS_JOBID}_$jobname

set MPI = "/opt/mpi/intel-12.1.6/openmpi-1.6.3/bin/mpirun -np $NPROC -hostfile $PBS_NODEFILE"
#set EXEC = "/opt/applic/vasp/bin/vasp-5.3.3-xe11-static-openmpi-1.6.3-gamma"
set EXEC = "/opt/applic/vasp/bin/vasp-5.3.3-xe11-static-openmpi-1.6.3-full"
#set EXEC = "/opt/applic/vasp/bin/vasp-5.3.3-xe11-static-openmpi-1.6.3-half"
setenv LD_LIBRARY_PATH /opt/intel/composer_xe_2011_sp1.13.367/compiler/lib/intel64:/opt/intel/composer_xe_2011_sp1.13.367/mkl/lib/intel64:/opt/mpi/intel-12.1.6/openmpi-1.6.3/lib:$LD_LIBRARY_PATH

if ( ! -d "$PBS_O_WORKDIR/$wdir") then
    echo "There is not $jobname directory" >> $log_file
    exit
endif

echo $jobname > $log_file
cat $PBS_NODEFILE >> $log_file
echo start >> $log_file
date >> $log_file

cd $PBS_O_WORKDIR/$wdir
$MPI $EXEC > $log_dir/$jobname.log
mv $PBS_O_WORKDIR/$PBS_JOBNAME.log $PBS_O_WORKDIR/$PBS_JOBNAME.out
echo end >> $log_file
date >> $log_file

