#!/bin/bash
#PBS -N H2O
#PBS -l nodes=1:ppn=16
#PBS -q dque
#PBS -l walltime=120:00:00
#PBS -e ./pbsvasp.err
#PBS -o ./pbsvasp.out

NPROC=`wc -l < $PBS_NODEFILE`
jobname=$PBS_JOBNAME
wdir=$jobname
log_dir=$PBS_O_WORKDIR
log_file=$log_dir/${PBS_JOBID}_$jobname

MPI="/opt/mpi/intel-12.1.6/openmpi-1.6.3/bin/mpirun -np $NPROC -hostfile $PBS_NODEFILE"
#set EXEC = "/opt/applic/vasp/bin/vasp-5.3.5-xe11-static-openmpi-1.6.3-gamma"
#EXEC="/opt/applic/vasp/bin/vasp-5.3.5-xe11-static-openmpi-1.6.3-full"
EXEC="/home/joonho/vasp.5.4.1/bin/vasp_gam_bas"
#set EXEC = "/opt/applic/vasp/bin/vasp-5.3.5-xe11-static-openmpi-1.6.3-half"
export LD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.13.367/compiler/lib/intel64:/opt/intel/composer_xe_2011_sp1.13.367/mkl/lib/intel64:/opt/mpi/intel-12.1.6/openmpi-1.6.3/lib:$LD_LIBRARY_PATH

if [ ! -d "$PBS_O_WORKDIR/$wdir"]; then
    echo "There is not $jobname directory" >> $log_file
    exit
fi

echo $jobname > $log_file
cat $PBS_NODEFILE >> $log_file
echo start >> $log_file
date >> $log_file

cd $PBS_O_WORKDIR/$wdir
if grep 'ISTART = 1' INCAR ; then
    sed -i 's/ISTART = 1/ISTART = 0/' INCAR
fi
if grep 'ICHARG = 0' INCAR ; then
    sed -i 's/ICHARG = 0/ICHARG = 2/' INCAR
fi
$MPI $EXEC > $log_dir/$jobname.log
mv $PBS_O_WORKDIR/$PBS_JOBNAME.log $PBS_O_WORKDIR/$PBS_JOBNAME.out
echo end >> $log_file
date >> $log_file

if grep 'ISIF = 3' INCAR ; then

    export i=1
    export gr=`grep E0 OSZICAR |wc -l`
    while [ $gr != 1 ]
    do
	cp POSCAR ${i}POSCAR
        cp XDATCAR ${i}XDATCAR
        cp OUTCAR ${i}OUTCAR
        #cp CHGCAR ${i}CHGCAR
        cp CONTCAR POSCAR
	if grep 'ISTART = 0' INCAR ; then
	    sed -i 's/ISTART = 0/ISTART = 1/' INCAR
	fi
	if grep 'ICHARG = 2' INCAR ; then
	    sed -i 's/ICHARG = 2/ICHARG = 0/' INCAR
	fi
	cp $PBS_O_WORKDIR/$PBS_JOBNAME.out ${i}$PBS_JOBNAME.out
	export i=`expr $i + 1`

	echo start >> $log_file
	$MPI $EXEC > $log_dir/$jobname.log
	mv $PBS_O_WORKDIR/$PBS_JOBNAME.log $PBS_O_WORKDIR/$PBS_JOBNAME.out
	echo end >> $log_file
	date >> $log_file
        export gr=`grep E0 OSZICAR |wc -l`
    done
	
fi

mv $PBS_O_WORKDIR/$PBS_JOBNAME.log $PBS_O_WORKDIR/$PBS_JOBNAME.out
    

