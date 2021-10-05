#!/bin/bash
#PBS -N MD-Q001-carba
#PBS -l nodes=2:ppn=16
#PBS -q dque
#PBS -l walltime=240:00:00
#PBS -e ./pbsvasp.err
#PBS -o ./pbsvasp.out

NPROC=`wc -l < $PBS_NODEFILE`
jobname=$PBS_JOBNAME
wdir=$jobname
log_dir=$PBS_O_WORKDIR
log_file=$log_dir/${PBS_JOBID}_$jobname

MPI="/opt/mpi/intel-12.1.6/openmpi-1.6.3/bin/mpirun -np $NPROC -hostfile $PBS_NODEFILE"
EXEC="/opt/applic/vasp/bin/vasp-5.3.3-xe11-static-openmpi-1.6.3-gamma"
#EXEC="/opt/applic/vasp/bin/vasp-5.3.3-xe11-static-openmpi-1.6.3-half"
#EXEC="/opt/applic/vasp/bin/vasp-5.3.3-xe11-static-openmpi-1.6.3-full"
#EXEC= # for V_run.sh
export LD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.13.367/compiler/lib/intel64:/opt/intel/composer_xe_2011_sp1.13.367/mkl/lib/intel64:/opt/mpi/intel-12.1.6/openmpi-1.6.3/lib:$LD_LIBRARY_PATH

if [ ! -d "$PBS_O_WORKDIR/$wdir"]; then
    echo "There is not $jobname directory" >> $log_file
    exit
fi

echo $jobname > $log_file
cat $PBS_NODEFILE >> $log_file
echo start >> $log_file
date >> $log_file


echo $EXEC >> $log_file
if [ -z "$EXEC" ]; then
    echo "EXEC is not defined" >> $log_file 	# blank|tap is working for -zero
    exit
fi

cd $PBS_O_WORKDIR/$wdir
if grep 'ISTART = 1' INCAR ; then
    sed -i 's/ISTART = 1/ISTART = 0/' INCAR
fi
if grep 'ICHARG = 0' INCAR ; then
    sed -i 's/ICHARG = 0/ICHARG = 2/' INCAR
fi
#if grep 'POTIM = 0.3' INCAR ; then
#    sed -i 's/POTIM = 0.3/POTIM = 0.1/' INCAR
#fi

$MPI $EXEC > $log_dir/$jobname.log
mv $PBS_O_WORKDIR/$PBS_JOBNAME.log $PBS_O_WORKDIR/$PBS_JOBNAME.out
echo end >> $log_file
date >> $log_file

### for no-cell relaxation, this part is skipped
if grep 'ISIF = 3' INCAR ; then

    export i=1
    export gr=`grep E0 OSZICAR |wc -l`
    while [ $gr != 1 ]
    do
	cp POSCAR ${i}POSCAR
        cp XDATCAR ${i}XDATCAR
        cp OUTCAR ${i}OUTCAR
	cp INCAR  ${i}INCAR
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
    

