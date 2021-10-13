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
jobname=$SLURM_JOB_NAME
wdir=$pdir/$jobname
partname=$SLURM_JOB_PARTITION
nodelist=$SLURM_JOB_NODELIST

logfile=${pdir}/${SLURM_JOB_ID}.${jobname}.${partname}
outfile=$pdir/$jobname.log
date >> $logfile
#echo "pdir $pdir, jobname $jobname, wdir $wdir" > $logfile
echo "HOSTNAME    JOB_NAME   Nproc" >> $logfile
echo "$partname  $jobname $SLURM_NTASKS " >> $logfile
echo "NODELIST: $nodelist"
cd $wdir

###### Modify INCAR
### change NPAR = NNode * npar_in_partition
if [ -z "$npar" ]; then
    if [ $partname == 'X1' -o $partname == 'X2' ]; then
        par=2
    else
        par=4
    fi
    npar=$(expr $SLURM_JOB_NUM_NODES \* $par )
fi
echo "npar = $npar" >> $logfile
sed -i "s/.*NPAR.*/NPAR = $npar/" INCAR
#npar=2
hmem=0
if [ $hmem -eq 1 ]; then
    if [ $partname == 'X1' ]; then
        nppn=2
    elif [ $partname == 'X2' ]; then
        nppn=2
    else
        nppn=4
    fi
    nproc=$(expr $SLURM_JOB_NUM_NODES \* $nppn )
    echo "high memory = $hmem; nproc per node $nppn for total $nproc " >> $logfile
fi
### change U-correction
ucorr=0
if [ $ucorr -eq 1 ]; then
    ldaul=$(python -m myvasp -j Ucorr -s $wdir/POSCAR -ind 0)
    ldauu=$(python -m myvasp -j Ucorr -s $wdir/POSCAR -ind 1)
    ldauj=$(python -m myvasp -j Ucorr -s $wdir/POSCAR -ind 2)
    echo ldaul, ldauu, ldauj >> $logfile
    sed -i "s/.*LDAUL.*/$ldaul/" INCAR
    sed -i "s/.*LDAUU.*/$ldauu/" INCAR
    sed -i "s/.*LDAUJ.*/$ldauj/" INCAR
fi

if [ $partname == 'X1' ]; then
    mpirun -genv I_MPI_FABRICS "shm:ofa" -np $SLURM_NTASKS  /TGM/Apps/VASP/bin/5.4.4/O2/NORMAL/vasp.5.4.4.pl2.O2.NORMAL.std.x > $outfile
else
    if [ $hmem -eq 0 ]; then
        mpirun -np $SLURM_NTASKS  /TGM/Apps/VASP/bin/5.4.4/O2/NORMAL/vasp.5.4.4.pl2.O2.NORMAL.std.x > $outfile
    else
        mpirun -np $nproc  /TGM/Apps/VASP/bin/5.4.4/O2/NORMAL/vasp.5.4.4.pl2.O2.NORMAL.std.x > $outfile
    fi
fi
date >> $logfile

mv $outfile $pdir/$jobname.out
