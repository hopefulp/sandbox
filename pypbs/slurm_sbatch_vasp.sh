#!/bin/bash
##SBATCH -J Grp              # job name, get jobname from dirname
##SBATCH -p X3               # X2(16), X3(22), X4(4)  
##SBATCH -N 6               # total number of nodesmpi tasks requested
##SBATCH -n 64               # total number of mpi tasks requested
###SBATCH --nodelist=n[053-068]

## HPC ENVIRONMENT
#. /etc/profile.d/SLURM.sh

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

###### Modify INCAR and mpinproc
### change NPAR = NNode * npar_in_partition
#hmem=1      # if 1, use mpirin nproc(nppn)
if [ $partname == 'X1' ]; then
    par=2; SLURM_CPUS_PER_NODE=8
elif [ $partname == 'X2' ]; then
    par=2; SLURM_CPUS_PER_NODE=12
elif [ $partname == 'X3' ]; then
    par=4; SLURM_CPUS_PER_NODE=20
    #if [ $hmem -eq 1 ]; then
    if [ $hmem ]; then
        par=2
    fi
elif [ $partname == 'X4' ]; then
    par=4; SLURM_CPUS_PER_NODE=24
    if [ $hmem  ]; then
        par=2
    fi
else    # if X5
    par=4; SLURM_CPUS_PER_NODE=32
fi
npar=$(expr $SLURM_JOB_NUM_NODES \* $par )

echo "XPARTITION    JOB_NAME   NNODE  Nproc" >> $logfile
echo "$partname      $jobname  $SLURM_JOB_NUM_NODES  $SLURM_NTASKS " >> $logfile
echo "ncpu per node: $SLURM_CPUS_PER_NODE in $partname" >> $logfile
echo "NODELIST: $nodelist"
cd $wdir

echo "npar = $npar" >> $logfile
sed -i "s/.*NPAR.*/NPAR = $npar/" INCAR
if [ $hmem ]; then
    mpiproc=$(expr $SLURM_JOB_NUM_NODES \* $SLURM_CPUS_PER_NODE / 2 )
    echo "high memory = $hmem; mpiproc $mpiproc per node is half of $SLURM_CPUS_PER_NODE" >> $logfile
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
### if not version, v6; passed via --export=v=5
if [[ -v v ]]; then
    version="5.4.4.pl2"
    vasp_bin="20241206_OLD_BIN/5.4.4.pl2/"
else
    version="6.4.2"
    vasp_bin="6.4.2/NORMAL/"
fi

vasp_dir="/TGM/Apps/VASP/VASP_BIN/${vasp_bin}"

### Default to "std" if exe is unset: exe = [gam|ncl|std] for 5.4.4, 
exe=${exe:-std}

if [[ $version == "5.4.4.pl2" ]]; then
    EXEC="vasp.${version}.${exe}.x"
else
    EXEC="vasp.${version}.normal.${exe}.x"
fi

if [ $hmem ]; then
    echo "mpirun -np $mpiproc  ${vasp_dir}${EXEC}" >> $logfile
    mpirun -np $mpiproc  ${vasp_dir}${EXEC} > $outfile
else
    echo "mpirun -np $SLURM_NTASKS  ${vasp_dir}${EXEC}" >> $logfile
    mpirun -np $SLURM_NTASKS  ${vasp_dir}${EXEC} > $outfile
fi
date >> $logfile

mv $outfile $pdir/$jobname.out

