#!/usr/bin/tcsh
#$ -cwd
#$ -N Anneal
#$ -pe numa 4
#$ -V

## bash is not working
## PBS_VARIABLE is not working

if ( ! $?tpr ) then
    echo "variable tpr is not defined"
    echo "use:: qsub -v tpr=mdname Mdrun.sh"
    exit(1)
endif    
    
/gpfs/opt/openmpi/bin/mpirun -np 4  /gpfs/home/joonho/local/gmx455d_mpi/bin/mdrun -s $tpr.tpr -c ${tpr}_final.gro -o $tpr.trr -e $tpr.edr -g $tpr.log
