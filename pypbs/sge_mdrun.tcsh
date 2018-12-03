#!/usr/bin/tcsh
#$ -cwd
#$ -pe numa 4
#$ -V

## bash is not working
## PBS_VARIABLE is not working

if ( ! $?tpr ) then
    echo "variable tpr is not defined"
    echo "use:: qsub -v tpr=mdname sge_mdrun.sh"
    exit(1)
endif    

if ( ! $?job ) then
    echo "variable job is not defined"
    exit(2)
endif

if ( $job == "md" ) then
    /gpfs/opt/openmpi/bin/mpirun -np 4  /gpfs/home/joonho/local/gmx455d_mpi/bin/mdrun -s $tpr.tpr -c ${tpr}_final.gro -o $tpr.trr -e $tpr.edr -g $tpr.log
else if ( $job == "rerun" ) then
    /gpfs/opt/openmpi/bin/mpirun -np 4  /gpfs/home/joonho/local/gmx455d_mpi/bin/mdrun -s $tpr.tpr -c ${tpr}_final.gro -o $tpr.trr -e $tpr.edr -g $tpr.log -rerun $trr
else
    echo "-v job should be one of md, rerun"
endif
