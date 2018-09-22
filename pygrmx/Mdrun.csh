#!/usr/bin/tcsh
#$ -cwd
#$ -N bglc_ann
#$ -pe numa 4

#EXE="/gpfs/opt/openmpi/bin/mpirun"
#MPICH_HOME="/gpfs/opt/openmpi"
#MPI_LIB="$MPICH_HOME/lib"
#MPI_INCLUDE="$MPICH_HOME/include"
#INTEL="/gpfs/opt/intel/compilers_and_libraries/linux"
#source $INTEL/bin/compilervars.sh intel64

setenv LD_LIBRARY_PATH "/gpfs/opt/intel/lib/intel64:/gpfs/opt/intel/mkl/lib/intel64"


/gpfs/opt/openmpi/bin/mpirun -np 4  /gpfs/home/joonho/local/gmx455d_mpi/bin/mdrun -s bglc_ann.tpr -c bglc_ann.gro -o bglc_ann.trr -e bglc_ann.edr -g bglc_ann.log
