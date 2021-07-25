#!/bin/bash
#SBATCH -J TEST            # job name
#SBATCH -o stdout.txt      # output and error file name (%j expands to 
##SBATCH -p X5              # name of cluster group 
##SBATCH -N 1               # total number of nodesmpi tasks requested
##SBATCH -n 12              # total number of mpi tasks requested

export VASP_COMMAND="mpirun -genv I_MPI_FABRICS "shm:ofa" -np $SLURM_NTASKS $HOME/bin/vasp"

neb-catlearn.py

