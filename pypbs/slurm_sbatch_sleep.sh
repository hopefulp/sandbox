#!/bin/bash
#SBATCH -J sleep              # job name, get jobname from dirname
#SBATCH -p X1               # X2(16), X3(22), X4(4)  
#SBATCH -N 1              # total number of nodesmpi tasks requested
#SBATCH -n 8               # total number of mpi tasks requested
#SBATCH --nodelist=n029
## HPC ENVIRONMENT
#. /etc/profile.d/SLURM.sh

sleep 720000000
