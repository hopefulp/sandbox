#PBS -l walltime=100:00:00
#PBS -l nodes=rho13:ppn=1
#PBS -N test05

#!/bin/tcsh

set NN = `cat $PBS_NODEFILE | wc -l`
sleep 7d

