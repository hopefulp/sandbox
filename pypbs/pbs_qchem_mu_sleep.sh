#PBS -N qchem
#PBS -q small
#PBS -l nodes=1:ppn=32
#PBS -l walltime=120:00:00

#!/bin/tcsh

set NN = `cat $PBS_NODEFILE | wc -l`
sleep 720000000
