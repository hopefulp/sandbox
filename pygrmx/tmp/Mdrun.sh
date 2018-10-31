#$ -cwd
#$ -N test
#$ -pe numa 4
#$ -V
#!/bin/bash

a=3

echo "a = $a" >> t.log


#/gpfs/opt/openmpi/bin/mpirun -np 4  mdrun -s bman_ann.tpr -c bman_ann.gro -o bman_ann.trr -e bman_ann.edr -g bman_ann.log
