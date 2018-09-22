#!/bin/bash
#$ -cwd
#$ -N bman_ann
#$ -pe numa 4
#$ -V

/gpfs/opt/openmpi/bin/mpirun -np 4  mdrun -s bman_ann.tpr -c bman_ann.gro -o bman_ann.trr -e bman_ann.edr -g bman_ann.log
