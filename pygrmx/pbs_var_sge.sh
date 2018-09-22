#!/bin/bash
#$ -cwd
#$ -N var_test
#$ -pe numa 4
#$ -V

#log_file=$PBS_JOBNAME
log_file="ne"
echo "pbs_jobname : $name" >> $log_file

