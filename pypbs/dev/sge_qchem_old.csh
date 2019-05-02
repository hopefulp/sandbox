#!/usr/bin/csh
#$ -cwd
#$ -pe numa 4
#$ -V

## bash is not working
## PBS_VARIABLE is not working


if ( ! $?job ) then
    echo "variable job is not defined"
    echo "Usage:: qsub -N jobname -v ver=[3.2s|3.2p|4.3s]  -v job=a(.in)  sge_qchem.sh"
    exit(2)
endif

### v.3.2 
###### serial
if ( $ver == "3.2s" ) then
    /gpfs/opt/qchem/bin/qchem  $job.in $job.out
###### parallel
else if ( $ver == "3.2p" ) then
    /gpfs/opt/qchem/bin/qchem -np 4 $job.in $job.out
### v.4.3
###### serial
else if ( $ver == "4.3s" ) then
    /gpfs/opt/qchem_trunk_4.3/bin/qchem $job.in $job.out
endif    




