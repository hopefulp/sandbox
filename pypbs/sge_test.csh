#!/usr/bin/csh
#$ -cwd
#$ -N sleep
#$ -V

## bash is not working
## PBS_VARIABLE is not working

set logfile = $SGE_O_WORKDIR/log.log

echo $SGE_O_WORKDIR >> $logfile

set dirname = `basename $SGE_O_WORKDIR`
echo $dirname >> $logfile

