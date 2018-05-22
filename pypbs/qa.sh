#!/bin/tcsh
if ($#argv == 0) then
 set fnw = biduri
else
 set fnw = $1
endif

qstat -nr | grep -A `qstat -n | wc -l` -B `qstat -n | wc -l` $fnw
echo ''
