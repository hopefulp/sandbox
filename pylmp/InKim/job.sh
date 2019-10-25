#!/bin/bash

nrun=`qstat -anu noische | grep ' R ' | wc -l`
nqueue=`qstat -anu noische | grep ' Q ' | wc -l`

date
qstat -anu noische
echo
echo $nrun 'jobs are running.'
echo $nqueue 'jobs are in the queue.'
echo
/qcfs/isty2e/scripts/qstat_detail -owner noische
echo
