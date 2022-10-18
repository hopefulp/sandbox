#!/usr/bin/csh

set dir = $1
set np = $2

set log_file = $dir.log
set queue_file = log.$dir
echo $dir > $queue_file
echo start >> $queue_file
date >> $queue_file

cd $dir
mpirun -np $np /gpfs/home/joonho/vasp.5.4.4/bin/vasp  > ../$log_file
cd ..
mv $log_file $dir.out
echo end >> $queue_file
date >> $queue_file

