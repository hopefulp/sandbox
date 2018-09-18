#!/bin/bash

dir1=cpo27_1co2_ads_
#dir1=cpo27_1co2_
#dir1=cpo27
#dir2a=CO2_sp
#dir2b=cell_sp

for Me in `cat metal.dir`
    do
	dir=$dir1${Me}_mag
	echo $dir
#	./sub_mk_dosdir.sh $dir ${dir}.extg $1 		# For DOS calculation
	./sub_mk_dosdir.sh $dir ${dir}_kb $1 $Me	# For pcharge calculation
    done


