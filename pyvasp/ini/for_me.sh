#!/bin/bash

qcvin=/qcfs/joonho/binvasp

for me in `more metal.dat`
    do
	echo $me
	#$qcvin/V_run_cont.sh ${me}-zz-fm ${me}-zz-fm1
	$qcvin/V_run_ini.sh ea $me nonmag 1 nm 
	#$qcvin/V_run_ini.sh ea $me FM 1 fm 
	#$qcvin/V_run_ini.sh py $me AFM 1 fm 
    done
