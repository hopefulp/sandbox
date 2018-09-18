#!/bin/bash

o_dir=$1
n_dir=$2
poscar_tag=$3

if [ $# == 2 ]; then
    echo "copy old_dir new_dir new.pos "
elif [ $# == 3 ]; then
    if [ $3 == "1" ]; then
    	echo "copy old_dir new_dir POSCAR " 
    elif [ $3 == "2" ]; then
	echo "copy old_dir new_dir CONTCAR "
    fi
else
    echo "Usage:: $0 o_dir n_dir poscar_tag[ |1|2]"
    exit
fi

qcvin=/qcfs/joonho/binvasp
qcvi=/qcfs/joonho/VaspINI

if [ -z $poscar_tag ]; then
    Poscar=$n_dir.pos
	echo "copy $Poscar $n_dir/POSCAR"
elif [ $poscar_tag -eq 1 ]; then
    Poscar=$o_dir/POSCAR
elif [ $poscar_tag -eq 2 ]; then
    Poscar=$o_dir/CONTCAR
fi

if [ -d $n_dir ]; then
    echo "There is $n_dir directory already"
    exit
else
    mkdir $n_dir
fi
	
cp $Poscar		$n_dir/POSCAR
cp $o_dir/POTCAR  	$n_dir
cp $o_dir/KPOINTS 	$n_dir
cp $o_dir/INCAR 	$n_dir


#$qcvin/changeline_pbs_$HOST.pl pbs-$HOST.csh $n_dir 

