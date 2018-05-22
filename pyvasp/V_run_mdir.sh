#!/bin/bash

job=$1			# tmag [transition of magnetism]  dos
o_dir=$2
n_dir=$3
metal=$4
inc=$5

if [ $# != 5 ]; then
    echo "Usage:: $0 job[tmag|dos] old_dir new_dir metal magnetism[FM|AFM for INCAR]"
    exit
fi

qcvin=/qcfs/joonho/binvasp
qcvi=/qcfs/joonho/VaspINI

position=$o_dir/CONTCAR		 # POSCAR
potcar=$o_dir/POTCAR
kpoints=$o_dir/KPOINTS

if [ -d $n_dir ]; then
    echo "There is $n_dir directory already"
    exit
else
    mkdir $n_dir
fi
	
cp $position 		$n_dir/POSCAR
cp $potcar	  	$n_dir/POTCAR
cp $kpoints	 	$n_dir/KPOINTS

if [ $job == "tmag" ]; then
    $qcvin/vrun_incar.pl $qcvi/Inc/incar.$inc $metal
    cp INCAR $n_dir
elif [ $job == "dos" ]; then
    cp $old_dir/WAVECAR 	$new_dir
    cp $old_dir/CHGCAR  	$new_dir
    cp ./INCAR.dos 	$new_dir/INCAR
fi

$qcvin/changeline_pbs_$HOST.pl pbs-$HOST.csh $n_dir 

