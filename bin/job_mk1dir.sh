#!/bin/bash

old=$1
new=$2

job=cp   # cp cont pchd		

ini_d=/qcfs/joonho/Vasp_ini
#incar=incar.hjk
#incar=incar.rp.d2c
#incar=incar.rp.d2cont
#incar=incar.bands


if [ ! -d $old ]; then
    echo "There is no $old dir"
    exit
fi

if [ -d $new ]; then
    echo "There is already $new dir"
    exit
fi

echo "copy $old $new"

mkdir $new

cp $old/CONTCAR 	$new/POSCAR
#cp $old/POSCAR	 	$new
#cp $new.pos		$new/POSCAR

cp $old/POTCAR	 	$new

#cp ./incar.bands  	$new/INCAR
#cp  ./$incar		$new/INCAR
cp $old/INCAR		$new


### For KPOINTS WAVCAR CHGCAR
if [ $job == "dos" -o $job == "pchd" ]; then
    if [ $job == "pchd" ]; then
	cp $old/KPOINTS $new
    else
        cp ./kp4.monk     	$new/KPOINTS
    fi
    cp $old/CHGCAR  	$new
    cp $old/WAVECAR 	$new
elif [ $job == "cont" ]; then
    cp $old/WAVECAR     $new
    cp $old/KPOINTS     $new
else
    cp $old/KPOINTS     $new
fi

/qcfs/joonho/bin/changeline_pbs_psi.pl pbs-psi.csh $new 1
