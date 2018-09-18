#!/bin/bash

old=$1
new=$2

if [ ! -e $old ]; then
    echo "There is no $old dir"
    exit
fi

if [ -e $new ]; then
    echo "There is already $new dir"
    exit
fi

mkdir $new

cp $old/POSCAR	 	$new/POSCAR
cp $old/POTCAR	 	$new/POTCAR
#cp ./incar.520.fine.d2  $new/INCAR
#cp ./incar.bands  	$new/INCAR
cp $old/INCAR		$new
cp $old/KPOINTS		$new
#cp ./kp2.monk     	$new/KPOINTS
#cp $old/CHGCAR  	$new
#cp $old/WAVECAR 	$new
