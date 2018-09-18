#!/bin/bash

#old=../vasp_Relax/$1
old=$1
new=$2
nnode=$3

if [ ! -e $old ]; then
    echo "There is no $old dir"
    exit
fi

if [ -e $new ]; then
    echo "There is already $new dir"
    exit
fi

mkdir $new

cp $old/POTCAR  		$new
cp $old/CONTCAR 		$new/POSCAR
cp ./incar.bands	  	$new/INCAR
#cp ./kp4.monk     		$new/KPOINTS
#cp ./incar.dos.1co2  	$new/INCAR
#cp ./kp1.dos     	$new/KPOINTS
cp $old/KPOINTS			$new
cp $old/WAVECAR 		$new
cp $old/CHGCAR			$new
~/bin/changeline_pbs.pl pbsfold-idft.csh $new nodes=$nnode


