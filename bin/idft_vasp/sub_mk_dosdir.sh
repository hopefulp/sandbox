#!/bin/bash

#old=../vasp_Relax/$1
old=$1
new=$2
nnode=$3
Me=$4

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
cp ./incar.bands.$Me	  	$new/INCAR
#cp $old/KPOINTS     		$new/KPOINTS
cp ./kp4.monk     		$new/KPOINTS
#cp ./incar.dos		  	$new/INCAR
#cp ./kp1.gamma	 	     	$new/KPOINTS
# for pcharge density calculation
cp $old/CHGCAR  		$new
cp $old/WAVECAR 		$new

~/bin/changeline_pbs.pl pbsfold-idft.csh $new nodes=$nnode


