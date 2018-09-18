#!/bin/bash

#old=../vasp_relax/$1
old=$1
new=$2
Me=$3
nnodes=$4
sys=$5 				# sys = cell or CO2 or 1CO2
func=$6				# lda or gga
INCAR=incar.520.fine.d2.cont		# for PBE + D2
#INCAR=incar.520.$func
#INCAR=incar.520.fine
if [ ! -e $old ]; then
    echo "There is no $old dir"
    exit
fi

if [ -e $new ]; then
    echo "There is already $new dir"
    exit
fi

mkdir $new

#cp Pot/${func}paw_${Me}OCH.pot  	$new/POTCAR	# for different Method, different POTCAR is needed
cp $old/POTCAR 				$new

cp $old/CONTCAR 			$new/POSCAR

~/bin/changeline_incar.pl $INCAR  $Me $sys
cp t.incar				$new/INCAR
#cp ./$INCAR				$new/INCAR

cp $old/KPOINTS				$new
#cp ./kp4.monk     			$new/KPOINTS
#cp ./kp1.gamma    			$new/KPOINTS

#cp $old/CHGCAR  			$new
#cp $old/WAVECAR 			$new

~/bin/changeline_pbs.pl pbsfold-idft.d2.csh $new nodes=$nnodes
