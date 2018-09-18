#!/bin/bash

old=$1
new=$2
old2=CO2_done_sp

#Me=$3

if [ ! -e $old ]; then
    echo "There is no $old dir"
    exit
fi

if [ -e $new ]; then
    echo "There is already $new dir"
    exit
fi

mkdir $new

# change POTCAR

#cp Pot/ldapaw_${Me}OCH.pot  	$new/POTCAR
cp $old/POTCAR			$new/POTCAR

# change INCAR and DOS

#cp $old/INCAR			$new/INCAR
#cp ./incar.520.lda.cont		$new/INCAR
cp incar.520.mol		$new/INCAR
cp $old/KPOINTS     		$new/KPOINTS
cp $old/POSCAR	 		$new/POSCAR
#cp $old/CHGCAR  		$new
#cp $old/WAVECAR 		$new

~/bin/changeline.pl pbsfold-idft.csh $new nodes=1
#cp $old/KPOINTS    		$new/KPOINTS

#~/bin/changeline_tmpl.pl incar.520.lda $3 cell

#cp Pot/ggapaw_${Me}OCH.pot  	$new/POTCAR
#~/bin/changeline_tmpl.pl incar.520.CO2.gga $3
#cp t.incar			$new/INCAR


#cp ./incar.dos.co2  	$new/INCAR
#cp ./incar.dos.1co2  	$new/INCAR
#cp ./kp1.dos     	$new/KPOINTS
