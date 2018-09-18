#!/bin/bash

job=new			# new cont contn dos

old_dir=$1
new_dir=$2

incar=./INCAR.test
position=$1/POSCAR #./ac-strand-frame.vas		# POSCAR
potcar=$1/POTCAR
kpoints=$1/KPOINTS


cd $PWD
	if [ -d $new_dir ]; then
	    echo "There is $new_dir directory already"
	    exit
	else
	    mkdir $new_dir
	fi
	
	cp $position 		$new_dir/POSCAR
	cp $potcar	  	$new_dir/POTCAR
	cp $kpoints	 	$new_dir/KPOINTS
	cp $incar		$new_dir/INCAR
	if [ $job == "contn" ]; then
	    cp $New_home/INCAR$incar_suff	$new_dir/INCAR
	elif [ $job == "cont" ]; then
	    cp $old_dir/WAVECAR 	$new_dir
            cp $old_dir/CHGCAR  	$new_dir
            cp $old_dir/INCAR.sp.cont $new_dir/INCAR
	elif [ $job == "dos" ]; then
	    cp $old_dir/WAVECAR 	$new_dir
	    cp $old_dir/CHGCAR  	$new_dir
	    cp $old_dir/INCAR.dos 	$new_dir/INCAR
	fi
    	/qcfs/joonho/bin/changeline_pbs_kdft.pl pbs-idft.csh $new_dir nodes=1 
