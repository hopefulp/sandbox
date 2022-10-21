#!/bin/bash

machine=g1

Old_home=$PWD
New_home=/qcfs/joonho/VASP_dat
new_suff=-dos2
old_suff=-conn

job=dos		# cont contn dos
incar_suff=.sp.afm

for dir_name in $@
    do
	old_dir=$dir_name$old_suff
	New_dir=$New_home/$dir_name$new_suff
	if [ -d $New_dir ]; then
	    echo "There is $New_dir directory already"
	    exit
	else
	    mkdir $New_dir
	fi
	
	cp $old_dir/CONTCAR 	$New_dir/POSCAR
	cp $old_dir/POTCAR  	$New_dir
	cp $old_dir/KPOINTS 	$New_dir
	if [ $job == "contn" ]; then
	    cp $New_home/INCAR$incar_suff	$New_dir/INCAR
	elif [ $job == "cont" ]; then
	    cp $old_dir/WAVECAR 	$New_dir
            cp $old_dir/CHGCAR  	$New_dir
            cp $New_home/INCAR.sp.cont $New_dir/INCAR
	elif [ $job == "dos" ]; then
	    cp $old_dir/WAVECAR 	$New_dir
	    cp $old_dir/CHGCAR  	$New_dir
	    cp $New_home/INCAR.dos 	$New_dir/INCAR
	fi
	cd $New_home
	    /qcfs/joonho/bin/changeline_pbs_kdft.pl pbs-kdft.csh $dir_name$new_suff nodes=1 $machine
	cd $Old_home
    done

