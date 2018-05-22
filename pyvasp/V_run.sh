#!/bin/bash

### Initialize  

declare -a job_list=( Usage drun drunsv cont contdir contsv contmag new newpack insert dos pchg ase vib )

pbs_file=$4

### if there is no argument, show usage
if [ $# -eq 0 ]; then
    fname=$(basename $0)
    #echo "Usage :: $fname  Job_type  1st_dir 2nd_dir [options]"
    echo "Usage :: $fname  Job_type  1st_dir 2nd_dir pbs_file[for run]"
    echo "Job_list == ${job_list[@]} "
    echo "Now dos vib ase [pbs_file running option] working"
    exit 1
fi

job=$1

### if job in argument is not specified in the Job_list, show Job_list

containsElement () {
    local e
    
    for e in "${@:2}"
	do
	    [[ "$e" == "$1" ]] && return 0
	done
    return 1
}

containsElement "$job"  "${job_list[@]}"
if [ $? -eq 1 ]; then
    echo "There is no job \"$job\" in the Job_list"
    echo "Job_list == ${job_list[@]}"
    exit 2
fi

i_dir=$2
#### Argument test

qcvin=/qcfs/joonho/binvasp
qcvi=/qcfs/joonho/VaspINI

if [ $# -lt 2 -o $job == "Usage" ]; then
    echo "Usage :: $0  job_type vasp_exe first_dir [options] "
    echo "      :: options "
    echo "       : drun	- no options"
    echo "       : cont	- vers# "
    echo "       : contsv new newpack 	- nei_dir key "
    echo "       : newpack key is the name of directory"
    echo "       : dos vib    		- old_dir "
    exit
elif [ $# -eq 2 ]; then
    if [ $job != "drun" ]; then
	echo "Error :: more options are needed for job $job"
	exit  4
    fi
fi

##### Directory test
if [ ! -d $i_dir ]; then
    echo "Error:: There does not exist $i_dir"
    exit 10
fi

#if [ $job == "ini" -o $job == "cpdir" -o $job == "dos"  ]; then
#    if [ -d $i_dir ]; then
#	echo "Error:: There exists $i_dir already"
#	exit
#    else
#    	mkdir $i_dir
#    fi
#fi

#mol=$3			# position from zz ac 
#metal=$3		# 
#mag=$4			# change magnetism, FM, AFM, NM
#incar_cont=$5

#dos_incar="$qcvi/Inc/incar.dos.algo.gauss"
dos_kpoints="$qcvi/Kpoints/kp222.dos"
kpoints_g="$qcvi/kp1.gamma"
pchg_incar="$qcvi/Inc/incar.bands"
tmp_gga="PE"

##### treat extra argument
w_dir=$i_dir
if [ $job == "ini" ]; then
    mol=$3				# for POSCAR
elif [ $job == "cont" ]; then
    ver=$3
elif [ $job == "contsv" -o $job == "contmag" ]; then
    n_dir=$3
#    metal=$4
elif [ $job == "new" -o $job == "newpack" ]; then
    newpack=$3
    key=${5:-$tmp_gga}
elif [ $job == "dos" -o $job == "vib" -o $job == "pchg" -o $job == "ase" -o $job == "contdir" ]; then
    n_dir=$3
    w_dir=$n_dir

    if [ $job == "pchg" ]; then
  	incar=$pchg_incar
  	#incar=${4:-$pchg_incar}
    fi
    if [ -d $w_dir ]; then
    	echo "Error:: There should be not $w_dir "
    	exit
    fi
    mkdir $w_dir
fi

###### FILE copy and run


if [ $job == "drun" ]; then
	:
elif [ $job == "ini" ]; then
    $qcvin/vrun_poscar.pl $qcvi/POSCAR.Fe.$mol $metal
    mv POSCAR $i_dir
    cp $qcvi/Pot/pawpbe_${pot_id}OCH.pot $i_dir/POTCAR
    cp $qcvi/kp1.gamma	$i_dir/KPOINTS
elif [ $job == "remod" ]; then
    echo "retry at $i_dir"

elif [ $job == "cont" ]; then
    mv $i_dir/POSCAR 	$i_dir/POSCAR_$ver
    cp $i_dir/CONTCAR 	$i_dir/CONTCAR_$ver
    cp $i_dir/CONTCAR 	$i_dir/POSCAR
    cp $i_dir/INCAR 	$i_dir/INCAR_$ver
    mv $i_dir.out $i_dir/$i_dir.$ver.out
    cd $i_dir
	$qcvin/changeline_incar.pl $job
        cd ..
    #echo "not running: modify INCAR for continous job and run V_run.sh"
elif [ $job == "contsv" -o $job == "contmag" ]; then
    mv $i_dir.out $i_dir
    cd $i_dir
	if [ -d $newdir ]; then
	    echo "there is sub-directory already, so stop"
	    exit 1
	else
    	    mkdir $newdir
	fi
      	mv * $newdir
	cd $newdir
	    cp POSCAR INCAR POTCAR KPOINTS CONTCAR ..
	    cd ..
	cp CONTCAR POSCAR
	if [ $job == "contmag" ]; then
	    cp WAVECAR CHGCAR INCAR KPOINTS $newdir
	    $qcvin/insert_mag.pl INCAR $metal 6 FMAll
	    cp t.incar INCAR
	fi
	cd ..
elif [ $job == "new" -o $job == "newpack" ]; then
    if [ ! -d $newpack ]; then
	echo "make dir $newpack and copy pbs-$HOST.sh"
    	exit
    fi
    n_dir=${i_dir}-$key
    V_mkdir.sh  $i_dir $newpack/$n_dir
    cd $newpack/$n_dir
	$qcvin/changeline_incar.pl "GGA=PE"
	    cd ..
	$qcvin/changeline_pbs_${HOST}_vasp.pl pbs-$HOST.sh $n_dir $exe_vasp
	cd ..
    
elif [ $job == "cpdir" ]; then
    cp $o_dir/CONTCAR	$i_dir/POSCAR
    cp $o_dir/POTCAR	$i_dir
    cp $o_dir/KPOINTS	$i_dir
elif [ $job == "dos" -o $job == "pchg" -o $job == "vib" -o $job == "contdir" ]; then
    cp $i_dir/CONTCAR   $w_dir/POSCAR
    cp $i_dir/POTCAR    $w_dir

    if [ $job == "pchg" ]; then 
        # it should copy from DOSCAR due to k-points
        # or use PROCAR in SCF for use of SCF
	$qcvin/incar_pchg.pl 	$i_dir/INCAR
	cp INCAR.pchg		$w_dir/INCAR
        cp $i_dir/KPOINTS   	$w_dir	
        cp $i_dir/WAVECAR   	$w_dir
	cp $i_dir/CHGCAR	$w_dir
    elif [ $job == "dos" ]; then	# copy CHGCAR for DOS cal
	$qcvin/incar_dos.pl 	$i_dir/INCAR
	cp INCAR.dos 		$w_dir/INCAR
	cp $dos_kpoints 	$w_dir/KPOINTS
        cp $i_dir/CHGCAR    	$w_dir
    elif [ $job == "vib" ]; then
	$qcvin/incar_vib.pl 	$i_dir/INCAR
	cp INCAR.vib		$w_dir/INCAR
	cp $qcvi/Kpoints/kp1.gamma $w_dir/KPOINTS
	cp $i_dir/CHGCAR	$w_dir
#	cp $i_dir/WAVECAR	$w_dir
    elif [ $job == "contdir" ]; then
	cp $kpoints_g       	$w_dir
	cp $i_dir/INCAR		$w_dir
	cp $i_dir/WAVECAR       $w_dir
	cp $i_dir/CHGCAR        $w_dir
	cd $w_dir
	    sed -i 's/ISTART = 0/ISTART = 1/' INCAR
	    sed -i 's/ICHARG = 2/ICHARG = 0/' INCAR
	cd ..
    fi
elif [ $job == "ase" ]; then
    cp $i_dir/CONTCAR $w_dir
    mkdir $w_dir/ini
    cp $i_dir/CONTCAR $w_dir/ini
fi
#echo "PBS file is $pbs"
if [ -z $pbs_file ]; then
    echo "not running now: input pbs-file"
else
    $qcvin/changeline_pbs.pl $pbs_file $w_dir
fi
### for direct PBS
#$qcvin/changeline_pbs_${HOST}_vasp.pl pbs-${HOST}.sh $w_dir $exe_vasp

