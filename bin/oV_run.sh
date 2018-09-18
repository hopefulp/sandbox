#!/bin/bash

### Initialize  

declare -a job_list=( Usage drun drunsv cont contsv contmag new newpack insert dos dosrun pchg ase)

### if there is no argument, show usage
if [ $# -eq 0 ]; then
    fname=$(basename $0)
    echo "Usage :: $fname  Job_type  vasp_exe 1st_dir [options]"
    echo "Job_list == ${job_list[@]} "
    echo "Now drun cont contmag are working"
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

exe_vasp=$2
i_dir=$3
#### Argument test

qcvin=/qcfs/joonho/binvasp
qcvi=/qcfs/joonho/VaspINI

if [ $# -lt 3 -o $job == "Usage" ]; then
    echo "Usage :: $0  job_type vasp_exe first_dir [options] "
    echo "      :: options "
    echo "       : drun	- no options"
    echo "       : cont	- vers# "
    echo "       : contsv new newpack 	- nei_dir key "
    echo "       : newpack key is the name of directory"
    echo "       : dos, dosrun   		- old_dir "
    exit
elif [ $# -eq 3 ]; then
    if [ $job != "drun" ]; then
	echo "Error :: more options are needed for job $job"
	exit  4
    fi
elif [ $# -eq 4 ]; then
    if [ $job == "cont" -o $job == "contsv" -o $job == "drunsv" -o $job == "dos" -o $job == "dosrun" -o $job == "pchg" -o $job == "ase" ]; then
	:
    else
	echo "Error :: $job needs more arguments"
	exit 5
    fi
elif [ $# -eq 5 ]; then
    if [ $job == "contmag" ]; then
	:
    else
	echo "Error :: $job needs more arguments"
	exit 6
    fi
fi

##### Directory test
if [ ! -d $i_dir ]; then
    echo "Error:: There does not exist $i_dir"
    exit 10
fi

#if [ $job == "ini" -o $job == "cpdir" -o $job == "dos" -o $job == "dosrun" ]; then
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

dos_incar="$qcvi/Inc/incar.dos.algo.gauss"
dos_kpoints="$qcvi/Kpoints/kp228.dos"
pchg_incar="$qcvi/Inc/incar.bands"
tmp_gga="RP"
##### treat extra argument
w_dir=$i_dir
if [ $job == "ini" ]; then
    mol=$4				# for POSCAR
elif [ $job == "cont" ]; then
    ver=$4
elif [ $job == "contsv" -o $job == "contmag" -o $job == "ase" ]; then
    n_dir=$4
    metal=$5
elif [ $job == "new" -o $job == "newpack" ]; then
    newpack=$4
    key=${5:-$tmp_gga}
elif [ $job == "dos" -o $job == "dosrun" -o $job == "pchg" ]; then
    n_dir=$4
    w_dir=$n_dir
    incar=${5:-$dos_incar}
    kpoints=${6:-$dos_kpoints}
    if [ $job == "pchg" ]; then
  	incar=${5:-$pchg_incar}
    fi
    if [ -d $n_dir ]; then
    	echo "Error:: There should be not $n_dir "
    	exit
    fi
    mkdir $n_dir
fi

###### FILE copy and run

pbs=pbs-$HOST-cell.sh

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
elif [ $job == "dos" -o $job == "dosrun" -o  $job == "pchg" ]; then
    cp $i_dir/CONTCAR   $n_dir/POSCAR
    cp $i_dir/POTCAR    $n_dir
    cp $i_dir/WAVECAR   $n_dir
    cp $i_dir/CHGCAR    $n_dir
    cp $incar		$n_dir/INCAR
    cp $kpoints		$n_dir/KPOINTS 
    if [ $job == "pchg" ]; then 
        # it should copy from DOSCAR due to k-points
        # or use PROCAR in SCF for use of SCF
        cp $i_dir/KPOINTS   $n_dir	# COPY from DOS cal.
    fi
    pbs=pbs-$HOST.sh
    if [ $job == "dos" ]; then
	exit 10
    fi
elif [ $job == "ase" ]; then
    mkdir $n_dir
    cp $i_dir/CONTCAR $n_dir
fi
#echo "PBS file is $pbs"
$qcvin/changeline_pbs_${HOST}_vasp.pl $pbs $w_dir $exe_vasp
### for direct PBS
#$qcvin/changeline_pbs_${HOST}_vasp.pl pbs-${HOST}.sh $w_dir $exe_vasp

