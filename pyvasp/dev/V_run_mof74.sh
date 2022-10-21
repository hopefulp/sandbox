#!/bin/bash

##  ini :: done,  

job=$1
w_dir=$2

if [ $# -eq 0 ]; then
    echo "Usage :: $0  job w_dir  metal mag incar-cont [ini:poscar|dos:o_dir|cont:version]"
    exit
elif [ $job == "rerun" ]; then
    if [ $# -eq 2 ]; then
	echo "OK"
    else
	echo "Usage :: rerun needs two arg's"
	exit
    fi
elif [ $# -lt 6 ]; then
    echo "Usage	:: $0  job w_dir  metal mag incar-cont [ini:poscar|dos:o_dir|cont:version]"
    echo "    	:: job == rerun ini cont cpdir dos"
    echo " 	::        rerun w_dir    "
    echo " 	::        ini   w_dir metal mag ini ac"
    echo "	:: 	  cont  w_dir metal mag wav 1 "
    echo "	::	  cpdir n_dir metal mag ini o_dir "
    echo "     	:: for job==ini POSCAR from POSCAR.Fe.mol mol= \"zz, ac, py, pa etc\" "
    echo "     	:: mag : incar.AFM or FM or NM"
    echo " 	:: INCAR continous == ini wav chg"
    exit
fi

qcvin=/qcfs/joonho/binvasp
qcvi=/qcfs/joonho/VaspINI

if [ $job == "ini" -o $job == "cpdir" -o $job == "dos" ]; then
    if [ -d $w_dir ]; then
	echo "Error:: There exists $w_dir already"
	exit
    else
    	mkdir $w_dir
    fi
elif [ $job == "cont" -o $job == "rerun" -o $job == "remod" ]; then
    if [ ! -d $w_dir ]; then
	echo "Error:: There does not exist $w_dir"
 	exit
    fi
    if [ $job == "rerun" ]; then
	$qcvin/changeline_pbs_$HOST.pl pbs-$HOST.csh $w_dir
        exit
    fi
else
    echo "There is not job identification"
    exit   
fi

#mol=$3			# position from zz ac 
metal=$3		# 
mag=$4			# change magnetism, FM, AFM, NM
incar_cont=$5

##### In the case of making new directory or running at the same directory
if [ $job == "ini" ]; then
    mol=$6				# for POSCAR
elif [ $job == "cont" ]; then
    ver=$6
elif [ $job == "dos" -o $job == "cpdir" ]; then
    o_dir=$6
    if [ ! -d $o_dir ]; then
	echo "Error:: There should be $o_dir "
        exit
    fi
fi

if [ $metal == "Sc" -o $metal == "Ca" ]; then
    pot_id=${metal}sv
else
    pot_id=${metal}
fi

###### FILE copy
if [ $job == "ini" ]; then
    $qcvin/vrun_poscar.pl $qcvi/POSCAR.Fe.$mol $metal
    mv POSCAR $w_dir
    cp $qcvi/Pot/pawpbe_${pot_id}OCH.pot $w_dir/POTCAR
    cp $qcvi/kp1.gamma	$w_dir/KPOINTS
elif [ $job == "remod" ]; then
    echo "retry at $w_dir"

elif [ $job == "cont" ]; then
    mv $w_dir/POSCAR 	$w_dir/POSCAR_$ver
    cp $w_dir/CONTCAR 	$w_dir/CONTCAR_$ver
    cp $w_dir/CONTCAR 	$w_dir/POSCAR
    cp $w_dir/INCAR 	$w_dir/INCAR_$ver
    mv $w_dir.out $w_dir/$w_dir.$ver.out
elif [ $job == "cpdir" ]; then
    cp $o_dir/CONTCAR	$w_dir/POSCAR
    cp $o_dir/POTCAR	$w_dir
    cp $o_dir/KPOINTS	$w_dir
elif [ $job == "dos" ]; then
    echo "Job is dos" 
else
    echo "Job identification is wrong so exit"
    exit
fi

### INCAR
if [ $mag == "NM" ]; then
    $qcvin/vrun_incar.sh $mag $metal $incar_cont prec rpd crelax log
elif [ $mag == "AFM" -o $mag == "FM" ]; then   
    $qcvin/vrun_incar.sh $mag $metal $incar_cont.spin prec rpd sp logall
############                                                   crelax
############                                                   opt
else
    echo "There is no $mag magnetism so exit"
    exit
fi

cp INCAR $w_dir

$qcvin/changeline_pbs_$HOST.pl pbs-$HOST.csh $w_dir



