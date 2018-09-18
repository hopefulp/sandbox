#!/bin/bash
## usage $0 job dir 

job=$1		# re cont
dir=$2
ver=$3

#### vrun_incar.sh variable
metal=$4
mag=$5
cont=$6

qcvi=/qcfs/joonho/VaspINI
qcvin=/qcfs/joonho/binvasp

if [ $# -ne 6 ]; then
    echo "Usage::$0 job dir version mag metal cont"
    exit
fi

if [ $job != "re" -a $job != "cont" -a $job != "remod" ]; then
    echo "Job should be \"re\" or \"cont\""
    echo "\"re\" to run at that directory; \"cont\" to cont from $qcvi"
    exit
fi

#### make INCAR
if [ $# -ge 4 ]; then
    if [ $mag == "AFM" -o $mag == "FM" -o $mag == "NM" ]; then
	if [ $cont == "ini" -o $cont == "wav" -o $cont == "chg" ]; then
    	    $qcvin/vrun_incar.sh $mag $metal $cont prec rpd crelax log
    	    echo "INCAR was made"
	else
	    echo "Continuation is wrong"
	    exit 
	fi
    else
	echo " Magnetization is wrong"
	exit
    fi
else
    echo "use INCAR"
fi

if [ ! -d $dir ]; then
    echo "There is no $dir directory "
    exit
else
    if [ $job == "re" ]; then
	echo "use directory "
    elif [ $job == "remod" ]; then
	echo "modify in running at dir"
	cp INCAR $dir/INCAR
    elif [ $job == "cont" ]; then
	if [ $ver == "" ]; then
	    echo "Error:: to continue, it need version number"
	    exit
	else
     	    cd $dir
 	        mv POSCAR POSCAR_$ver
 	        cp CONTCAR CONTCAR_$ver
 	        cp CONTCAR POSCAR
 	        cp INCAR INCAR_$ver
            cd ..
 	    cp INCAR $dir/INCAR
 	    mv $dir.out $dir/$dir.$ver.out
	fi
    fi
fi

$qcvin/changeline_pbs_$HOST.pl pbs-$HOST.csh $dir 

