#!/bin/bash
## usage $0 job dir version

job=$1		# re cont
dir=$2
ver=$3
qcvi=/qcfs/joonho/VaspINI
qcvin=/qcfs/joonho/binvasp

if [ $job != "re" -a $job != "cont" ]; then
    echo "Job should be \"re\" or \"cont\""
    echo "\"re\" to run at that directory; \"cont\" to cont from $qcvi"
    echo "Usage::$0 job dir version"
    exit
fi

incar=$qcvi/INCAR

if [ ! -d $dir ]; then
    echo "There is no $dir directory "
    exit
else
    if [ $job == "re" ]; then
	echo "use directory "
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
 	    cp $incar $dir/INCAR
 	    mv $dir.out $dir/$dir.$ver.out
	fi
    fi
fi

$qcvin/changeline_pbs_$HOST.pl pbs-$HOST.csh $dir 

