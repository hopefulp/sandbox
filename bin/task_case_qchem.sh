#!/bin/bash
### the number of arguments is one of 1, 2, 4 in the order of job_kind, file_w_list, file_frag1ix, file_frag2ix
#pre_dir=/qcfs/joonho/work/

job1="opt"
method="no"
### string for total energy except MP2
string="Total energy"

##########  0          1            2     3    4   5    6     7   8    9       10         11
job_list=( ene_opt  ene_opt_long ene_sp )

if [ $# -lt 1 ]; then
    echo "Usage::$0 job"
    exit
fi

Job=${1:-$job1}
Method=${2:-$method}
job=ene_$Job

echo $job

for file in `ls *out`  	#`ls $frag1 -d1`; `cat $Sys.txt`
    do
	### for file
	case $job in
	    "${job_list[0]}")
		grep "Final energy" $file
		;;
	    "${job_list[1]}")
		if [ $Method == "mp2" ]; then
		    grep "TOTAL ENERGY" $file
		else
		    grep "$string" $file
		fi
		;;
	    "${job_list[2]}")
		if [ $Method == "mp2" ]; then
		    grep "RIMP2         total energy" $file
		else
		    grep "$string" $file
		fi
		;;
	    "${job_list[2]}")
		if [ -f $w_dir.out ]; then
		    grep F= $w_dir.out
		elif [ -f $w_dir.log ]; then
		    echo " There is $w_dir.log "
		    grep F= $w_dir.log
		fi
		;;
	    *)
	   	echo "the keyword should be one of \" ${job_list[@]} \" "
		echo "present setting is no-job frag1ix=$frag1 $Sys-$frag2 " 
		exit
		;; 
	esac
		    if [ $? -eq 0 ]; then
			echo $file
		    fi
    done
