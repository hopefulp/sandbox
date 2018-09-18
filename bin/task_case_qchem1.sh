#!/bin/bash
### the number of arguments is one of 1, 2, 4 in the order of job_kind, file_w_list, file_frag1ix, file_frag2ix
#pre_dir=/qcfs/joonho/work/

job1="opt"
method="no"
### string for total energy except MP2
string="Total energy"

##########  0          1            2     3    4   5    6     7   8    9       10         11
job_list=( ene_opt  ene_opt_long ene_sp bsse)

if [ $# -lt 1 ]; then
    echo "Usage::$0 job"
    exit
fi

Job=${1:-$job1}
Method=${3:-$method}
file=$2

job=$Job

mp2_string="RIMP2         total"

echo $job

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
	    "${job_list[3]}")
		if [ $Method == "mp2" ]; then
		    echo "$mp2_string"
		    grep "$mp2_string" $file/*out | awk '{print $6}'  | perl -ne 'END {print "\n"} chomp($_); print "$_\t";'
		else
		    grep "$string" $file/*out | awk '{print $10}'  | perl -ne 'END {print "\n"} chomp($_); print "$_\t";'
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
