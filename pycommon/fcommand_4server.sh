#!/bin/bash

### sandbox directories on each server, N.B. in sub shell, ENV is not inherited

server_dir="x1813a01@nurion.ksc.re.kr:/scratch/x1813a01"

if [ $# -eq 0 ]; then
    echo "Usage:: $0 [job_index]  [python_dir] \$2[next option]"
    echo "job index: [0]run_python [1]scp_kisti "
    exit
fi

### job      0           1      2     3    4  for CASE
job=( "python" "scp_kisti" )
### find job type $1 and shift
if [[ $1 =~ [0-9] ]]; then
    j=${job[$1]}
else
    if [[ $1 =~ "py" ]]; then
        j=${job[0]}
    else
        j=${job[1]}
    fi
fi
shift
#echo Job is $j in server 
### find python dir if different from pycommon
dir_name="pycommon"
if [[ $1 =~ ".py" ]]; then
    dir=$dir_name
else
    dir=${1:-$dir_name}
    shift
fi
SBB=${SB}${dir}     # $SB is exported at alias.sh
#echo $SBB

f=$1
#echo $f
case $j in
    ### run python which have different directories for python script
    "python")
        echo "python $SBB/$f ${@:2:4}"
        echo "Do you want to run? [y/n]"
        read y_or_n
        if [ $y_or_n == "y" ]; then
            python $SBB/$f ${@:2:4}
        fi
        ;;
    ### SCP file to KISTI
    "scp_kisti")
        echo "scp $f $server_dir"
        ;;
    ### default::
    *)
        echo "Default:: Not so many jobs now"
        ;;
esac
