#!/bin/bash

### sandbox directories on each server, N.B. in sub shell, ENV is not inherited
#hostname=$(hostname)
#if [ $hostname -eq "chi" ]; then
#    SB=/home/joonho/sandbox_gl
#elif [ $hostname -eq "login" ]; then
#    SB=/gpfs/home/joonho/sandboxg
#else
#    SB=/home01/x1813a01/sandbox
#fi

server_dir="x1813a01@nurion.ksc.re.kr:/scratch/x1813a01"

if [ $# -eq 0 ]; then
    echo "Usage:: $0 [python_dir] job_index \$2[next option]"
    echo "job index: [0]scp_kisti "
    exit
fi

dir_name="pycommon"
if [[ $1 =~ [a-z] ]]; then
    dir=${1:-$dir_name}
    shift
else
    dir=$dir_name
fi

SBB=${SB}${dir}
echo $SBB

### job 0        1      2     3    4  for CASE
job=( "scp_kisti" "run_python" )


j=${job[$1]}
f=$2

case $j in
    ### SCP file to KISTI
    "scp_kisti")
        echo "scp $f $server_dir"
        ;;
    "run_python")
        echo "python $SBB/$2 ${@:3:3}"
        ;;
    ### default::
    *)
        echo "Not so many jobs now"
        ;;
esac
