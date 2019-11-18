#!/bin/bash

server_dir="x1813a01@nurion.ksc.re.kr:/scratch/x1813a01"
SBB=/gpfs/home/joonho/sandboxg/pycommon/
### job 0        1      2     3    4  for CASE
job=( "scp_kisti" "run_python" )

if [ $# -eq 0 ]; then
    echo "Usage:: $0 job_index \$2[next option]"
    echo "job index: [0]scp_kisti "
    exit
fi

j=${job[$1]}
f=$2

case $j in
    ### SCP file to KISTI
    "scp_kisti")
        echo "scp $f $server_dir"
        ;;
    "run_python")
        echo "python $SBB/$2 $3 $4"
        ;;
    ### default::
    *)
        echo "Not so many jobs now"
        ;;
esac
