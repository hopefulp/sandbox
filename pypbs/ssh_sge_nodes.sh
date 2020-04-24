#!/bin/bash

nodes=$(qstat -f | sed -e '/---/d' -e '/adus/d' | awk '/@/ { print $1 }' | awk -F@ '{print $2}'  )
#echo $nodes

jobs=( ls rm mkdir vasp qchem ln )
job=${1:-"vasp"}
#echo $job
if [ $job == "ls" -o $job == "mkdir" -o $job == "rm" ]; then
    dir=$2
    if [ -z $dir ]; then
        echo "input directory as $2"
        exit 1
    fi
fi

if [ $# -eq 0 ]; then
    echo "This runs $0 $job"
fi

njob=0

for node in $nodes; do
    echo $node
    case $job in
        "vasp")
            ssh $node ps aux | grep $job
            ;;
        "ls")
            ssh $node ls $dir
            if [ $? -eq 0 ]; then
                njob=$(expr $njob + 1)
            fi
            ;;
        "rm")
            for f in $(ssh $node ls $dir); do
                ssh $node rm -r $dir/$f
                done
            ;;
        "mkdir")
            ssh $node mkdir $dir
            if [ $? -eq 0 ]; then
                njob=$(expr $njob + 1)
            fi
            ;;
        "ln")
            echo "ssh $node ln -s /gpfs/home/joonho /home/joonho"
            #ssh $node ln -s /gpfs/home/joonho /home/joonho
            ;;
        *)
            echo "Job is not in case"
            ;;
        esac
    done
echo "$njob job succeeded"
echo "NB: only avaiable nodes are checked"
echo "Usage:: $0 vasp{default}"
echo "possible processes are ( ${jobs[@]} )"
