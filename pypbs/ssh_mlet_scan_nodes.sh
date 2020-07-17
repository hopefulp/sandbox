#!/bin/bash

nodes=$(qstat -f | sed -e '/---/d' -e '/adus/d' | awk '/@/ { print $1 }' | awk -F@ '{print $2}'  )
#echo $nodes

jobs=( ls rm mkdir qchem ln ps)

### Use default at the moment
job_default='ps'
pname='python'
if [ $# -eq 0 ]; then
    echo "Usage:: $0 job=[ ${jobs[@]} ] arg2"
    echo "Use default: job '$job_default' process_name '$pname'"
    #exit 1
fi


job=${1:-$job_default}


#echo $job
if [ $job == "ls" -o $job == "mkdir" -o $job == "rm" -o $job == "ps" ]; then
    arg2=${2:-$pname}
    #if [ -z $arg2 ]; then
    #    echo "input 2nd argument for $job as \$2"
    #    exit 2
    #fi
fi

njob=0

for node in $nodes; do
    echo $node
    case $job in
        "ps")
            ssh $node ps aux | grep $arg2
            ;;
        "ls")
            ssh $node ls $arg2
            if [ $? -eq 0 ]; then
                njob=$(expr $njob + 1)
            fi
            ;;
        "rm")
            for f in $(ssh $node ls $arg2); do
                ssh $node rm -r $arg2/$f
                done
            ;;
        "mkdir")
            ssh $node mkdir $arg
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
