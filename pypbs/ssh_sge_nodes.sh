#!/bin/bash

nodes=$(qstat -f | sed -e '/---/d' -e '/adus/d' | awk '/@/ { print $1 }' | awk -F@ '{print $2}'  )
#echo $nodes

jobs=( mkdir vasp qchem )
job=${1:-"vasp"}
#echo $job

for node in $nodes; do
    echo $node
    #case
    ssh $node ps aux | grep $job
    done

echo "NB: only avaiable nodes are checked"
echo "Usage:: $0 vasp{default}"
echo "possible processes are ( ${jobs[@]} )"
