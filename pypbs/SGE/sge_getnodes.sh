#!/bin/bash

nodes=( node23 node22 node21 node20 node16 node15 node11 )

i=0
d_nnode=1
nnode=${1:-$d_nnode}

for node in ${nodes[@]}
    do
        num=${node:4:2}
        echo node$num
        Lmatch=0
        for run_node in `qstat | awk '/@/ {print $8}' | cut -d@ -f2`
            do
                if [ $run_node == $node ]; then
                    Lmatch=1
                    break
                fi
            done
        if [ $Lmatch -eq 0 ]; then
            #echo qsub -N amplong$num -pe numa 36 -q skylake@$node $SB/pypbs/sge_sleep.csh
            st="qsub -N amplong$num -pe numa 36 -q skylake@$node $SB/pypbs/sge_sleep.csh"
            read -p  "$st | Will you run ? [enter/no]" quest
            if [ -z $quest ]; then
                #qsub -N amplong$num -pe numa 36 -q skylake@$node $SB/pypbs/sge_sleep.csh
                $st
            fi
            i=`expr $i + 1`
            if [ $i -eq $nnode ]; then
                break
            fi
        fi
    done

