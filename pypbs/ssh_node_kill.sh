#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage:: $0 pname(to be killed) node_list"
    exit 1
fi

pname=$1
shift

for node in $@; do
    echo $node
    echo "ssh $node ps aux | grep $pname | awk '{print \$2}'"
    ssh $node ps aux | grep $pname | awk '{print $2}'
    for proc in $(ssh $node ps aux | grep $pname | awk '{print $2}'); do
        echo "ssh $node kill $proc"
        ssh $node kill $proc
        done
    done

