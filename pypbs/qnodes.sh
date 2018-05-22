#!/bin/bash

echo " to make known_hosts key for login to nodes without passwd"
echo " first empty .ssh/known_hosts, then, run this script "

host=$(hostname)
id=$(whoami)

if [ $host == "psi" ]; then
    fnode=44
elif [ $host == "rho" ]; then
    fnode=20
elif [ $host == "idft" ]; then
    fnode=22
elif [ $host == "kdft" ]; then
    fnode=26
fi


for i in $(seq 1  $fnode)
  do
    if [ $i -lt 10 ]; then
	node=${host}0$i
    else
	node=${host}$i
    fi

    echo $node
    #ssh $node ls /scratch
    ssh $node ls /scratch/$id
    if [ $? -ne 0 ]; then
	    ssh $node mkdir /scratch/$id
	    echo Made /scratch/$id
    fi
    ssh $node ls /scratch/$id/QCLOC
    if [ $? -ne 0 ]; then
    	ssh $node mkdir /scratch/$id/QCLOC
	    echo Made /scratch/$id/QCLOC
    else
	    ssh $node rm -rf /scratch/$id/QCLOC/qchem*
    fi
    ssh $node ls /scratch/$id/QCSCR
    if [ $? -ne 0 ]; then
    	ssh $node mkdir /scratch/$id/QCSCR
	    echo Made /scratch/$id/QCSCR
    else
	    ssh $node rm -rf /scratch/$id/QCLOC/qchem*
    fi
  done
