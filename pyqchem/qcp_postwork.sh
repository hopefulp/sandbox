#!/bin/bash

job_list=( clean cleanall cleanscr getmol xyzmol )

if [ $# -lt 1 ]; then
    echo "Usage::$0 job dir's"
    echo "  job = ${job_list[*]}"
    echo "subroutine::getmol: qc_m_getmol_all.sh -> qgeo_mo.pl"
    exit
fi

#echo $HOST
#exit

Job=$1

host=$(hostname)
id=$(whoami)

if [ $host == "psi" ]; then
    fnode=44
elif [ $host == "rho" ]; then
    fnode=20
elif [ $host == "kdft" ]; then
    fnode=26
fi

case $Job in
    "${job_list[0]}")
        rm *.e* *.o[0-9]* PI[1-9]*
        ;;
    "${job_list[1]}")
	    rm *.e* *.o[0-9]* *in *mol *xyz [0-9]*.$HOST* *.fchk
	    ;;
    "${job_list[2]}")
	    for i in $(seq 1 $fnode)
            do
                if [ $i -lt 10 ]; then
                    node=${host}0$i
                else
                    node=${host}$i
                fi

                echo $node
                ssh $node ls /scratch/$id
                ssh $node ls /scratch/$id/QCSCR
                ssh $node rm -rf /scratch/$id/QCSCR/*
                done
        ;;
    "${job_list[3]}")
        qcp_getmol_all.sh
        ;;
    "${job_list[4]}")
        for x in `ls *.mol *.xyz`
            do
                xyz22mol.pl $x
            done
        ;;
    *)
        echo "there is no $Job option in $job_list"
        exit
        ;;
esac
