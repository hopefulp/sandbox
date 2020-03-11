#!/bin/bash

hname=`hostname`
#echo $hname

software=qchem  # modify here 

qjob=dname
inf=input_file

fin=${1:-$inf}
np=${2:-24}
qjobname=${3:-$qjob}

if [ $hname == 'login' ]; then
    if [ $# -ne 3 ]; then
        echo "$hname requires 3 arguments: input_file np qname"
        exit 1
    fi
    echo "qsub_server.py sge qchem -j $qjobname -i $fin -n $np"
    echo "want to run?"
    read var
    if [ $var ]; then
        qsub_server.py sge qchem -j $qjobname -i $fin -n $np
    fi
elif [ $hname == 'chi' ]; then
    if [ $# -ne 2 ]; then
        echo "$hname requires 2 arguments: input_file np "
        exit 1
    fi
    echo "qsub_server.py chi qchem -i $fin -n $np"
    echo "want to run?"
    read var
    if [ $var ]; then
        qsub_server.py chi qchem -i $fin -n $np
    fi
fi
