#!/bin/bash

hname=`hostname`
#echo $hname

software=qchem  # modify here 

qjob=dname
inf=input_file

qjobname=${1:-$qjob}
fin=${2:-$inf}
np=${3:-24}
if [ $hname == 'login' ]; then
    echo "qsub_server.py sge qchem -j $qjobname -i $fin -n $np"
    echo "want to run?"
    read var
    if [ $var ]; then
        qsub_server.py sge qchem -j $qjobname -i $fin -n $np
    fi
fi


