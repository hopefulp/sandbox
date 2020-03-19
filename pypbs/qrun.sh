#!/bin/bash

### Default values
de_qjob=dname
de_inf=input_file
de_np=16
de_mem=""

hname=`hostname`
if [ $hname == "chi" ]; then
    narg=2
    fin=${1:-$de_inf}
    np=${2:-$de_np}
elif [ $hname == "login" ]; then
    narg=3
    qjobname=${1:-$de_qjob}
    fin=${2:-$de_inf}
    np=${3:-$de_np}
    mem=${4:-$de_mem}
fi
echo memory $mem G
if [ $# -lt $narg ]; then
    echo "server \"$hname\" requires $narg arguments: [qname] input_file np "
    exit 1
fi

### Modify here
software=qchem  

if [ $hname == 'login' ]; then
    if [ -z $mem ]; then
        echo "qsub_server.py sge qchem -j $qjobname -i $fin -n $np "
    else
        echo "qsub_server.py sge qchem -j $qjobname -i $fin -n $np -m $mem"
    fi
    read -p "want to run? [enter for yes/any]" var
    #if [ -z "$var" -o $var == 'y' ]; then
    if [ -z $var ]; then
        if [ -z $mem ]; then
            qsub_server.py sge qchem -j $qjobname -i $fin -n $np
        else
            qsub_server.py sge qchem -j $qjobname -i $fin -n $np -m $mem
        fi
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
