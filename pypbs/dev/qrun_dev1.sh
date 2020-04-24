#!/bin/bash

### Modify here
lsoftware=( qchem amp )
software=${lsoftware[1]}

###### Default values
### common
de_qjob=dname
de_inf=OUTCAR
de_np=16
### qchem
de_mem="2G"     # for default
### AMP
d_hl="10 10"
d_el=0.001
d_fl=0.01
d_data="0 1500 2000"
d_job=tr
hname=`hostname`

if [ $hname == "chi" ]; then
    fin=${1:-$de_inf}
    np=${2:-$de_np}
elif [ $hname == "login" ]; then
    qjobname=${1:-$de_qjob}
    fin=${2:-$de_inf}
    np=${3:-$de_np}
    mem=${4:-$de_mem}
    hl=${5:-$d_hl}
    el=${6:-$d_el}
    fl=${7:-$d_fl}
    data=${8:-$d_data}
    ampjob=${9:-$d_job}
fi
echo memory $mem G

if [ $# -eq 0 ]; then
    echo "Usage:: $0 software arg2 arg3 ...."
    echo "software:: [${lsoftware[@]}] is necessary"
    #echo "all args has default values: should add args one by one"
    echo "args in the order:: qname,fname,np,mem,hl,el,data-mine"
    echo "AMP: qrun.sh amp N1000 OUTCAR 16 5 '10 10' 0.001 '0 1500 1500 2000' "
    echo "Qchem: qrun.sh qchem PNP5 a.in 16 3"
    exit 1
fi


if [ $hname == 'login' ]; then
    if [ $software == 'qchem' ]; then
        echo "args: qname=$qjobname fname=$fin np=$np mem=$mem"
        if [ -z $mem ]; then
            echo "qsub_server.py qchem -qj $qjobname -i $fin -n $np "
        else
            echo "qsub_server.py qchem -qj $qjobname -i $fin -n $np -m $mem"
        fi
        read -p "want to run? [enter for yes/any]" var
        #if [ -z "$var" -o $var == 'y' ]; then
        if [ -z $var ]; then
            if [ -z $mem ]; then
                qsub_server.py qchem -qj $qjobname -i $fin -n $np
            else
                qsub_server.py qchem -qj $qjobname -i $fin -n $np -m $mem
            fi
        fi
    elif [ $software == 'amp' ]; then
        echo "args: qname=$qjobname fname=$fin np=$np mem=$mem hl=$hl el=$el data=$data"
        echo "qsub_server.py amp -qj $qjobname -i $fin -n $np -hl $hl -el $el -fl $fl -di $data -m $mem"
        read -p "want to run? [enter/no]" var
        if [ -z $var ]; then
            qsub_server.py amp -qj $qjobname -i $fin -n $np -hl $hl -el $el -fl $fl -di $data -m $mem
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
