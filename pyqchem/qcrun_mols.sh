#!/bin/bash

#rem=eda_pbed

#if [ $# -ne 0 ]; then
#  rem=$1
#fi

pbs=pbs_qchem1_$HOST.sh
bin=/qcfs/joonho/bin

if [ ! -e $pbs ]; then
    echo "There is no pbs file of $pbs"
    exit
fi

case $HOST in
    "kdft")
	ppn=12
    	;;
    "psi")
	ppn=16
	;;
    *)
	echo "Error:: There is no host in the host list"
	exit
	;;
esac

for x in `ls *.mol`
  do
    fmol=$( basename $x | cut -d'.' -f1 )
    #cat rem.eda_pbed $f_pre.mol > $f_pre.in
    $bin/changeline_pbs_QChem.pl $pbs $fmol  $ppn
  done

