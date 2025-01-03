#!/bin/bash

list_all=( CHG CHGCAR CONTCAR DOSCAR EIGENVAL IBZKPT OSZICAR  OUTCAR PCDAT REPORT vasprun.xml WAVECAR XDATCAR )
list_L1=( CHG CHGCAR CONTCAR DOSCAR EIGENVAL IBZKPT OSZICAR  OUTCAR PCDAT REPORT vasprun.xml WAVECAR XDATCAR ) 
if [ $# -lt 2 ]; then
    echo "Usage:: $0 dir_name option[All|all|L1|L2]"
    echo "     All => ${list_all[@]}"
    exit 0
fi

dir=$1
level=$2
if [ $level == "All" -o $level == "all" ]; then
    cd $dir
    for file in "${list_all[@]}"
	do
	    rm $file
  	done
    cd ..
else
    for file in $*
    	do
    	    find . -name $file |  xargs echo "rm"
      	done
fi


