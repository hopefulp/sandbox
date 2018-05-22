#!/bin/bash

echo " $0 :: removes files after finding "

#pd=${1:-$PWD}

#echo Working directory is $pd
#cd $pd

function delete {
    if [ -e CONTCAR ]; then
         echo "Here is $1-th vasp dir $2"
         if [ $3 == "all" ]; then
             rm AECCAR* CHG* WAVE* DOSCAR OUTCAR PROCAR
         elif [ $3 == "W" ]; then
             rm AECCAR* CHG* DOSCAR OUTCAR PROCAR
         elif [ $ == "WC" ]; then
             rm AECCAR* CHG DOSCAR OUTCAR PROCAR
         fi
     fi
}

dir=$1
shift
fulldir=$PWD/$dir

echo "$dir, $fulldir"
for file  in $*
    do
	find $fulldir -name $file -exec rm -rf {} \;
	#find $fulldir -name $file 
	
    done
   

