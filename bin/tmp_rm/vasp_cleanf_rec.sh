#!/bin/bash

echo " $0 removes VASP files from 1st upto 2nd sub-directory "

tag=tet
pd=${1:-$PWD}
tclear=${2:-$tag}
echo Working directory is $pd
cd $pd

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

isub=0
for dir in `ls -1`
do
    if [ -d $dir ]; then
	#echo $dir
	cd $dir
	    isub=`expr $isub + 1`
	    sdir[1]=$dir
	    if [ -e CONTCAR ]; then
		delete $isub $dir $tclear
	    fi
	    for sub in `ls -1`
	    do
		if [ -d $sub ]; then
		    cd $sub
			dirs=$dir/$sub
			sdir[2]=$dirs
		    	isub=`expr $isub + 1`
		    	delete $isub $dirs $tclear
			for sub2 in `ls -1`
			do
			    if [ -d $sub2 ]; then
				cd $sub2
				    dirs=$dir/$sub/$sub2
				    isub=`expr $isub + 1`
				    delete $isub $dirs $tclear
				 cd ..
				 isub=`expr $isub - 1`
			    fi
			done
			cd ..
		    isub=`expr $isub - 1`
	    	fi
	    done
	    cd ..
	isub=`expr $isub - 1`
    fi
done
   
if [ $tclear == "tet" ]; then
    echo "Usage:: $0 .[directory] option[all|W|CW]"
    echo "     :: \"all\" clears AECCAR* CHG* WAVE* DOSCAR OUTCAR PROCAR"
    echo "        \"W\" leaves WAVECAR"
    echo "        \"CW\" leaves WAVECAR CHGCAR"
fi   

