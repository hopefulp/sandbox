#!/bin/bash

dir=$1

#for dir in $@
#    do
	cd $dir
	    gnu_3f.sh  $3 $2  Fermi.dat
	    gnu_3f.sh  $4 $2  Fermi.dat
	    gnu_3f.sh  $5 $2  Fermi.dat
	    gnu_3f.sh  $6 $2  Fermi.dat
	    gnu_3f.sh  $7 $2  Fermi.dat
	    gnu_3f.sh  $8 $2  Fermi.dat
	    cd ..
#    done
