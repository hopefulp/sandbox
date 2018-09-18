#!/bin/bash

dir=$1

#for dir in $@
#    do
	cd $dir
	    gnu_2f.sh  $3 $2  $dir
	    gnu_2f.sh  $4 $2  $dir
	    gnu_2f.sh  $5 $2  $dir
	    gnu_2f.sh  $6 $2  $dir
	    gnu_2f.sh  $7 $2  $dir
	    gnu_2f.sh  $8 $2  $dir
	    cd ..
#    done
