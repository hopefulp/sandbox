#!/bin/bash

for f in *
    do
	if [ -d $f ]; then
	    cd $f
	        if [ -e "AECCAR0" ]; then
		    echo "rm $f/AECCAR0"
		    echo "rm $f/AECCAR1"
		    echo "rm $f/AECCAR2"
		fi
		if [ -e "CHG" ]; then
		    echo "rm $f/CHG"
		    echo "rm $f/CHGCAR"
		fi
		cd ..
	fi
    done
		 
