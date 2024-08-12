#!/bin/bash

#if [ $# -lt 3 ]; then
#    echo "Usage:: $0 a l dirnames"
#    echo "        a= atom index"
#    echo "        l= angular quantum number"
#    echo "        input dirnames"
#fi

#a=$1
#l=$2

atoms=1	# all

for dir in $@
    do
	cd $dir
	    if [ $atoms == "all" ]; then
	  	dosall.pl DOSCAR
	  	dosall.pl a=1
	  	dosall.pl a=2 # l=$l 
	  	dosall.pl a=3
	  	dosall.pl a=4
	  	dosall.pl a=5 # l=$l 
	  	dosall.pl a=6 # l=$l 
	  	dosall.pl a=7
	  	dosall.pl a=8
	  	dosall.pl a=9 # l=$l 
	  	dosall.pl a=10 # l=$l 
	  	dosall.pl a=11
	  	dosall.pl a=12 # l=$l 
#	  	dosall.pl a=10 a=33:35 a=74:79
	  	dosall.pl a=49
	  	dosall.pl a=49:50 a=135  
	    else
		dosall.pl a=12 l=2
	        dosall.pl a=49
		dosall.pl DOSCAR
			
	    fi

#	    if [ $l == 1 ]; then
#	  	dosall.pl a=$a l=$l m=0
#	  	dosall.pl a=$a l=$l m=1
#	  	dosall.pl a=$a l=$l m=2
#	    fi
	  cd ..
    done
