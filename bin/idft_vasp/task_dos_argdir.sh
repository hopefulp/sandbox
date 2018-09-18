#!/bin/bash

for dir in $@
    do
	cd $dir
	  ../dosall.pl 
	  ../dosall.pl a=1
	  ../dosall.pl a=1 l=0
	  ../dosall.pl a=1 l=1
	  ../dosall.pl a=3
	  ../dosall.pl a=3 l=0
	  ../dosall.pl a=3 l=1
	  cd ..
    done
