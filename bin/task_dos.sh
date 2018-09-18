#!/bin/bash

bin=/qcfs/joonho/bin

##### old version
dosfile=dosall.pl
#dosfile=dosall_newf.pl

	# zz: pair: 57(O21) 43(O7); 58 (O22) 
	# ac: pair: 39(O) 52(O)
        # zz,ac: 9(Fe-up) 1(Fe-down) 5(Fe-up) 57, 43 (O-bridge) 167(C-ac) 373(H-Benzene) 
	       # 9-197, 300 (C-ac)
a=9
l=2

for dir in $@
    do
	echo $dir
	cd $dir
	    #$bin/$dosfile 
	    $bin/$dosfile a=$a 
	    $bin/$dosfile a=$a l=$l
	    $bin/$dosfile a=$a l=$l s=1
	    $bin/$dosfile a=$a l=$l s=-1
	    $bin/$dosfile a=$a l=$l m=0
	    $bin/$dosfile a=$a l=$l m=0 s=1
	    $bin/$dosfile a=$a l=$l m=0 s=-1
	    $bin/$dosfile a=$a l=$l m=1
	    $bin/$dosfile a=$a l=$l m=1 s=1
	    $bin/$dosfile a=$a l=$l m=1 s=-1
	    $bin/$dosfile a=$a l=$l m=2
	    $bin/$dosfile a=$a l=$l m=2 s=1
	    $bin/$dosfile a=$a l=$l m=2 s=-1
	    if [ $l -ge 2 ]; then
	    	$bin/$dosfile a=$a l=$l m=3
	    	$bin/$dosfile a=$a l=$l m=3 s=1
	    	$bin/$dosfile a=$a l=$l m=3 s=-1
	    	$bin/$dosfile a=$a l=$l m=4
	    	$bin/$dosfile a=$a l=$l m=4 s=1
	    	$bin/$dosfile a=$a l=$l m=4 s=-1
	    fi
	  cd ..
    done

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
