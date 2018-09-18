#!/bin/bash

#dir=$1
dir_a=cpo27
#dir_b=_cell_d2_r2
dir_b=_CO2_d2_r2
#for Me in `cat metal.dir`
#    do
#	dir=${dir_a}${Me}$dir_b
#	rm -r $dir
#    	cd $dir
    	for x in `ls -1 PARCHG*  `
  	    do
     	    	echo  $x
      	    	mv $x $x.vas
  	    done

#    	cd ..
#	grep F= $dir.out
#    done
