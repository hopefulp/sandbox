#!/bin/bash

sys=cell			# cell or CO2
dir1=cpo27
dir2=_${sys}_done		# _cell_done or _CO2_done

func=gga	 		# gga or lda
w_type=relax			# CO2sp, relax, CO2opt
w2_type=relax2

for Me in `cat metal.dir`
    do
	old_dir=$func$Me$w_type 			# dir = cpo27Ca_cell_done
	new_dir=$func$Me$w2_type
	echo "OLD: $old_dir  New: $new_dir"
	./job_mkdir_cont_sub.sh $old_dir $new_dir 
	~/bin/changeline.pl pbsfold-idft.csh $new_dir nodes=$1
    done


