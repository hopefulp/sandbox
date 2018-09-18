#!/bin/bash

if [ $# -ne 1 ]; then
    echo "stopped for ERROR: it needs one argument: number of nodes"
    exit 
fi

func=gga 		# gga or lda
old=relax2		# CO2sp, relax, CO2opt
w_type=co2_opt
for Me in `cat metal.dir`
    do
	old_dir=$func$Me$old			# dir = cpo27Ca_cell_done
	new_dir=$func$Me$w_type
	mkdir $new_dir
	echo "OLD: $old_dir  New: $new_dir"
	cp $old_dir/CONTCAR MOF.pos
	./insert_co2opt.pl MOF.pos CONTCAR.6co2.dirt	# this makes MOF.pos.add6co2
	cp MOF.pos.add6co2 	$new_dir/POSCAR
	cp $old_dir/POTCAR 	$new_dir
	cp $old_dir/KPOINTS	$new_dir
	~/bin/changeline_tmpl.pl incar.520.$func.6co2 $Me CO2
	cp t.incar $new_dir/INCAR
	~/bin/changeline.pl pbsfold-idft.csh $new_dir nodes=$1
    done


