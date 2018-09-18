#!/bin/bash

if [ $# -ne 1 ]; then
    echo "stopped for ERROR: it needs one argument: number of nodes"
    exit 
fi

#func=gga 		# gga or lda
#old=relax2		# CO2sp, relax, CO2opt
#w_type=co2_opt

### adding 1 CO2 in the long distance from metal
old=../vasp_relax/cpo27
old2=_cell_done
new=cpo27_1co2_

### removing 5 CO2 from MOF-CO2 to gain MOF-1CO2 complex
old=../vasp_relax/cpo27
old2=_CO2_done
new=cpo27_1co2_ads_

for Me in `cat metal1.dir`
    do
	old_dir=$old$Me$old2			# dir = cpo27Ca_cell_done
	new_dir=$new$Me
	mkdir $new_dir
	echo "OLD: $old_dir  New: $new_dir"
	cp $old_dir/CONTCAR MOF.pos
#	./insert_co2.pl MOF.pos CONTCAR.1co2.dirt	# this makes CONTCAR.new
	./ext_rm_nco2a.pl MOF.pos 			# this makes MOF.pos.abs1co2
	cp MOF.pos.ads1co2	$new_dir/POSCAR
	cp $old_dir/POTCAR 	$new_dir
	cp $old_dir/KPOINTS	$new_dir
	~/bin/changeline_tmpl.pl incar.sp.spin.ini $Me 1CO2
	cp t.incar $new_dir/INCAR
	~/bin/changeline.pl pbsfold-idft.csh $new_dir nodes=$1
    done


