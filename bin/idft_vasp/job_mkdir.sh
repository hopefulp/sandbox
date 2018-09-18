#!/bin/bash

#dir1=cpo27			# ../vasp_relax in the "job_mkdir_sub"
#dir2=_${sys}_done		# _cell_done or _CO2_done
#dir3=_${sys}_d2
#dir_root=/qcfs/joonho/backup/vasp_relax/
dir_root=
dir_a=cpo27
#dir_b=_cell_d2
#dir_c=_cell_d2_r2
dir_b=_CO2_d2
dir_c=_CO2_d2_r2
#dir_b=_CO2_rgga
#dir_c=_CO2_rgga1

sys=CO2		# CO2 or cell or 1CO2
func=gga

for Me in `cat metal1.dir`
    do
	old_dir=$dir_root$dir_a$Me$dir_b				# dir = cpo27Ca_cell_done
	new_dir=$dir_a$Me$dir_c
	echo "OLD: $old_dir  New: $new_dir"
	./job_mkdir_sub.sh $old_dir $new_dir $Me  $1 $sys $func
    done


