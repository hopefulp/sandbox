#!/bin/bash
# at the first sub-directory there is no "POSCAR" file
# scan first directory then copy

#old_dir=/qcfs/biduri/mof-74/hydrocarbon/from-exp
old_dir=/qcfs/joonho/backup/PBE_MOF_CO2/PBED2

if [ -n "$1" ]; then
    old_dir=$1
fi

new_dir=`pwd`

echo "new_dir=$new_dir"

for sub_dir in `ls $old_dir`
    do
	work_dir=$new_dir/$sub_dir
	mkdir $work_dir
	
	# the file name holds the present directory "./"
	for file in POSCAR CONTCAR INCAR KPOINTS out
	    do
		cd $old_dir/$sub_dir
		find . -name "*${file}*" > $work_dir/test.txt
	
		cd $work_dir
		### to read text file line-by-line
		while IFS= read -r NAME 
    	    	    do
			NAME=${NAME/.\//}
  			cp -v "$old_dir/$sub_dir/$NAME" "${NAME//\//__}"
    	    	    done < "test.txt"
	    done
    done


