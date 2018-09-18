#!/bin/bash
# at the first sub-directory there is no "POSCAR" file
# scan first directory then copy

old_dir=/qcfs/biduri/mof-74/hydrocarbon/from-exp

#if [ -n $1 ]; then
#    old_dir=$1
#fi
list="POSCAR CONTCAR INCAR KPOINTS .out"

new_dir=`pwd`

echo "new_dir=$new_dir"

for sub_dir in `ls $old_dir`
    do
	work_dir=$new_dir/$sub_dir
	mkdir $work_dir
	
	cd $old_dir/$sub_dir
	# the file name holds the present directory "./"
	find . -name "out" > $work_dir/test.txt
	
	cd $work_dir

	while IFS= read -r NAME 
    	    do
		NAME=${NAME/.\//}
  		cp -v "$old_dir/$sub_dir/$NAME" "${NAME//\//__}"
    	    done < "test.txt"
    done


