#!/bin/bash
# .vas is for POSCAR read in VESTA
# usage: ./cp_find.sh pre_dir obj_dir filename
# e.g. : ./cp_find.sh . $qj/VASP CONTCAR
# 

dir=$1
bak=$2
file=$3
#ext=.vas

for f in `find $dir -name $file`
    do
	newf=`dirname $f | cut -c3-`	#cut and write after 3rd column
	fname=$newf-$file
	cp $f $bak/$fname$ext

    done

