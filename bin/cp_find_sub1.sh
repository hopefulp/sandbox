#!/bin/bash
# .vas is for POSCAR read in VESTA
# usage: ./cp_find_sub1.sh pre_dir obj_dir filename
# e.g. : ./cp_find_sub1.sh . $qj/VASP_dat CONTCAR
# 

dir=$1
bak=$2
file=$3
#ext=.vas

if [ $# -lt 3 ]; then
    echo "Usage:: ./cp_find_sub1.sh pre_dir obj_dir filename"
    echo "e.g. :: ./cp_find_sub1.sh . $qj/VASP_dat CONTCAR"
    echo "note :: the present directory should be ./"
    exit
fi

for f in `find $dir -name $file`
    do
	newf=`dirname $f | cut -c3-`	#cut and write after 3rd column
	fname=$newf-$file
	cp $f $bak/$fname$ext

    done

